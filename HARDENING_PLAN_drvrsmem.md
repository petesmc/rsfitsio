# Hardening Plan — `shared_mem` driver (`src/drvrsmem.rs`)

Status: proposed. Run in a copy of the repo.

Scope: the SysV-IPC shared-memory FITS driver, gated
`#[cfg(all(feature = "shared_mem", not(target_os = "windows")))]`. It is a
line-by-line port of cfitsio's `drvrsmem.c` and preserves C idioms (global
`static mut` tables, raw `memcpy`, `sscanf`, `union semun`). The goal is to
remove language-level undefined behavior and the ported logic bugs that have a
memory-safety dimension, **without changing observable behavior or the C ABI**.

## Ground rules

- Keep the `#[no_mangle] extern "C"` symbol names, signatures, and return codes
  identical. Callers (including `src/bin/smem`) must still link and behave the same.
- Preserve the fcntl lock-file + SysV semaphore protocol exactly — that is the
  cross-process contract and other processes / the C library rely on it.
- After each phase: `cargo build --features shared_mem`,
  `cargo test --features shared_mem -- --test-threads=1` (segments are
  process-global with fixed names `h0..h14`; parallel tests interfere), and
  `cargo clippy --features shared_mem`.
- Do not attempt to run the shared_mem tests under Miri — they call real
  `shmget`/`semget`/`fcntl` and need `-Zmiri-disable-isolation` and a live kernel;
  treat Miri as out of scope for this module.
- Land phases as separate commits so each is reviewable and revertible.

## Baseline to capture before starting

- [ ] `cargo build --features shared_mem` clean.
- [ ] `cargo test --features shared_mem -- --test-threads=1` — record pass/fail
      list; this is the regression oracle for every phase.
- [ ] Note that `src/drvrsmem.rs:1` has `#![allow(static_mut_refs)]`; the end goal
      is to be able to remove that allow (or scope it down) once Phase 1 lands.

---

## Phase 0 — Make the tests a merge gate before touching anything (do first)

**Why.** The 20 existing tests in `src/drvrsmem.rs` (`:1987-2529`) cover the
*serial* path well and pass 20/20 with `--test-threads=1`, but:
- they never run in CI — `.github/workflows/build.yml` runs `cargo test` with no
  `--features shared_mem`, and the module is feature-gated, so 0 shmem tests run;
- run in parallel they collide on the fixed global names `h0..h14` and the
  process-global driver state — confirmed: failures plus a panic at
  `shared_free` (`:1308`, `SHARED_LT[idx].p.unwrap()`), i.e. exactly the Tier-1
  race / Tier-2 panic this plan removes;
- nothing exercises the actual cross-process sharing (see the new
  "Cross-process test" section after Phase 2).

We want the serial suite gating merges *before* refactoring, so every later
phase has a green oracle. **Do NOT** globalize `--test-threads=1` in CI — that
would serialize and potentially mask genuine concurrency issues in the rest of
the suite. Scope the serialization to the shmem tests only.

**Change — serialize the shmem tests among themselves (dependency-free).**
- [ ] Add a module-local gate in the `drvrsmem` test module:
      ```rust
      // Serializes the shmem tests w.r.t. EACH OTHER only; the rest of the suite
      // still runs in parallel. Recover on poison so one failure doesn't cascade.
      static SHMEM_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
      fn shmem_guard() -> std::sync::MutexGuard<'static, ()> {
          SHMEM_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner())
      }
      ```
- [ ] Make each `#[test]` take the guard as its first line:
      `let _g = shmem_guard();` — hold it for the whole test body so no two shmem
      tests overlap. (Equivalent to the `serial_test` crate's `#[serial]`, without
      adding a dependency — this repo minimizes deps.)
- [ ] Now the tests pass under a normal parallel `cargo test --features
      shared_mem --lib` (they serialize via the mutex; everything else stays
      parallel). Verify locally both ways.

**Change — wire into CI as a separate job/step.**
- [ ] In `.github/workflows/build.yml`, on the Linux and macOS jobs (NOT Windows —
      the module is `not(target_os = "windows")`), add a step after the existing
      test step:
      ```yaml
      - name: Run shared_mem driver tests
        run: cargo test --features shared_mem --lib drvrsmem
      ```
      No `--test-threads=1` needed thanks to the mutex gate. Keep the existing
      feature-off `cargo test` step exactly as-is so the main suite stays parallel
      and unchanged.
- [ ] Add a teardown safety net so leaked SysV segments don't fail re-runs on a
      persistent runner (see the cleanup note in the Cross-process section):
      `ipcs -m` before/after for visibility, and `- if: always()` cleanup that
      removes this suite's segments (a unique `SHMEM_LIB_KEYBASE`, below, makes
      that targetable).

**Isolation lever (use throughout Phase 0 and the cross-process tests).** The
driver derives its IPC keys and lock-file name
(`/tmp/.shmem-lockfile.<keybase>.<maxseg>`, `:415`) from `SHMEM_LIB_KEYBASE` /
`SHMEM_LIB_MAXSEG` (`:373`). Setting a unique `SHMEM_LIB_KEYBASE` per CI job (or
per test binary) gives that run a private segment namespace, so concurrent jobs
on a shared runner — and any unrelated shmem user on the box — cannot collide.

**Verify.** `cargo test --features shared_mem --lib drvrsmem` green in parallel;
main suite unaffected; CI job added and passing on Linux + macOS.

---

## Phase 1 — Global-table representation: kill the `&mut [T]` aliasing (Tier 1)

**Problem.** `static mut SHARED_GT: &mut [SHARED_GTAB]`
(`src/drvrsmem.rs:168`), populated via `slice::from_raw_parts_mut` at ~`:482`
and `:494`, is a Rust exclusive reference into a SysV segment that **other
processes concurrently mutate by design**. `&mut` asserts unique, non-aliased
access; that assumption is false here, so every access is UB at the model level
regardless of threads. `SHARED_GT_PTR` (`:169`) already holds the raw pointer.

**Why atomics, not volatile.** The obvious fix is "raw pointer +
`read_volatile`/`write_volatile`," but volatile only means "don't optimize the
access away" — it carries **no** concurrency semantics. The fields in
`SHARED_GTAB` are genuinely modified by other processes, so the correct model is
**atomics**: `AtomicI32` load/store with `Ordering::Relaxed`. `AtomicI32` is
lock-free on every platform this builds for, has the same size/layout as `c_int`
(`#[repr(C)]`-compatible), and gives tear-free cross-process integer access.
The real mutual exclusion still comes from the fcntl lock-file + semaphores;
the atomics just make the individual field reads/writes well-defined instead of
UB. This is a change from the earlier "volatile" sketch — atomics are strictly
better here.

Caveat to record in a comment: the C library accesses these same fields as plain
`int`. Rust-atomic ↔ C-plain-int on one address is, pedantically, a data race in
the abstract machine; in practice it is fine because the fcntl/semaphore locks
serialize real access and `Relaxed` only requires tear-free integer load/store.

> ### BUG-1 — global-table cleanup is dead code (port regression, leaks a SysV segment)
>
> **Confirmed, with reproducer. Not present in cfitsio.** Fixed as part of the
> `SHARED_GT_PTR.is_null()` rewrite above, but tracked explicitly because it is a
> real resource leak that diverges from C behavior.
>
> `shared_cleanup` guards the entire global-table detach/delete block with
> ([drvrsmem.rs:270](src/drvrsmem.rs#L270)):
> ```rust
> if !SHARED_GT.len() == 0 {   // parses as (!SHARED_GT.len()) == 0
> ```
> `.len()` is `usize` and `!` is **bitwise** NOT (binds tighter than `==`), so:
> - attached → `len==16` → `!16usize` is huge → `!= 0` → **false**
> - empty    → `len==0`  → `!0usize == usize::MAX` → `!= 0` → **false**
>
> The condition is **always false**, so the block never runs: the global-table
> segment is never `shmdt`'d, never `IPC_RMID`'d, and the whole-file lock is never
> released. The C original ([cfitsio drvrsmem.c:112](drvrsmem.c#L112)) guards this
> with `if (NULL != shared_gt)` and, when this is the last process with no
> segments left, actually `shmctl(shared_gt_h, IPC_RMID, 0)`s the table.
>
> **Impact:** every process that attaches the global table under a given
> `SHMEM_LIB_KEYBASE` leaks one ~448-byte SysV segment (16 × `SHARED_GTAB`) at
> exit; over time this exhausts the kernel `shmmni` limit and causes key
> collisions on re-runs. Observed directly by the cross-process test — its
> teardown needs an `ipcrm -M <keybase>` backstop *today* purely to paper over
> this bug; once BUG-1 is fixed that backstop should become unnecessary (keep it
> as belt-and-suspenders, but add an assertion that nothing leaked).
>
> **Fix:** `if !SHARED_GT_PTR.is_null()` (the faithful translation of
> `NULL != shared_gt`). Add a regression test: run a full create→delete cycle in
> a child process under a private keybase, then assert (via `ipcs`/`shmget`
> probe) that the keybase's global-table segment is gone.

**Change.**
- [ ] Stop storing a slice. Keep only `SHARED_GT_PTR: *mut SHARED_GTAB` plus a
      length (reuse `SHARED_MAXSEG`). Remove `static mut SHARED_GT: &mut [...]`.
- [ ] Change the integer fields of `SHARED_GTAB` (`:101`) — `sem`, `semkey`,
      `key`, `handle`, `size`, `nprocdebug` — from `c_int` to
      `std::sync::atomic::AtomicI32`. Leave `attr: c_char` as-is (see below).
      `AtomicI32` is `#[repr(C)]`/`#[repr(transparent)]` over `i32`, so segment
      layout and the C ABI are unchanged. `SHARED_GTAB` must NOT derive
      `Copy`/`Clone` after this (atomics aren't `Copy`) — access fields in place.
- [ ] Introduce a field-access helper that indexes the raw pointer and never
      forms a `&mut` into the segment. Returning a `&` to the atomic is fine
      (shared refs to atomics that others mutate are sound; `&mut` is not — never
      hand out `&mut`). Take the base pointer as a parameter so it composes with
      the Phase 2 `SharedState` (which owns the pointer):
      `unsafe fn gt<'a>(base: *mut SHARED_GTAB, idx: usize) -> &'a SHARED_GTAB { &*base.add(idx) }`
      then `gt(base, idx).key.load(Relaxed)` / `.store(v, Relaxed)`. Before
      Phase 2 lands, `base` is `SHARED_GT_PTR`; after, it is `st.gt_ptr`.
- [ ] `attr: c_char` is the odd one out (byte, not `i32`, so no `AtomicI8`
      needed unless you want uniformity). It IS written cross-process
      (`shared_set_attr` `:1437`, `shared_malloc` `:980`). Simplest sound option:
      make it `AtomicU8` too. If left as plain `c_char`, document that it is only
      ever written under the fcntl write-lock and read under a lock, so no torn
      access — but prefer `AtomicU8` for consistency.
- [ ] Replace every `SHARED_GT[i].field` read/write (~40 sites; grep `SHARED_GT`)
      with `gt(i).field.load/store(..., Relaxed)`.
- [ ] Replace `SHARED_GT.is_empty()` / `SHARED_GT.len()` init checks with
      `SHARED_GT_PTR.is_null()`. **This fixes BUG-1 (see callout below)** — the
      `!SHARED_GT.len() == 0` guard at `:270` must become `!SHARED_GT_PTR.is_null()`
      ("table is attached"), matching C's `if (NULL != shared_gt)`. Add a
      regression test that asserts no global-table segment leaks after cleanup.
- [ ] Where the code copies a whole `SHARED_GTAB` (if any), replace with
      per-field atomic loads into a plain local struct; you can't `read_volatile`
      a struct of atomics.

**Layout note.** `SHARED_GTAB` (`:101`) mixes `c_int` fields with a trailing
`c_char`. Swapping `c_int`→`AtomicI32` and `c_char`→`AtomicU8` preserves size,
alignment, and offsets (atomics are `repr(transparent)` over their integer), so
the segment stays byte-compatible with what the C side wrote. Do **not** change
field order or `repr`.

**Verify.** Build + single-threaded tests green. Behavior byte-identical. If you
have the C `cfitsio`/`smem` build handy (see the ground-truth tooling note),
create a segment with one and inspect/delete it with the other to confirm the
on-segment layout still interoperates.

---

## Phase 2 — One in-process lock over the local state (Tier 1)

**Problem.** `SHARED_LT` (`Vec` in `static mut`, `:170`), `SHARED_FD`,
`SHARED_INIT_CALLED`, `SHARED_GT_PTR`, and `static mut COUNTER` in
`shared_get_hash` (`:841`) are read/written with no intra-process
synchronization. Two threads in any `smem_*`/`shared_*` entry point race = UB.
The `unsafe impl Send for SHARED_LTAB` (`:120`) and the "run single-threaded"
test comment acknowledge but do not contain this.

**Change (pick one, prefer A).**
- [ ] **A — coarse global mutex.** Wrap all process-local state in a single
      `SharedState` (drafted below) behind a `static STATE: Mutex<SharedState>`.
      Each public entry point locks once at the top. This is the smallest change
      that makes the intra-process race go away. Keep the fcntl/semaphore calls
      exactly where they are — they handle the *cross*-process side; the mutex is
      purely intra-process.
- [ ] **B — if a full struct refactor is too invasive**, at minimum move
      `COUNTER` (`:841`) and the init/attach path under a lock, and document the
      remaining single-thread requirement loudly. (Weaker; only do this if A is
      blocked.)

**Drafted `SharedState` (option A).** Note `gt_ptr` stays a raw pointer (Phase 1
made the *fields* atomic; the pointer itself is process-local state and belongs
in here). The mutex protects the process-local bookkeeping; the atomics + fcntl
+ semaphores protect the cross-process segment.

```rust
use std::sync::Mutex;

/// All process-local shared-memory-driver state, guarded by one mutex.
/// Cross-process state lives in the SysV segments themselves (Phase 1 atomics)
/// and is serialized by the fcntl lock-file + semaphores, NOT by this mutex.
struct SharedState {
    fd: c_int,                 // was SHARED_FD  — global lock-file descriptor
    gt_ptr: *mut SHARED_GTAB,  // was SHARED_GT_PTR — attached global table
    gt_h: c_int,               // was SHARED_GT_H — global-table shmem handle
    lt: Vec<SHARED_LTAB>,      // was SHARED_LT  — per-segment local table
    kbase: c_int,              // was SHARED_KBASE
    maxseg: c_int,             // was SHARED_MAXSEG
    range: c_int,              // was SHARED_RANGE
    init_called: bool,         // was SHARED_INIT_CALLED
    debug: bool,               // was SHARED_DEBUG
    create_mode: c_int,        // was SHARED_CREATE_MODE
    hash_counter: c_int,       // was the `static mut COUNTER` in shared_get_hash
}

// SAFETY: the raw pointers are only dereferenced while the mutex is held, and
// the memory they point at is synchronized cross-process via fcntl + atomics.
unsafe impl Send for SharedState {}

impl SharedState {
    const fn new() -> Self {
        SharedState {
            fd: SHARED_INVALID,
            gt_ptr: ptr::null_mut(),
            gt_h: SHARED_INVALID,
            lt: Vec::new(),
            kbase: 0,
            maxseg: 0,
            range: 0,
            init_called: false,
            debug: false,
            create_mode: 0o666,
            hash_counter: 0,
        }
    }
}

static STATE: Mutex<SharedState> = Mutex::new(SharedState::new());
```

(`Mutex::new` is `const` since Rust 1.63, so no `OnceLock`/lazy init needed —
this satisfies the 1.88 MSRV.)

**Entry-point shape (handles the re-entrancy trap).** Several functions call
each other internally — `shared_malloc` → `shared_init`; `shared_free` →
`shared_validate` → `shared_mux`; `smem_remove` → `smem_open`; every op →
`shared_check_locked_index` → `shared_init`. A single non-reentrant `Mutex`
would self-deadlock. Resolve by splitting each function into a public locking
wrapper and a `_locked` core that assumes the guard is already held and only
ever calls other `_locked` cores:

```rust
#[cfg_attr(not(test), unsafe(no_mangle))]
pub unsafe extern "C" fn shared_malloc(size: c_long, mode: c_int, newhandle: c_int) -> c_int {
    // Recover the guard rather than .unwrap()/.expect() it: lock poisoning
    // genuinely CAN happen (another thread panicked mid-op), and a panic here
    // would unwind across the extern "C" boundary (UB); C never panics on a lock.
    let mut st = STATE.lock().unwrap_or_else(|e| e.into_inner());
    shared_malloc_locked(&mut st, size, mode, newhandle)
}

// Internal: assumes `st` is locked; calls only *_locked helpers. Never locks.
fn shared_malloc_locked(st: &mut SharedState, size: c_long, mode: c_int, newhandle: c_int) -> c_int {
    /* ... former body, with SHARED_GT[i]/SHARED_LT[i]/SHARED_* replaced by
       st.gt_ptr-based atomic accessors and st.lt[i] / st.field ... */
}
```

Mechanical porting rules for each body:
- `SHARED_LT[i]`            → `st.lt[i]`
- `SHARED_MAXSEG` etc.      → `st.maxseg` etc.
- `SHARED_GT[i].field`      → `gt(st.gt_ptr, i).field.load/store(Relaxed)` (Phase 1)
- internal call `shared_x(..)` → `shared_x_locked(st, ..)`
- the `static mut COUNTER` in `shared_get_hash` → `st.hash_counter`

**Do this before wiring the mutex:** map the internal call graph so you know
which functions need a `_locked` sibling. Grep for internal callers:
`grep -nE 'shared_(init|mux|demux|map|validate|attach|free|lock|unlock|destroy_entry|clear_entry|check_locked_index|get_free_entry|realloc)\b' src/drvrsmem.rs`
and `grep -n 'smem_open' src/drvrsmem.rs`. Everything reachable from a public
entry point needs the `_locked` form.

**Poisoning (recover the guard, don't panic):** if a thread panics while holding
the lock, `STATE.lock()` returns `Err`. Poisoning genuinely *can* occur, so this
is not an "impossible" case — **always recover the guard** with
`.lock().unwrap_or_else(|e| e.into_inner())` and continue. Never `.unwrap()`,
`.expect()`, or `panic!` on the lock: matching C behavior means returning
`SHARED_*` codes, and panicking out of an `extern "C"` fn is UB (the very thing
Phase 4 removes). Since Phase 4 wraps every body in `catch_unwind`, poisoning
should be rare regardless. Apply this identical pattern at every lock site (all
public entry points and `shared_cleanup`).

**When `.expect()` IS appropriate (elsewhere):** the rule above is specific to
the lock, whose failure is recoverable. For conditions that are *genuinely
impossible* given an established invariant — not merely "shouldn't happen" —
prefer `.expect("why this cannot fail")` over silent handling: it documents the
invariant and, if it ever fires, flags a real logic bug rather than masking it.
The test is strict: use `.expect()` only when a failure would mean the code's own
invariants are already broken (so C would be equally broken / UB there too), NOT
for anything an untrusted caller, a racy cross-process write, or the environment
can trigger — those must return a `SHARED_*` code. See Phase 4 for applying this
to the existing `try_into().unwrap()` sites.

**`atexit(shared_cleanup)` interaction:** `shared_cleanup` will also need to lock
`STATE`. It runs at process exit on the exiting thread; taking the mutex there is
fine unless the process aborted mid-operation while holding it (then it's
poisoned — recover the guard as above). Keep the Phase 4 `catch_unwind` around
its body regardless.

**Internal RAII helper (folded in from the reviewed snippet).** The other
agent's `SystemVSharedMemory` idea is worth keeping as a *private* attach/detach
guard used inside `shared_malloc`/`shared_attach`/`shared_map`, where today the
"shmat → on any error: shmdt + shmctl(RMID) + continue" ladders (`:926-958`,
`:1224-1231`) are hand-unwound and easy to leak. It is NOT a replacement for
`shared_free`/`shared_destroy_entry` — the *destroy* decision depends on the
semaphore process-count and the `PERSIST` attr and must stay there. Keep it
minimal and detach-only:

```rust
/// RAII guard for an attached SysV segment: detaches on drop. Does NOT destroy
/// (IPC_RMID) — destruction is decided by shared_free via the semaphore count.
struct AttachedSeg(*mut BLKHEAD);
impl AttachedSeg {
    /// Returns None on shmat failure (compared against the explicit -1 sentinel,
    /// see Phase 5b — do not rely on `SHARED_INVALID as *mut`).
    unsafe fn attach(handle: c_int) -> Option<Self> {
        let p = libc::shmat(handle, ptr::null(), 0) as *mut BLKHEAD;
        if p == (-1isize as *mut BLKHEAD) { None } else { Some(AttachedSeg(p)) }
    }
    fn as_ptr(&self) -> *mut BLKHEAD { self.0 }
    /// Give up ownership without detaching (segment stays attached in the table).
    fn into_raw(self) -> *mut BLKHEAD { let p = self.0; std::mem::forget(self); p }
}
impl Drop for AttachedSeg {
    fn drop(&mut self) { unsafe { libc::shmdt(self.0 as *const c_void); } }
}
```

Use `into_raw()` on the success path (when the pointer is stored into
`st.lt[idx].p`), and let `drop` clean up on the early-return error paths. This is
optional polish for Phase 2 — land the mutex first, then refactor the ladders if
time allows.

**Watch for:** re-entrancy (covered above via `_locked` split) and holding the
mutex across the blocking `fcntl(F_SETLKW)` calls in `shared_mux`. Holding the
process-local mutex while blocked on a cross-process file lock serializes all
threads behind one segment's waiter — acceptable for correctness, but note it;
if it becomes a throughput problem, that's a later refinement, not a Phase 2
concern.

**Verify.** Build + serial suite green (Phase 0 gate). Then add one new
in-process concurrency test that spawns 2–3 threads each
creating/writing/reading/deleting *distinct* handles concurrently, asserting no
panic and correct data. Have the test own its threads internally and take the
Phase 0 `shmem_guard()` so it doesn't overlap the other shmem tests.

**Scope note (important):** this thread test validates the **Phase 2 mutex**
only — it does NOT validate the driver's fcntl/semaphore locking protocol.
POSIX `fcntl` record locks (`shared_mux`/`shared_demux`, `:651`) and `SEM_UNDO`
process-counting (`:764`) are **owned per process**, so two threads in one
process bypass them entirely (thread B's `F_SETLKW` succeeds even while thread A
"holds" it). Real lock-protocol validation requires separate processes — see the
next section.

---

## Cross-process test — the only thing that validates real sharing (after Phase 2)

**Why it can't be a thread test.** As noted above, the driver's entire
synchronization mechanism is process-scoped by design. Threads share the
process's single attachment and (post-Phase-2) are serialized by the in-process
mutex, so a thread test can never exercise `fcntl` contention, the semaphore
process-count, `SEM_UNDO` cleanup on process death, or the master/slave
global-table creation race. Those only happen across real processes.

**Approach — spawn a helper binary (robust, CI-friendly).** Prefer a small
test-only binary driven via `std::process::Command` over `libc::fork()`: the
Rust test harness is multi-threaded, and `fork()` in a multi-threaded process
only permits async-signal-safe calls before `exec` — fragile. An `exec`'d helper
gets a clean address space and genuinely separate IPC state.

- [ ] Add a feature-gated helper bin, e.g. `src/bin/shmem_helper/main.rs`, gated
      `#[cfg(all(feature = "shared_mem", not(windows)))]`, exposing one op per
      invocation with a byte pattern that's easy to assert:
      - `shmem_helper create <hN> <len> <seed>` — create segment, fill `len`
        bytes with a deterministic pattern derived from `seed`, close (leaving it
        PERSIST so the parent can read it).
      - `shmem_helper read <hN> <len> <seed>` — open, read `len` bytes, compare to
        the same pattern; exit 0 on match, non-zero on mismatch/error.
      - `shmem_helper delete <hN>` — `smem_remove` for teardown.
      It must honor `SHMEM_LIB_KEYBASE`/`SHMEM_LIB_MAXSEG` from the environment
      (the driver already reads them at init, `:373`) so the parent can put the
      child in a private namespace.
- [ ] Add `tests/test_drvrsmem_xproc.rs` (gated the same way). Each test:
      1. picks a **unique `SHMEM_LIB_KEYBASE`** (e.g. base + a per-test offset)
         and passes it to every child via `Command::env`, isolating this test
         from parallel tests, other CI jobs, and stray segments on the runner;
      2. spawns `create` in one child, waits for exit 0;
      3. spawns `read` in a *second* child, asserts exit 0 — proving the bytes
         written by process A are visible to process B through the segment;
      4. optionally: hold the segment locked in one child (a `hold <hN> <ms>`
         op) and assert a second child's `SHARED_NOWAIT` open returns
         `SHARED_AGAIN`, to actually exercise cross-process fcntl contention;
      5. **teardown that always runs** (even on assertion failure — use a drop
         guard or `std::panic::catch_unwind` + rethrow): `delete <hN>`, and as a
         backstop remove any segment under this test's keybase.
- [ ] Locate the helper binary at test time via the `CARGO_BIN_EXE_shmem_helper`
      env var that Cargo sets for integration tests, so you don't hard-code paths.

**Cleanup discipline (must-have).** SysV segments persist past process exit until
`IPC_RMID`; a leaked PERSIST segment collides on key re-use and slowly exhausts
the kernel `shmmni` limit on a persistent runner. Every cross-process test must
delete its segments in an always-run teardown, and CI should have an
`if: always()` backstop (`ipcs -m` for visibility; targeted removal by this
run's keybase). The unique-keybase-per-test rule makes both the isolation and
the cleanup targetable.

**Ordering.** Write the helper bin any time, but only enable the cross-process
*assertions* after Phases 1–2 land — before that, the concurrent global-table
access is itself UB and the results are meaningless. Add this section's tests to
the Phase 0 CI job once they're green.

**Verify.** `create` in one process, `read` in another, bytes match; the
optional lock-contention test observes `SHARED_AGAIN`; no segments leak
(`ipcs -m` clean after the run).

---

## Phase 3 — Remove the `&[c_char]` → `&mut [c_char]` cast (Tier 1)

**Problem.** `smem_remove` casts an immutable slice to mutable at
`src/drvrsmem.rs:1812-1814` (self-flagged "SAFETY: Absolutely none") only to
pass it to `smem_open`, which never writes it — it just `sscanf_d`s `h%d`.

**Change.**
- [ ] Change `smem_open`'s first parameter from `&mut [c_char]` to `&[c_char]`
      (`:1673`). Confirm the body only reads (`sscanf_d` reads; the rest uses the
      parsed `h`). Update `driverhandle` out-param stays `&mut c_int`.
- [ ] Delete the `slice::from_raw_parts_mut` cast in `smem_remove`; pass
      `filename` directly.
- [ ] Fix call sites: `smem_open` is referenced in the driver registration
      (`src/cfileio.rs:4760`) via `Some(smem_open)` — the registration takes a
      function pointer; confirm the expected signature there. If the driver
      table's `open` slot types the buffer as `&mut`, either keep the registered
      shim `&mut` and have it call an inner `&[..]` parser, or update the table
      type if it's local to this driver. Check `fits_register_driver`'s expected
      `open` fn type before changing the public signature.
- [ ] Update `shared_getaddr` (`:1576`) which calls `smem_open(&mut segname, ...)`
      — it can pass `&segname`.

**Verify.** Build + tests green. `test_delete_shmem` and `test_smem_remove`
specifically exercise this path.

---

## Phase 4 — Contain panics across the FFI boundary (Tier 2)

**Problem.** Public fns are `extern "C"`; a panic unwinding out of them is UB.
The body uses `.try_into().unwrap()`, `SHARED_LT[idx].p.unwrap()`, and panicking
slice indexing throughout. Most `p.unwrap()` sites are guarded by a prior
`shared_check_locked_index`, but the Tier-1 races can violate those invariants.

**Change.**
- [ ] Triage each `try_into().unwrap()` by what can make it fail:
      - **Untrusted / racy / environmental input** (values derived from a caller
        argument, a cross-process atomic read, or `getenv`) → return the
        function's `SHARED_*` error code; never `unwrap`.
      - **Genuinely impossible given an invariant** (e.g. converting
        `SHARED_MAXSEG` — already validated `> 0` and bounded at init — to
        `usize`) → use `.expect("SHARED_MAXSEG is validated positive at init")`.
        This documents the invariant; if it ever fires it's a real logic bug, and
        C would be equally broken at that point.
      Don't blanket-`unwrap`, and don't reflexively convert every one to an error
      return either — pick per the two cases above.
- [ ] Convert `SHARED_LT[idx].p.unwrap()` on I/O paths (`smem_read` `:1907`,
      `smem_write` `:1963`, `shared_getaddr` `:1580`) to a checked pattern that
      returns `SHARED_BADARG`/`SHARED_IPCERR` / null on `None`.
- [ ] Wrap each public `extern "C"` entry point body in `catch_unwind`, mapping a
      caught panic to the function's error sentinel (`SHARED_IPCERR` for
      `c_int`-returning, `ptr::null()` for `SHARED_P`-returning). Note
      `catch_unwind` requires the closure's captures be `UnwindSafe`; raw
      pointers and the mutex guard may need `AssertUnwindSafe` — acceptable here
      because on panic we abandon the operation and return an error.
- [ ] `shared_cleanup` (`:220`, registered via `atexit` at `:534`) is the most
      important one to wrap — it runs at process teardown, touches all global
      state, and does `print!`; a broken-pipe panic there is UB. Wrap its body in
      `catch_unwind` and swallow.

**Verify.** Build + tests green. Add a test that triggers an error path (e.g.
open a non-existent handle, read beyond EOF — `test_read_beyond_eof`,
`test_smem_read_beyond_eof` already exist) and confirm it returns a code rather
than aborting.

---

## Phase 5 — Ported-from-C logic bugs with a safety dimension (Tier 3)

### 5a. Dead `nbytes < 0` checks + real overflow risk
- [ ] `smem_read` (`:1889`) and `smem_write` (`:1939`) check `nbytes < 0` on a
      `usize` — always false, dead. Remove.
- [ ] The live risk is `SHARED_LT[idx].seekpos + nbytes as c_long` (`:1899`,
      `:1943`) overflowing `c_long` and defeating the EOF bound before the raw
      `memcpy` (`:1905`, `:1962`). Use `checked_add` (as `i64`/`c_long`) and
      return `SHARED_BADARG` on overflow **before** the bound comparison.

### 5b. `SHARED_INVALID as *mut BLKHEAD` shmat-failure sentinel
- [ ] `shmat` returns `(void*)-1` on failure; compared at `:933`, `:1100`,
      `:1227` against `SHARED_INVALID as *mut BLKHEAD`. Correctness relies on
      `-1i32 as *mut T` sign-extending to `0xFFFF…FFFF` rather than
      zero-extending. Replace with an explicit sentinel:
      `const SHM_FAILED: *mut c_void = usize::MAX as *mut c_void;` (or
      `!0usize as *mut _`) and compare against that after casting. Do the same
      for the `SHARED_GTAB`/`BLKHEAD` variants. This removes reliance on int→ptr
      cast semantics.

### 5c. `shared_getaddr` dangling pointer + leaked lock
- [ ] `shared_getaddr` (`:1580`) returns an interior pointer with `smem_close(i)`
      commented out (`:1582`): the segment stays locked/attached (lock leak — the
      `smem` binary loops it over 16 ids and never frees) and the returned
      pointer dangles if another process reallocs (realloc swaps the SysV segment
      and detaches, `:1249`).
- [ ] This one is a **behavior/contract** question, not a mechanical fix: cfitsio
      leaves it this way deliberately (the caller is expected to hold the lock).
      Do **not** silently change semantics. Instead: add a doc-comment stating the
      lifetime/lock contract explicitly, and note the leak in `src/bin/smem`.
      Decide with the maintainer whether to change the return contract before
      touching it. Leave as-is + documented unless told otherwise.

**Verify.** Build + tests green; overflow/EOF tests still pass.

---

## Phase 6 — Optional cleanups (only if time allows, lowest priority)

- [ ] Replace hand-rolled `memcpy` + pointer arithmetic (`:1905`, `:1962`,
      `:1242`) with `slice::from_raw_parts(_mut)` + `copy_from_slice` once bounds
      are proven — gives length-checked copies in debug builds. Behavior must stay
      identical; verify byte output unchanged.
- [ ] Pin the `union semun` / variadic `semctl` ABI assumption (`:128`, `:746`)
      with a comment referencing the target platforms.
- [ ] Once Phases 1–2 land, try removing `#![allow(static_mut_refs)]` (`:1`) and
      fix or narrowly re-scope whatever remains.

---

## Explicit non-goals

- Not migrating SysV IPC → POSIX shm/`mmap`. That changes the cross-process wire
  format and breaks compatibility with the C library and existing segments.
- Not changing the `shmem://hN` URL grammar, driver registration order, error
  codes, or `smem` CLI behavior.
- Not enabling the feature by default or claiming it "fully works" — per
  CLAUDE.md it is a known-incomplete optional feature.

## Suggested commit sequence

1. Phase 1 — raw-pointer + volatile global table.
2. Phase 2 — single mutex over local state (+ `*_locked` inner fns).
3. Phase 3 — `smem_open` takes `&[c_char]`, drop the unsound cast.
4. Phase 4 — `catch_unwind` wrappers + checked error returns.
5. Phase 5a/5b — overflow guard + explicit shmat sentinel.
6. Phase 5c — documentation-only (or deferred pending maintainer decision).
7. Phase 6 — optional, separate commits.

Each commit: `cargo build --features shared_mem`,
`cargo test --features shared_mem -- --test-threads=1`, `cargo fmt`,
`cargo clippy --features shared_mem`.
