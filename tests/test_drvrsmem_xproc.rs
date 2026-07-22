//! Cross-process tests for the shared-memory driver.
//!
//! SKETCH (see HARDENING_PLAN_drvrsmem.md). These validate the ONE thing a
//! thread test cannot: that data written by one process is visible to another
//! through a SysV segment, and that the fcntl lock protocol actually contends
//! across processes. The driver's fcntl record locks and SEM_UNDO counting are
//! process-scoped, so only separate processes exercise them.
//!
//! Strategy: spawn the `shmem_helper` binary (exec'd, clean address space — no
//! fork-in-a-threaded-harness hazard) via Command. Each test uses a UNIQUE
//! `SHMEM_LIB_KEYBASE` so it can't collide with parallel tests, other CI jobs,
//! or stray segments on the runner. Teardown always runs.
//!
//! Enable only after Phases 1–2 land — before that the driver's global-table
//! access is itself UB and results are meaningless. Run with:
//!   cargo test --features shared_mem --test test_drvrsmem_xproc

#![cfg(all(feature = "shared_mem", not(windows)))]

use std::process::Command;

/// Path to the helper binary, injected by Cargo for integration tests.
const HELPER: &str = env!("CARGO_BIN_EXE_shmem_helper");

/// Run one helper op in its own process, in the given key namespace.
fn run(keybase: i32, args: &[&str]) -> i32 {
    let status = Command::new(HELPER)
        .args(args)
        .env("SHMEM_LIB_KEYBASE", keybase.to_string())
        .env("SHMEM_LIB_MAXSEG", "16")
        .status()
        .expect("failed to spawn shmem_helper");
    status.code().unwrap_or(-1)
}

/// Always-run teardown: delete the segment even if an assertion panicked.
///
/// Note: deleting the FITS *data* segment (via the helper) does NOT remove the
/// per-keybase *global table* segment (key == keybase), which the driver only
/// unlinks under a last-attacher race it usually loses across short-lived helper
/// processes. So we also `ipcrm -M <keybase>` as a backstop — otherwise each
/// test leaks one 448-byte SysV segment and eventually hits the `shmmni` limit.
struct Teardown {
    keybase: i32,
    handle: &'static str,
}
impl Drop for Teardown {
    fn drop(&mut self) {
        // Best-effort; ignore results (may already be gone).
        let _ = run(self.keybase, &["delete", self.handle]);
        let _ = Command::new("ipcrm")
            .args(["-M", &self.keybase.to_string()])
            .status();
    }
}

#[test]
fn write_in_one_process_read_in_another() {
    // Unique namespace for this test (base + offset). Keep offsets distinct
    // across tests in this file so parallel test threads don't share a keybase.
    let keybase = 24_010_001;
    let handle = "h0";
    let seed = "1234";
    let _teardown = Teardown { keybase, handle };

    // Process A creates + fills the segment, then exits (segment is PERSIST).
    assert_eq!(
        run(keybase, &["create", handle, seed]),
        0,
        "writer process failed"
    );

    // Process B — a genuinely separate process — reads and verifies the bytes.
    assert_eq!(
        run(keybase, &["read", handle, seed]),
        0,
        "reader process did not see the writer's data across the segment"
    );
}

#[test]
fn wrong_seed_is_detected() {
    // Sanity check the oracle: reading with a different seed must fail (exit 3),
    // proving the read genuinely compares cross-process bytes.
    let keybase = 24_020_001;
    let handle = "h0";
    let _teardown = Teardown { keybase, handle };

    assert_eq!(run(keybase, &["create", handle, "1"]), 0);
    assert_eq!(
        run(keybase, &["read", handle, "999"]),
        3,
        "reader should have reported a mismatch"
    );
}

// TODO (after Phase 2, when a `hold <hN> <ms>` op is added to the helper):
// cross-process lock contention — process A holds the segment RW-locked while
// process B attempts a SHARED_NOWAIT open and must observe SHARED_AGAIN. This
// is the test that actually exercises the fcntl protocol between processes.
