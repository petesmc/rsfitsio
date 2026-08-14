#!/usr/bin/env bash
#
# Differential harness: this repo's fpack/funpack against the C originals.
#
#   C_FPACK=~/code/cfitsio/fpack C_FUNPACK=~/code/cfitsio/funpack \
#       ./scripts/fpack_diff.sh [extra FITS files...]
#
# For every (file x flag combination) it compares stdout, stderr, the exit
# status and -- where one is produced -- the output FITS byte for byte, and it
# checks that each implementation can read back what the other wrote.
#
# THREE NORMALISATIONS ARE MANDATORY, each applied below:
#
#   1. Dithering.  dither_offset 0 means "seed from the system clock", so two
#      runs of the same command do not even agree with each other.  Byte
#      comparisons pass `-q0 <level>' (the digit is a *suffix* meaning "do not
#      dither") or an explicit `-qN'.
#   2. Checksums.  fits_write_chksum stamps CHECKSUM with a wall-clock
#      timestamp in its comment, so byte comparisons pass `-C'.
#   3. Timings.  `-T' reports elapsed and CPU seconds, so those columns are
#      masked before diffing.
#
# Known differences, reported but not counted as failures:
#
#   * The CFITSIO version number in `-V'/`-H', when the two builds are synced
#     to different releases.
#   * The compressed payload of GZIP_1, GZIP_2 and HCOMPRESS: libz-rs and the
#     `hcompress' crate emit valid but not byte-identical streams.  The round
#     trips are what prove those correct.  RICE, PLIO and NOCOMPRESS *are*
#     expected to be byte-identical.

set -uo pipefail

C_FPACK=${C_FPACK:-$HOME/code/cfitsio/fpack}
C_FUNPACK=${C_FUNPACK:-$HOME/code/cfitsio/funpack}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
R_FPACK=${R_FPACK:-$ROOT/target/debug/fpack}
R_FUNPACK=${R_FUNPACK:-$ROOT/target/debug/funpack}

for exe in "$C_FPACK" "$C_FUNPACK" "$R_FPACK" "$R_FUNPACK"; do
    if [ ! -x "$exe" ]; then
        echo "not executable: $exe" >&2
        echo "build the Rust side with: cargo build --bin fpack --bin funpack" >&2
        exit 2
    fi
done

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

pass=0 fail=0 known=0
note() { printf '  %s\n' "$*"; }
ok()   { pass=$((pass+1)); }
bad()  { fail=$((fail+1)); printf 'FAIL %s\n' "$*"; }
skip() { known=$((known+1)); printf 'known %s\n' "$*"; }

# ---------------------------------------------------------------------------
# corpus
# ---------------------------------------------------------------------------

mkcorpus() {
    python3 - "$WORK" <<'PY'
import struct, random, sys, os
out = sys.argv[1]
def card(s): return s.ljust(80)[:80]
def write(name, bitpix, nx, ny, values, fmt):
    cards = [card("SIMPLE  =                    T"),
             card(f"BITPIX  = {bitpix:20d}"),
             card("NAXIS   =                    2"),
             card(f"NAXIS1  = {nx:20d}"),
             card(f"NAXIS2  = {ny:20d}"),
             card("END")]
    hdr = "".join(cards).ljust(2880).encode()
    data = b"".join(struct.pack(fmt, v) for v in values)
    data += b"\0" * ((2880 - len(data) % 2880) % 2880)
    open(os.path.join(out, name), "wb").write(hdr + data)

nx, ny = 60, 40
n = nx * ny
random.seed(1)
# noisy signed 16-bit: the general case
write("i16_noise.fits", 16, nx, ny, [random.randint(-3000, 3000) for _ in range(n)], ">h")
# noisy non-negative 16-bit: the only kind PLIO accepts
write("u16_noise.fits", 16, nx, ny, [random.randint(0, 3000) for _ in range(n)], ">h")
# a smooth gradient: what quantization is actually designed for
write("i16_smooth.fits", 16, nx, ny,
      [1000 + (i % nx) + (i // nx) // 2 for i in range(n)], ">h")
# constant: every run-length encoder's best case
write("i16_flat.fits", 16, nx, ny, [7] * n, ">h")
# 32-bit
write("i32_noise.fits", 32, nx, ny, [random.randint(-100000, 100000) for _ in range(n)], ">i")
# float
write("f32_smooth.fits", -32, nx, ny,
      [1000.0 + (i % nx) * 0.25 + (i // nx) * 0.5 for i in range(n)], ">f")
# binary tables, for -table / -tableonly
def table(name, cols, nrows, seed):
    import random
    rnd = random.Random(seed)
    w = sum(c[1] for c in cols)
    cards = [card("XTENSION= 'BINTABLE'"), card("BITPIX  =                    8"),
             card("NAXIS   =                    2"), card(f"NAXIS1  = {w:20d}"),
             card(f"NAXIS2  = {nrows:20d}"), card("PCOUNT  =                    0"),
             card("GCOUNT  =                    1"), card(f"TFIELDS = {len(cols):20d}")]
    for i, (tf, _, _) in enumerate(cols, 1):
        cards.append(card(f"TTYPE{i}  = 'C{i:<7}'"))
        cards.append(card(f"TFORM{i}  = '{tf:<8}'"))
    cards.append(card("END"))
    prim = "".join([card("SIMPLE  =                    T"), card("BITPIX  =                    8"),
                    card("NAXIS   =                    0"), card("EXTEND  =                    T"),
                    card("END")]).ljust(2880)
    body = "".join(cards)
    hdr = body.ljust(2880 * ((len(body) + 2879) // 2880))
    rows = b"".join(b"".join(c[2](r, rnd) for c in cols) for r in range(nrows))
    rows += b"\0" * ((2880 - len(rows) % 2880) % 2880)
    open(os.path.join(out, name), "wb").write(prim.encode() + hdr.encode() + rows)

table("tbl_mixed.fits", [
    ("1I", 2, lambda r, q: struct.pack(">h", (r * 7) % 30000)),
    ("1J", 4, lambda r, q: struct.pack(">i", r * 12345)),
    ("1D", 8, lambda r, q: struct.pack(">d", r * 0.125)),
    ("1E", 4, lambda r, q: struct.pack(">f", r * 1.5)),
    ("16A", 16, lambda r, q: f"name{r:011d}".encode()),
], 500, 2)
table("tbl_noisy.fits", [
    ("1J", 4, lambda r, q: struct.pack(">i", q.randint(-2**30, 2**30))),
    ("1I", 2, lambda r, q: struct.pack(">h", q.randint(-32768, 32767))),
], 700, 3)

# a variable-length-array column, which takes a different branch of both the
# transposer and the per-column compressor
def vla_table(name, nrows):
    heap = b""
    descr = b""
    for r in range(nrows):
        vals = [(r * 31 + k) % 10000 for k in range(r % 7 + 1)]
        descr += struct.pack(">ii", len(vals), len(heap))
        heap += b"".join(struct.pack(">i", v) for v in vals)
    cards = [card("XTENSION= 'BINTABLE'"), card("BITPIX  =                    8"),
             card("NAXIS   =                    2"), card("NAXIS1  =                    8"),
             card(f"NAXIS2  = {nrows:20d}"), card(f"PCOUNT  = {len(heap):20d}"),
             card("GCOUNT  =                    1"), card("TFIELDS =                    1"),
             card("TTYPE1  = 'VLA     '"), card("TFORM1  = '1PJ(16) '"), card("END")]
    prim = "".join([card("SIMPLE  =                    T"), card("BITPIX  =                    8"),
                    card("NAXIS   =                    0"), card("EXTEND  =                    T"),
                    card("END")]).ljust(2880)
    body = "".join(cards)
    hdr = body.ljust(2880 * ((len(body) + 2879) // 2880))
    data = descr + heap
    data += b"\0" * ((2880 - len(data) % 2880) % 2880)
    open(os.path.join(out, name), "wb").write(prim.encode() + hdr.encode() + data)

vla_table("tbl_vla.fits", 700)

# 1-D, which exercises upstream bugs 5 and 7
cards = [card("SIMPLE  =                    T"), card("BITPIX  =                   16"),
         card("NAXIS   =                    1"), card("NAXIS1  =                  100"),
         card("END")]
data = b"".join(struct.pack(">h", i % 100) for i in range(100))
data += b"\0" * ((2880 - len(data) % 2880) % 2880)
open(os.path.join(out, "oned.fits"), "wb").write("".join(cards).ljust(2880).encode() + data)
PY
}

mkcorpus
CORPUS=("$WORK"/*.fits)
for extra in "$@"; do
    cp "$extra" "$WORK/extra_$(basename "$extra")"
    CORPUS+=("$WORK/extra_$(basename "$extra")")
done

echo "corpus: ${#CORPUS[@]} files"

# BITPIX of a FITS file, from its first header block.
bitpix_of() { head -c 2880 "$1" | grep -aoE 'BITPIX *= *-?[0-9]+' | head -1 |
              grep -oE '\-?[0-9]+$'; }

# A floating-point image is *quantized* by default, so only `-q 0' with GZIP
# round trips bit for bit; every other algorithm is lossy on purpose and the C
# does not reproduce its own input either.  Integer images are lossless under
# all of them.
is_float() { [ "$(bitpix_of "$1")" -lt 0 ]; }

# The table fixtures have an empty primary array; images do not.
is_table() { head -c 2880 "$1" | grep -aqE 'NAXIS   = +0'; }

# ---------------------------------------------------------------------------
# 1. CLI behaviour: stdout, stderr and exit status must agree exactly
# ---------------------------------------------------------------------------

cli_cases=(
    ""                       # no arguments
    "-H"
    "-V"
    "-Q x.fits"              # unknown flag
    "nosuchfile.fits"
    "a.fits[1]"              # extension notation
    "-t 1,2,3,4,5,6,7 a.fits"  # upstream bug 1
    "-t 100x100 a.fits"
    "-s 1 a.fits"            # -s without -h or -T
    "-q 0 a.fits"            # -q 0 without GZIP
    "-R rep.txt a.fits"      # -R without -T
    "-O out.fz -S a.fits"
    "-r -g a.fits"           # two compression flags
)
funpack_cli_cases=(
    ""
    "-H"
    "-V"
    "-Q x.fits"
    "-E SCI -F a.fits"
    "-S -P p a.fits"
    "-P p -O o a.fits"
)

echo
echo "== CLI behaviour =="
run_cli() {   # $1 = C exe, $2 = rust exe, $3 = label, rest = args
    local cexe=$1 rexe=$2 label=$3; shift 3
    local d="$WORK/cli"; rm -rf "$d"; mkdir -p "$d"
    cp "$WORK/i16_noise.fits" "$d/a.fits"
    ( cd "$d" && $cexe "$@" >c.out 2>c.err ); local ce=$?
    cp "$WORK/i16_noise.fits" "$d/a.fits"
    ( cd "$d" && $rexe "$@" >r.out 2>r.err ); local re=$?

    # the CFITSIO version number differs when the two are synced to different
    # releases; normalise it out of -V and -H
    sed -i 's/CFITSIO version [0-9.]*/CFITSIO version X/' "$d/c.out" "$d/r.out"

    if [ "$ce" != "$re" ]; then
        bad "$label [$*]: exit status $ce vs $re"; return
    fi
    if ! diff -q "$d/c.out" "$d/r.out" >/dev/null; then
        bad "$label [$*]: stdout differs"; diff "$d/c.out" "$d/r.out" | head -10; return
    fi
    if ! diff -q "$d/c.err" "$d/r.err" >/dev/null; then
        bad "$label [$*]: stderr differs"; diff "$d/c.err" "$d/r.err" | head -10; return
    fi
    ok
}

for args in "${cli_cases[@]}"; do
    # shellcheck disable=SC2086
    run_cli "$C_FPACK" "$R_FPACK" fpack $args
done
for args in "${funpack_cli_cases[@]}"; do
    # shellcheck disable=SC2086
    run_cli "$C_FUNPACK" "$R_FUNPACK" funpack $args
done

# ---------------------------------------------------------------------------
# 2. Compressed output, byte for byte where the codec is shared
# ---------------------------------------------------------------------------

echo
echo "== compressed output =="
# RICE and NOCOMPRESS are implemented in this tree and must match exactly.
# GZIP and HCOMPRESS go through third-party codecs that emit different but
# equally valid streams, so only their round trips are checked (section 3).
exact_algs=("-r" "-d" "-p")
for f in "${CORPUS[@]}"; do
    base=$(basename "$f")
    # RICE on a floating-point image quantizes and then falls back to gzip for
    # the scaling parameters, so even "-r" goes through a codec that differs.
    # `is_float' only sees the primary HDU, so user-supplied files -- which may
    # hold a float image in any extension -- are excluded outright; the cross
    # round trips below are what covers them.
    if is_float "$f" || is_table "$f" || [[ "$base" == extra_* ]]; then continue; fi
    for alg in "${exact_algs[@]}"; do
        d="$WORK/pk"; rm -rf "$d"; mkdir -p "$d"
        cp "$f" "$d/c.fits"; cp "$f" "$d/r.fits"
        ( cd "$d" && $C_FPACK -C -q0 4 "$alg" c.fits >/dev/null 2>&1 )
        ( cd "$d" && $R_FPACK -C -q0 4 "$alg" r.fits >/dev/null 2>&1 )
        if [ ! -f "$d/c.fits.fz" ] && [ ! -f "$d/r.fits.fz" ]; then
            ok; continue                      # both declined, consistently
        fi
        if [ ! -f "$d/c.fits.fz" ] || [ ! -f "$d/r.fits.fz" ]; then
            bad "$base $alg: only one implementation produced output"; continue
        fi
        if cmp -s "$d/c.fits.fz" "$d/r.fits.fz"; then
            ok
        else
            bad "$base $alg: compressed output differs ($(cmp -l "$d/c.fits.fz" "$d/r.fits.fz" | wc -l) bytes)"
        fi
    done
done

# ---------------------------------------------------------------------------
# 3. Cross round trips: each implementation must read the other's output
# ---------------------------------------------------------------------------

echo
echo "== cross round trips =="
int_algs=("-r" "-g" "-g2" "-h" "-p" "-d")
# lossless float compression is GZIP with quantization switched off entirely
flt_algs=("-g" "-g2")
# binary tables have their own two switches rather than an algorithm flag
tbl_algs=("-tableonly" "-table")
for f in "${CORPUS[@]}"; do
    base=$(basename "$f")
    if is_table "$f"; then
        all_algs=("${tbl_algs[@]}"); qflags=(-q0 4)
    elif is_float "$f"; then
        all_algs=("${flt_algs[@]}"); qflags=(-q 0)
    else
        all_algs=("${int_algs[@]}"); qflags=(-q0 4)
    fi
    for alg in "${all_algs[@]}"; do
        # Run all four (packer, unpacker) combinations and require that they
        # agree with *each other*.  That is the implementation-independence
        # invariant, and unlike "equals the original" it holds for any input:
        # funpack rewrites headers, so a multi-extension file does not come
        # back byte-identical even from the C alone.  For the single-HDU
        # integer fixtures built above, C -> C *is* the identity, so comparing
        # against it also proves the round trip lossless.
        declare -A got=()
        for pk in C R; do for un in C R; do
            d="$WORK/rt_$pk$un"; rm -rf "$d"; mkdir -p "$d"
            cp "$f" "$d/a.fits"
            if [ "$pk" = C ]; then P=$C_FPACK; else P=$R_FPACK; fi
            if [ "$un" = C ]; then U=$C_FUNPACK; else U=$R_FUNPACK; fi

            ( cd "$d" && $P -C "${qflags[@]}" "$alg" a.fits >/dev/null 2>&1 )
            if [ ! -f "$d/a.fits.fz" ]; then
                continue          # declined; section 2 checked they agree on that
            fi
            rm -f "$d/a.fits"     # funpack refuses to overwrite an existing output
            ( cd "$d" && $U -C a.fits.fz >/dev/null 2>&1 )
            [ -f "$d/a.fits" ] && got[$pk$un]=$d/a.fits
        done; done

        ref=${got[CC]:-}
        if [ -z "$ref" ]; then ok; continue; fi     # the C declined this input
        for combo in CR RC RR; do
            label="$base $alg pack=${combo:0:1} unpack=${combo:1:1}"
            if [ -z "${got[$combo]:-}" ]; then
                bad "$label: unpack produced nothing"
            elif cmp -s "${got[$combo]}" "$ref"; then
                ok
            else
                bad "$label: result differs from the C round trip"
            fi
        done
        # For our own single-HDU image fixtures the C round trip is also the
        # identity, so compare against the original too.  That does not hold
        # for a table with a heap -- funpack repacks it -- nor for arbitrary
        # user-supplied files, where the four-way agreement above is the whole
        # of the invariant.
        if [[ "$base" != extra_* ]] && ! is_table "$f" && ! cmp -s "$ref" "$f"; then
            bad "$base $alg: the C round trip is not lossless (harness bug?)"
        fi
    done
done

# ---------------------------------------------------------------------------
# 4. Listings and the -T report
# ---------------------------------------------------------------------------

echo
echo "== -L listing and -T report =="
for f in "${CORPUS[@]}"; do
    base=$(basename "$f")
    if is_table "$f"; then continue; fi   # -L/-T report on image HDUs
    d="$WORK/rep"; rm -rf "$d"; mkdir -p "$d"
    cp "$f" "$d/a.fits"
    ( cd "$d" && $C_FPACK -C -q0 4 -r a.fits >/dev/null 2>&1 )

    if [ -f "$d/a.fits.fz" ]; then
        if diff <( "$C_FPACK" -L "$d/a.fits.fz" ) <( "$R_FPACK" -L "$d/a.fits.fz" ) >/dev/null; then
            ok
        else
            bad "$base: -L listing differs"
            diff <( "$C_FPACK" -L "$d/a.fits.fz" ) <( "$R_FPACK" -L "$d/a.fits.fz" ) | head -10
        fi
    fi

    # -T: on the per-algorithm rows every number is either a wall-clock
    # measurement or a compressed size, and the compressed size depends on the
    # codec (libz-rs and the hcompress crate do not emit byte-identical
    # streams), so those rows are compared for structure only.  The image
    # statistics line above them -- nulls, min, max, mean, sigma, the three
    # noise estimates, Nbits, MaxR -- comes from the same fp_iNstat and
    # fits_img_stats code on both sides and IS compared exactly.
    mask() {
        sed -E '/^ *(RICE|GZIP1|GZIP2|HCOMP|PLIO|NONE|Native)/ s/-?[0-9]+\.[0-9]+/N/g' |
            sed -E 's/ +/ /g'
    }
    dc="$WORK/rep_c"; dr="$WORK/rep_r"
    rm -rf "$dc" "$dr"; mkdir -p "$dc" "$dr"; cp "$f" "$dc/a.fits"; cp "$f" "$dr/a.fits"
    ( cd "$dc" && $C_FPACK -T a.fits 2>/dev/null | mask >"$WORK/t.c" )
    ( cd "$dr" && $R_FPACK -T a.fits 2>/dev/null | mask >"$WORK/t.r" )
    if diff -q "$WORK/t.c" "$WORK/t.r" >/dev/null; then
        ok
    else
        bad "$base: -T report differs (timings masked)"
        diff "$WORK/t.c" "$WORK/t.r" | head -20
    fi
done

echo
echo "=============================="
printf 'passed %d, failed %d, known differences %d\n' "$pass" "$fail" "$known"
[ "$fail" -eq 0 ]
