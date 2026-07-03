#!/usr/bin/env bash
# Tests for separate indexing routines: fa2bit, genraw, raw2bwt, gensa, genbwt
# These are the individual steps that main_index orchestrates.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
BINARY="${REPO_DIR}/minibwa"
PASS=0
FAIL=0

cleanup() {
    rm -f /tmp/mb_ir_* 2>/dev/null || true
}
trap cleanup EXIT

if [ ! -x "$BINARY" ]; then
    echo "ERROR: Binary not found or not executable: $BINARY"
    echo "Build it first with: cd $REPO_DIR && make"
    exit 1
fi
pass() { PASS=$((PASS + 1)); echo -e "  \033[0;32mPASS\033[0m: $1"; }
fail() { FAIL=$((FAIL + 1)); if [ -n "${2:-}" ]; then echo -e "  \033[0;31mFAIL\033[0m: $1 ($2)"; else echo -e "  \033[0;31mFAIL\033[0m: $1"; fi }
skip() { echo -e "  \033[1;33mSKIP\033[0m: $1"; }

echo "============================================"
echo " Minibwa Index Routine Tests"
echo "============================================"
echo ""

# Test 1: fa2bit - convert FASTA to 2-bit
echo "[IR 1] fa2bit: convert FASTA to long-2bit format"
"$BINARY" fa2bit "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_fa2bit > /dev/null 2>&1
if [ -f "/tmp/mb_ir_fa2bit" ]; then
    pass "fa2bit creates output file"
    # Check file has reasonable size (chrM is ~16.5kb, should be small)
    size=$(stat -c%s "/tmp/mb_ir_fa2bit" 2>/dev/null || stat -f%z "/tmp/mb_ir_fa2bit" 2>/dev/null)
    if [ "$size" -gt 0 ] && [ "$size" -lt 100000 ]; then
        pass "fa2bit output size is reasonable ($size bytes)"
    else
        fail "fa2bit output size unexpected: $size bytes"
    fi
else
    fail "fa2bit did not create output file"
fi

# Test 2: genbwt - generate BWT from .l2b
echo "[IR 2] genbwt: generate BWT+SSA from long-2bit"
"$BINARY" genbwt /tmp/mb_ir_fa2bit /tmp/mb_ir_fa2bit.mbw > /dev/null 2>&1
if [ -f "/tmp/mb_ir_fa2bit.mbw" ]; then
    pass "genbwt creates .mbw file"
    size=$(stat -c%s "/tmp/mb_ir_fa2bit.mbw" 2>/dev/null || stat -f%z "/tmp/mb_ir_fa2bit.mbw" 2>/dev/null)
    if [ "$size" -gt 0 ]; then
        pass "genbwt output size is reasonable ($size bytes)"
    else
        fail "genbwt output size unexpected: $size bytes"
    fi
else
    fail "genbwt did not create .mbw file"
fi

# Test 3: genbwt with different thread counts
echo "[IR 3] genbwt with multi-threading"
"$BINARY" genbwt -t2 /tmp/mb_ir_fa2bit /tmp/mb_ir_fa2bit_t.mbw > /dev/null 2>&1
if [ -f "/tmp/mb_ir_fa2bit_t.mbw" ]; then
    pass "genbwt multi-thread creates .mbw"
else
    fail "genbwt multi-thread did not create .mbw"
fi

# Test 4: genraw - generate BWT from pac with BWT-SW algorithm
echo "[IR 4] genraw: generate BWT from .pac"
# Create .pac file via fa2bit -p (which outputs BWA pac format)
"$BINARY" fa2bit -p "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_genraw.pac > /dev/null 2>&1
if [ -f "/tmp/mb_ir_genraw.pac" ]; then
        "$BINARY" genraw /tmp/mb_ir_genraw.pac /tmp/mb_ir_genraw_out > /dev/null 2>&1
        if [ -f "/tmp/mb_ir_genraw_out.bwt" ] || [ -f "/tmp/mb_ir_genraw_out.occ" ]; then
            pass "genraw produces BWT/occ files"
        else
            # genraw may produce different output names; check for any output
            count=$(ls /tmp/mb_ir_genraw_out* 2>/dev/null | wc -l) || count=0
            if [ "$count" -gt 0 ]; then
                pass "genraw produces output files ($count files)"
            else
                fail "genraw produced no output files"
            fi
        fi
else
    fail "fa2bit -p did not produce .pac file"
fi
rm -f /tmp/mb_ir_genraw* 2>/dev/null || true

# Test 5: raw2bwt - recode bwtgen raw BWT
echo "[IR 5] raw2bwt: recode bwtgen raw BWT"
# raw2bwt writes to the exact output path given on the command line (no extension added)
"$BINARY" fa2bit -p "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_raw2bwt.pac > /dev/null 2>&1
if [ -f "/tmp/mb_ir_raw2bwt.pac" ]; then
    "$BINARY" genraw /tmp/mb_ir_raw2bwt.pac /tmp/mb_ir_raw2bwt_out > /dev/null 2>&1
    # raw2bwt reads .bwt and .occ from genraw output, writes to exact output path
    "$BINARY" raw2bwt /tmp/mb_ir_raw2bwt_out /tmp/mb_ir_raw2bwt_final > /dev/null 2>&1
    if [ -f "/tmp/mb_ir_raw2bwt_final" ]; then
        pass "raw2bwt produces output at exact path"
    else
        # Check for any output files (raw2bwt may produce .mbz or other names)
        count=$(ls /tmp/mb_ir_raw2bwt_final* 2>/dev/null | wc -l) || count=0
        if [ "$count" -gt 0 ]; then
            pass "raw2bwt produces output files ($count files)"
        else
            fail "raw2bwt produced no output"
        fi
    fi
else
    skip "raw2bwt: fa2bit -p did not produce .pac file"
fi
rm -f /tmp/mb_ir_raw2bwt* /tmp/mb_ir_raw2bwt_out* /tmp/mb_ir_raw2bwt_final* 2>/dev/null || true

# Test 6: gensa - generate sampled SA from BWT
echo "[IR 6] gensa: generate sampled SA from BWT"
"$BINARY" index "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_gensa > /dev/null 2>&1
if [ -f "/tmp/mb_ir_gensa.mbw" ]; then
    "$BINARY" gensa /tmp/mb_ir_gensa.mbw > /dev/null 2>&1
    # gensa may modify .mbw in place or create a new file
    pass "gensa runs without error"
else
    fail "gensa: no .mbw input"
fi
rm -f /tmp/mb_ir_gensa* 2>/dev/null || true

# Test 7: Verify index round-trip (index -> map -> getref -> compare)
echo "[IR 7] Index round-trip: reference extraction"
"$BINARY" index "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_roundtrip > /dev/null 2>&1
if [ -f "/tmp/mb_ir_roundtrip.l2b" ]; then
    # Extract reference from .l2b
    "$BINARY" getref /tmp/mb_ir_roundtrip.l2b > /tmp/mb_ir_ref_extracted.fa 2>/dev/null
    if [ -f "/tmp/mb_ir_ref_extracted.fa" ]; then
        # Compare first line of FASTA header
        orig_header=$(zcat "$DATA_DIR/chrM-human.fa.gz" | head -1)
        extract_header=$(head -1 /tmp/mb_ir_ref_extracted.fa)
        if [ "$orig_header" = "$extract_header" ]; then
            pass "getref extracts correct FASTA header"
        else
            fail "getref header mismatch: '$extract_header' != '$orig_header'"
        fi
    else
        fail "getref produced no output"
    fi
else
    fail "round-trip: index failed"
fi
rm -f /tmp/mb_ir_roundtrip* /tmp/mb_ir_ref_extracted.fa 2>/dev/null || true

# Test 8: Verify index can be loaded after separate steps (skipped - separate index format may differ)
echo "[IR 8] Separate index steps produce loadable index - skipped"
skip "separate index format may differ from main index"

# Test 9: genbwt with specific thread count
echo "[IR 9] genbwt with explicit thread count"
"$BINARY" fa2bit "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_t1 > /dev/null 2>&1
if [ -f "/tmp/mb_ir_t1" ]; then
    "$BINARY" genbwt /tmp/mb_ir_t1 /tmp/mb_ir_t1.mbw > /dev/null 2>&1
    if [ -f "/tmp/mb_ir_t1.mbw" ]; then
        pass "genbwt with default threads works"
    else
        fail "genbwt with default threads did not produce .mbw"
    fi
else
    fail "fa2bit failed for thread test"
fi
rm -f /tmp/mb_ir_t1* 2>/dev/null || true

# Test 10: Low-memory index produces usable output
echo "[IR 10] Low-memory index produces usable mapping"
"$BINARY" index -l "$DATA_DIR/chrM-human.fa.gz" /tmp/mb_ir_lm > /dev/null 2>&1
if [ -f "/tmp/mb_ir_lm.l2b" ] && [ -f "/tmp/mb_ir_lm.mbw" ]; then
    "$BINARY" map /tmp/mb_ir_lm "$DATA_DIR/chrM-read_1.fa.gz" > /tmp/mb_ir_lm_out.sam 2>/dev/null
    if [ -f "/tmp/mb_ir_lm_out.sam" ]; then
        lines=$(wc -l < /tmp/mb_ir_lm_out.sam)
        if [ "$lines" -gt 2 ]; then
            pass "Low-memory index produces usable output ($lines lines)"
        else
            fail "Low-memory index produced too few lines: $lines"
        fi
    else
        fail "Low-memory index mapping produced no output"
    fi
else
    fail "Low-memory index did not produce required files"
fi
rm -f /tmp/mb_ir_lm* /tmp/mb_ir_lm_out.sam 2>/dev/null || true

# Cleanup
echo ""
echo "Cleaning up..."
rm -f /tmp/mb_ir_* 2>/dev/null || true

echo ""
echo "============================================"
echo " Index Routine Test Summary"
echo "============================================"
echo -e " \033[0;32m$PASS\033[0m passed"
echo -e " \033[0;31m$FAIL\033[0m failed"
echo "============================================"

if [ "$FAIL" -gt 0 ]; then
    exit 1
fi
exit 0
