# BS-seq Ground-Truth Position Validation

## Purpose

Add ground-truth position validation to `test/test-bisulfite.sh` so that simulated
bisulfite-reads mapping results are compared against known-good expected positions,
focusing on the SAM `RNAME` (chromosome) and `POS` (leftmost position) fields.

## Problem

The current BS-seq test suite (`test-bisulfite.sh`) generates simulated reads from
known positions on the reference (chrM, ~16.5 kb) using deterministic formulas, then
validates that the mapper produces SAM output with correct headers and line counts.
However, it never verifies that the reads actually map to the positions they were
generated from — a critical correctness check.

## Solution

Add a new test group `[17] Ground-Truth Position Validation` to `test-bisulfite.sh`
that validates mapped read positions against the known generation positions.

## Design

### 1. New helper function

```bash
validate_bs_positions() {
    # Args: sam_file ref_len read_len strand step offset expected_count
    # Computes expected positions from the generation formula,
    # parses SAM POS/RNAME fields, compares against expected.
    # Reports per-read pass/fail, returns summary counts.
}
```

The function:
- Computes expected positions using the same formula as the generator:
  - Forward: `pos = (i * step + offset) % ref_len`
  - Reverse: `pos = (i * step + offset) % ref_len`
- Parses SAM data lines (non-header, tab-delimited)
- Matches read names to expected positions
- Compares RNAME against `chrM`
- Compares POS against expected (exact match, tolerance = 0)
- Counts unmapped reads (flag 4) separately
- Returns pass/fail counts via global variables

### 2. New test cases (5 tests)

| Test | Description | Reads Validated |
|------|-------------|-----------------|
| [17.1] Forward SE position check | Validate 100 C-to-T converted reads map to expected positions | 100 |
| [17.2] Reverse SE position check | Validate 50 G-to-A converted reads map to expected positions | 50 |
| [17.3] PE position check | Validate R1 and R2 positions for 50 simulated PE pairs | 100 |
| [17.4] Cross-mode consistency | Verify regular and --meth modes agree on chr/pos for overlapping reads | ~50 |
| [17.5] Batch position validation | Validate positions for large batch (200 reads) and mixed-strand batch | 200+ |

### 3. Position tolerance

**Exact match (±0 bp)**. The reference is small (chrM ~16.5 kb), reads are generated
from known positions, and bisulfite-converted reads should map back to their source
positions. Tolerance can be increased later if needed.

### 4. Edge cases

- **Unmapped reads** (flag 4): Counted as skipped, not failed.
- **Reverse strand**: SAM POS is the leftmost alignment coordinate (already adjusted
  by the mapper). No flag inversion needed for comparison.
- **Boundary reads**: Reads truncated by the mapper at reference boundaries are
  excluded from exact matching and reported as notes.

### 5. Implementation details

- Add `validate_bs_positions()` near the top of `test-bisulfite.sh` (after existing
  helper functions).
- Append test group `[17]` after the existing `[16]` group.
- Reuse existing simulated read files (no regeneration needed).
- Each test calls `validate_bs_positions()` with the appropriate parameters.
- Summary line includes ground-truth validation pass/fail counts.

### 6. Files changed

- `test/test-bisulfite.sh` — add helper function + test group (~120 lines)

## Success Criteria

- All 5 new tests pass on the CI runner.
- Total BS-seq test count increases from ~64 to ~69.
- Ground-truth validation catches intentional position mismatches (verified via
  a deliberate test with wrong expected position).
