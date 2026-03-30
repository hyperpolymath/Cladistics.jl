# TEST-NEEDS: Cladistics.jl

## Current State

| Category | Count | Details |
|----------|-------|---------|
| **Source modules** | 2 | 1,197 lines |
| **Test files** | 1 | 606 lines, 161 @test/@testset |
| **Benchmarks** | 0 | None |

## What's Missing

- [ ] **Performance**: No benchmarks for phylogenetic tree computation
- [ ] **Error handling**: No tests for malformed taxa, cyclic trees

## FLAGGED ISSUES
- **161 tests for 2 modules = 80.5 tests/module** -- excellent
- **0 benchmarks** -- tree algorithms should have scaling benchmarks

## Priority: P3 (LOW)

## FAKE-FUZZ ALERT

- `tests/fuzz/placeholder.txt` is a scorecard placeholder inherited from rsr-template-repo — it does NOT provide real fuzz testing
- Replace with an actual fuzz harness (see rsr-template-repo/tests/fuzz/README.adoc) or remove the file
- Priority: P2 — creates false impression of fuzz coverage
