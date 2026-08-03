# r19_gcca

Campaign [`03_naming_the_coordinates`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
SVD on the concatenation is variance-weighted, so the largest block (qi=1, n=14,901) owns every
factor -- drop-top-block stability was 0.25-0.70. A SHARED coordinate is not the direction with
the most total variance, it is the direction best represented in EVERY block.
That is MAXVAR generalized CCA: eigenvectors of sum_b P_b, P_b = projector onto block b's column
space. Block size enters only through the projector's rank, not through variance.
Same attack as before: drop the biggest block and see if the coordinate survives.
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
