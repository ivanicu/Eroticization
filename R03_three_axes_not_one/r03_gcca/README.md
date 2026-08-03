# R03·r03 — gcca

Part of [`R03_three_axes_not_one`](../README.md).

## What this sub-round asks

```
SVD on the concatenation is variance-weighted, so the largest block (qi=1, n=14,901) owns every
factor -- drop-top-block stability was 0.25-0.70. A SHARED coordinate is not the direction with
the most total variance, it is the direction best represented in EVERY block.
That is MAXVAR generalized CCA: eigenvectors of sum_b P_b, P_b = projector onto block b's column
space. Block size enters only through the projector's rank, not through variance.
Same attack as before: drop the biggest block and see if the coordinate survives.
```

## Result

Top-level [`README.md`](../../README.md). This file records the design.
