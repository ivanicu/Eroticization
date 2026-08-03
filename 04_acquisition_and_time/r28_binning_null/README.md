# r28_binning_null

Campaign [`04_acquisition_and_time`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
The killer confound for TEMPO: the onset bins are coarse (0-4, 5-6 ... 19-25, 26+). Two
categories whose population means are both ~15.5 will land in the SAME bin for many people
purely mechanically, while categories 3 years apart rarely share one. That alone makes residual
onset correlation decay with mean-onset distance -- with no developmental tempo at all.

Parametric bootstrap null: generate onsets from person effect + category effect + iid noise
(NO tempo structure, by construction), push them through the SAME binning and the SAME
missingness pattern, and run the identical TEMPO regression. If the null reproduces t=4.6 the
finding is binning. If it gives t~0, tempo is real.
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
