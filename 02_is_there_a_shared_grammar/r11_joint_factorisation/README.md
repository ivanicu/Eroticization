# r11_joint_factorisation

Campaign [`02_is_there_a_shared_grammar`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
ITER 1. The transferable subspace is real but unnamed. Name it from the data, not from an LLM
(the corpus that would label these options also contains the kink taxonomy = leakage).

Method: one joint low-rank factorization across ALL blocks at once, on OBSERVED cells only
(the survey is a gated tree, so "not asked" != "disliked"; imputing 0 would conflate them).
Within-block double-centering first, so factors are profile SHAPE, not propensity.
Held-out cells give an honest score. Extreme loadings NAME each factor.
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
