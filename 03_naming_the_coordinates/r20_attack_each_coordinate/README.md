# r20_attack_each_coordinate

Campaign [`03_naming_the_coordinates`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
Attack each GCCA coordinate with its OWN matched control: find every block contributing to its
top-25 loadings, DELETE THEM ALL, refit, and ask whether the coordinate is still there.
Drop-one-block was too weak -- coord 4 reads as self/other but its loadings are the 7-block
fluid family, which shares one option template. If the coordinate is a template artifact it
cannot survive deleting the template.
Also: is coord 4 just the POWER composite from iteration 0 re-derived?
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
