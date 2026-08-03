# r27_unit_vs_tempo

Campaign [`04_acquisition_and_time`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
ITER 4. Two live readings of the onset structure, and they predict DIFFERENT things.

  UNIT ACQUISITION : categories sharing a GCCA coordinate were acquired together.
      predictor = coordinate-loading similarity between categories.
  DEVELOPMENTAL TEMPO : a person's maturational schedule shifts categories that the POPULATION
      places at similar ages together, regardless of any coordinate.
      predictor = -|mean_onset_i - mean_onset_j|  (population arrival-time distance)

Both regressed on within-person onset similarity, with preference similarity controlled.
Prediction matrix: unit -> coordinate term survives, tempo term dies. tempo -> the reverse.
Both -> decompose. Neither -> the iter-3 residual structure is something else again.
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
