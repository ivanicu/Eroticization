# r26_onset_beyond_preference

Campaign [`04_acquisition_and_time`](../README.md) · run with `python run.py` from anywhere.

## What this round asks

```
Two remaining alternatives to 'acquired as a unit':
 (T) topical near-synonymy -- secretions/dirtiness, vore/bestiality are the same thing twice.
     Control: partial out lexical overlap of the category NAMES (no LLM, no embedding).
 (R) onset structure is redundant with preference structure and adds nothing.
     Control: residualise the onset-similarity matrix on the preference-similarity matrix and
     test whether the RESIDUAL still has non-random structure. Redundant => nothing left.
Already noticed: secretions+abnormal-body has r_onset=+0.149 with r_pref=-0.019, so the two
matrices are demonstrably not collinear. Quantify it.
```

## Result

The number, its scope and its controls are in the top-level [`README.md`](../../README.md). This file records the design; the ledger records the finding.
