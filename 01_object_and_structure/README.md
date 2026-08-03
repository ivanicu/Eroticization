# the object and its structure

**Rounds r01–r06, 6 of the project's 52.**

What the release actually is. The item-level data is not the 68 category ratings; it is 101 multiselect columns exploding to 1,332 options. Entry to every block is gated on a parent rating at P=0.99, which is undocumented and constrains every later design.

| round | directory |
|---|---|
| `r01` | [`r01_schema`](r01_schema) |
| `r02` | [`r02_column_kinds`](r02_column_kinds) |
| `r03` | [`r03_explode_multiselects`](r03_explode_multiselects) |
| `r04` | [`r04_gated_tree`](r04_gated_tree) |
| `r05` | [`r05_reliability_ceilings`](r05_reliability_ceilings) |
| `r06` | [`r06_option_census`](r06_option_census) |

Findings, intervals and caveats for every round live in the top-level [`README.md`](../README.md); this file is a table of contents, not a second account of the results — one home per fact.
