# Paper 1: frozen aggregate outputs

These are the outputs of paper 1's registered analyses that `reproduce_registered.py` compares against
(`--step19`, `--step20`, `--step25`). They were frozen when the analyses ran (steps 19 and 20, September 2026;
step 25, amendment A10, 3 October 2026) and are deposited byte for byte in the repository and its archived
releases (the Python distributions do not include them). `SHA256SUMS` lists every file:

```bash
cd repro/paper1/frozen && sha256sum -c SHA256SUMS
```

## What is here

| Folder | Files | Contents |
|---|---|---|
| `step19/` | 21 | Step 19 (E1–E4): result table (`gw19_results.csv`), propensity balance, the 1,000 E4 matched control sets (one row per set), the plasmode summaries and their replicate files (one row per replicate and test), D1 results and D1 plasmode replicates, and the run reports |
| `step20/real/` | 5 | Step 20 (E5, amendment A5): result table, decision table, run report, and the testability summaries (counts per exposure set) |
| `step20/cal/` | 55 | Step 20 calibration (amendment A6 and the as-written run): design constants and generator coefficients per exposure, plasmode summaries and replicate files, the NC_E5a and NC_E5c control sets, balance tables, manifests and run logs |
| `step25/` | 8 | Amendment A10: the E2-design plasmode (design and replicate files), the 1,000 E2-matched control sets and the 2,000 replicates of the 5-Mb genomic-block bootstrap |

Every file is an aggregate. A row is a result, a covariate's balance statistic, a simulation replicate, a
control set or a bootstrap replicate; no row describes a gene. The D1 replicate files carry `yhash`, a hash of the
simulated outcome vector, which the script uses to confirm that a regenerated replicate is the same draw.

## What is not here

- The three inputs: the U_gw covariate table, the SMR/HEIDI table and the E5 exposure-set table. They are gene
  level and derived from MetaBrain, whose terms do not allow redistribution; rebuild them as `../README.md`
  describes. Their frozen SHA-256 values are in `nomibench.registered.FROZEN_SHA256`.
- The freeze manifests that record local file paths, the step-25 `private/` folder (a per-gene table) and the
  figures.

## Screening

Before deposit, every file was checked for Ensembl gene and transcript identifiers, rsIDs, gene symbols,
identifier columns, local file paths, machine or user names and e-mail addresses; none was found.
`tests/test_frozen_deposit.py` repeats the identifier and checksum checks, and `tests/test_governance_cli.py`
scans the CSV, JSON and text files for gene identifiers and private paths.

## Provenance checks

- The 55 calibration files match the hashes in the step-20 freeze manifest (28 September 2026), and the
  testability summaries match the same manifest.
- The step-19 tables (`gw19_results.csv`, `gw19_plasmode.csv`, `gw19_balance.csv`, `gw19_e4_sets.csv`,
  `d1_plasmode.csv`, `d1_results.csv`) and the step-20 tables (`e5_results.csv`, `e5_decisions.csv`) match the
  hashes recorded in paper 1's results freeze.
- The step-25 files came from the step-25 output folder, which was hash-checked when it was written (3 October
  2026); the freeze records only the summary tables' hashes, so for these files the evidence is R10, which
  regenerates their first replicates and control sets from the registered seeds.

With these files and the three rebuilt inputs, `reproduce_registered.py --require-frozen` passes all 28 criteria
in paper 1's step-20 environment (`docs/VALIDATION.md`).
