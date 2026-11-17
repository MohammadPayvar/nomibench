# Changelog

## 1.0.1 (2026-11-17)

A patch release. The registered numbers of paper 1 are unchanged: in paper 1's step-20 environment the reproduction
gives the same 28 passing criteria, line for line, under 1.0.0 and 1.0.1 (`docs/VALIDATION.md`, §0).

**Exact reading of tables (`nomibench.io`, new)**
- `read_table` and `read_smr` parse floats exactly (`float_precision="round_trip"`). pandas' default parser
  rounds long fixed-notation numbers, and SMR writes P values that way: in paper 1's SMR table it misread 24 of
  11,224 `p_SMR` values (4 as 0; all below 6×10⁻¹¹) and 13 `p_HEIDI` values (2 as 0). No gene's outcome changed
  under the BH or the Bonferroni rule. `reproduce_registered.py` and `build_covariates.py` now use these readers;
  the command line uses `read_table`.
- **Missing identifiers.** Identifier columns are read as text without automatic NA conversion. A value that is
  empty or one of `NA`, `N/A`, `NaN`, `null`, `None`, `-` or `.` (case-insensitive, after stripping white space)
  is missing. Such rows raise by default or are dropped and counted (`missing_ids="drop"`). Duplicated
  identifiers raise unless `unique=False`. `read_smr` also checks that `p_SMR` and `p_HEIDI` lie in [0, 1].
- Command line: `--id-column` and `--missing-ids`. `.txt` and `.smr` files (optionally gzipped) are now read as
  tab-separated, like `.tsv`.

**Multi-start Firth fits**
- `firth_fit_multistart`: the fit from zero, plus 30 starts perturbed from it by Gaussian noise of 0.5, 2 and
  5 Fisher standard errors (10 at each size); a start counts if it converges to a stationary point (largest modified score below 1e-6); the estimate
  is the stationary point with the highest penalised likelihood, with ties within 1e-6 going to the zero start.
  This is the rule paper 1 registered in amendment A15. With a unimodal penalised likelihood, the usual case,
  the result is the zero-start fit bit for bit.
- `g1` and `evaluate` use it by default (`starts=30`, `start_seed=0`) and report `starts`, `n_stationary` and
  `zero_start_best`. The penalised LRT uses the better of two restricted fits (from the estimate and from zero).
  `starts=0` gives the single start of paper 1's steps 19–25 (`nomibench.registered.STARTS`), which the
  reproduction uses. Command line: `evaluate --starts`.
- Plasmode replicates, matched control sets and the block bootstrap keep the single zero start, as registered.
- The rule lowers, but does not remove, the risk of reporting a lower local maximum: on the bimodal test fixture,
  30 starts reach the higher maximum for 18 of 60 seeds.

**Tests and repository**
- New tests: exact parsing and the missing-identifier rule (`tests/test_io.py`); the multi-start rule on a
  bimodal and a unimodal design and the selection rule itself (`tests/test_firth_multistart.py`, with
  `tests/data/firth_bimodal.csv` and its generator `tests/reference/make_bimodal_fixture.py`); the deposited
  frozen outputs (`tests/test_frozen_deposit.py`).
- `repro/paper1/frozen/`: paper 1's frozen aggregate outputs of steps 19, 20 and 25 (89 files, no gene-level
  values), with `SHA256SUMS`, so the reproduction's comparison files are public.
- The linter's rule set is stated explicitly (ruff 0.16 enlarged its defaults).
- Metadata: author, ORCID, affiliation and the repository's new address.

## 1.0.0 (2026-10-06)

First release, accompanying paper 1. The estimators are the registered analysis code, ported without changing
their arithmetic, and reproduce the registered numbers (`docs/VALIDATION.md`).

Changes from the registered scripts:

**Convergence and failures**
- `firth_fit` reports `converged`, `n_iter` and the largest modified score, and warns when its iteration cap
  binds. The default cap is 1,000; the registered analyses used 100 (`nomibench.registered.MAXIT`).
- `firth_profile_ci` judges convergence only at the restricted fit nearest each bound, and reports `saturated`
  when a fitted probability there is below 1e-15 or above 1 - 1e-15 (complete separation; treat the bound as
  unbounded). An unbracketed bound is NaN, as before.
- `cfl_fit` reports convergence (cap unchanged at 200).
- Every plasmode estimator is wrapped so that a failed fit returns NaN for that replicate; NaN counts as a
  non-rejection, as registered. Summaries report the failure rate, mark a scenario with half or more failed
  fits "FAILED FITS", and give the profile false-exclusion rate both ways (failure as non-exclusion and as
  exclusion).
- `fit_ps` raises if the propensity model separates the exposure; its c-statistic is NaN, not an error, when a
  group is empty.

**Decisions**
- The plasmode's equivalence check stores the event itself (Wald upper bound below the margin). The registered
  step-19 code stored it as a 0/1 pseudo-P, which its report read the wrong way round (paper 1, step-28 erratum).
- `gate_rule(lo, power, calibrated)` applies the registered evaluability condition (calibrated, and power at
  OR 2.7 of at least 0.80) and returns "uninformative" otherwise. `coprimary_reading` accepts the gate in any
  case and rejects unknown values, so a typo cannot silently skip the downgrade.

**Inputs and reproducibility**
- One-call entry points (`evaluate`, `calibrate`, `matched_controls`), result objects and a command line.
- Gene tables are validated (0/1 coding, completeness, finite covariates) before any fit; the command line prints
  a one-line error instead of a traceback.
- Every simulation, set draw, bootstrap and fold split requires a seed. A `SeedSequence` passed as a root is
  copied, so reusing it reproduces the same streams.
- `calibrate` refuses N2, H and W1 without locus columns, reports `calibrated = None` unless both N0 and N1 were
  simulated with G1, and also requires at most 5% failed G1 fits under each null (`max_failed`; not part of
  the registered rule). It builds every requested estimator's inputs (G1b quintiles, G1c and H designs, W1
  strata). Locus identifiers of any type are factorised in sorted order, so integer ids 0..K-1 keep their values.
- `evaluate` reports H or W1 as "not estimable", with NaN, when the data cannot identify them.
- The `linear` argument defaults to no linear covariates everywhere; paper 1's choice is
  `nomibench.registered.TIER_P_LINEAR`.
- The plasmode's Wald equivalence check uses the exact normal quantile by default (`PlasmodeDesign.z_equiv`);
  step 19 used 1.959964 (`nomibench.registered.Z_STEP19`).
- The bootstrap's percentile interval uses exactly the 0.025 and 0.975 quantiles at level 0.95, as registered.
- D1's gradient boosting runs single-threaded, so results do not depend on the core count.
