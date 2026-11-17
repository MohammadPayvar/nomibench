"""1.0.1: the multi-start rule for Firth fits (paper 1, amendment A15).

``firth_bimodal.csv`` is a simulated sparse design (800 genes, 8 events, events concentrated at one end of the
propensity score) whose penalised likelihood has two strict local maxima: penalised log-likelihood -21.153461 at
ln OR 1.193055, which a single start from zero finds, and -17.441079 at ln OR 1.424082. Both were confirmed with
an independent optimiser (scipy BFGS on the penalised likelihood; the numerical Hessian is negative definite at
each point).

The rule reduces, but does not remove, the risk of reporting a lower maximum: on this fixture 30 perturbed starts
reach the higher maximum for 18 of 60 seeds (0-59), which is why the tests fix the seed. A fit that found more than
one stationary point says so (``n_stationary``, ``zero_start_best``).
"""
import os

import numpy as np
import pandas as pd
import pytest

from conftest import DATA
from nomibench.estimators import g1, g1_design
from nomibench.firth import _select_start, firth_fit, firth_fit_multistart, firth_lrt


@pytest.fixture(scope="module")
def bimodal():
    d = pd.read_csv(os.path.join(DATA, "firth_bimodal.csv"), float_precision="round_trip")
    return g1_design(d.lps.values, d.e.values), d.y.values, d


def test_unimodal_returns_the_zero_start_fit_bit_for_bit():
    d = pd.read_csv(os.path.join(DATA, "firth_moderate.csv"), float_precision="round_trip")
    X = np.c_[np.ones(len(d)), d[["x1", "x2", "x3"]].values]
    f0 = firth_fit(X, d.y.values)
    fm = firth_fit_multistart(X, d.y.values, seed=1)
    assert np.array_equal(fm.beta, f0.beta) and fm.loglik == f0.loglik
    assert fm.zero_start_best and fm.n_stationary == 1 and fm.n_starts == 31 and not fm.failed


def test_bimodal_keeps_the_higher_maximum(bimodal):
    X, y, _ = bimodal
    f0 = firth_fit(X, y)
    fm = firth_fit_multistart(X, y, seed=24)
    assert f0.loglik == pytest.approx(-21.153461, abs=1e-6) and f0.beta[1] == pytest.approx(1.193055, abs=1e-6)
    assert fm.loglik == pytest.approx(-17.441079, abs=1e-6) and fm.beta[1] == pytest.approx(1.424082, abs=1e-6)
    assert fm.n_stationary == 2 and not fm.zero_start_best and fm.max_abs_score < 1e-6


def test_g1_default_is_multistart_and_starts_zero_is_registered(bimodal):
    _, y, d = bimodal
    single = g1(y, d.e.values, d.lps.values, starts=0)
    multi = g1(y, d.e.values, d.lps.values, start_seed=24)
    assert single.log_or == pytest.approx(1.193055, abs=1e-6) and single.starts == 0
    assert multi.log_or == pytest.approx(1.424082, abs=1e-6)
    assert multi.starts == 30 and multi.n_stationary == 2 and not multi.zero_start_best
    lo, hi = multi.profile[0.95]
    assert lo < multi.log_or < hi


def test_lrt_zero_start_never_lowers_the_restricted_maximum(bimodal):
    X, y, _ = bimodal
    fm = firth_fit_multistart(X, y, seed=24)
    a = firth_lrt(X, y, 1, full=fm)
    b = firth_lrt(X, y, 1, full=fm, zero_start=True)
    assert b.stat <= a.stat + 1e-9


def test_seeded_and_deterministic(bimodal):
    X, y, _ = bimodal
    a, b = firth_fit_multistart(X, y, seed=5), firth_fit_multistart(X, y, seed=5)
    assert np.array_equal(a.beta, b.beta)


def test_selection_rule():
    # (penalised log-likelihood, stationary); the zero start is first
    assert _select_start([(-10.0, True), (-9.0, True)]) == (1, 2, False)        # higher maximum wins
    assert _select_start([(-10.0, True), (-10.0 + 5e-7, True)]) == (0, 1, False)  # a tie goes to the zero start
    assert _select_start([(-10.0, True), (-9.0, False)]) == (0, 1, False)       # non-stationary starts ignored
    assert _select_start([(-10.0, False), (-11.0, True)]) == (1, 1, False)      # zero start not stationary
    assert _select_start([(-10.0, False), (-9.0, False)]) == (0, 0, True)       # no stationary point: failed
