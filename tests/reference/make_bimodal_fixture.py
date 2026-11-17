"""Regenerate ``tests/data/firth_bimodal.csv`` (nomibench 1.0.1).

A sparse simulated design (no real genes) whose Firth penalised likelihood has two strict local maxima. It came
from a search over random sparse designs with events concentrated at the high end of the propensity score; seed 24
gave the case below. Run from the repository root:

    python tests/reference/make_bimodal_fixture.py [OUT]

The design and the two maxima were confirmed with scipy's BFGS on the penalised likelihood (the numerical Hessian
of the negative penalised likelihood is positive definite at each point); see ``tests/test_firth_multistart.py``.
"""
import sys

import numpy as np
import pandas as pd

SEED = 24


def make(seed=SEED):
    rng = np.random.default_rng(seed)
    n = rng.choice([400, 800, 1500])
    lps = rng.normal(size=n) * rng.choice([1, 2, 3])
    e = rng.binomial(1, 1 / (1 + np.exp(-(lps - 1))))
    k = rng.choice([4, 6, 8, 12])
    order = np.argsort(lps)
    y = np.zeros(n)
    pool = order[-rng.integers(k, 6 * k):]
    y[rng.choice(pool, size=min(k, len(pool)), replace=False)] = 1
    if rng.random() < 0.5:
        y[rng.choice(n, 2, replace=False)] = 1
    return pd.DataFrame(dict(lps=lps, e=e, y=y))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "tests/data/firth_bimodal.csv"
    make().to_csv(out, index=False, float_format="%.17g")
