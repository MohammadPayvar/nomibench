"""Reading gene tables and SMR output: exact numbers and an explicit rule for missing identifiers (1.0.1).

**Exact float parsing.** pandas' default C parser converts long fixed-notation numbers to the nearest double only
up to about 17 significant characters; beyond that it truncates. SMR writes P values in fixed notation, so a value
such as ``0.0000000000000000000000000000009043744`` can be read as 0 or with a wrong mantissa. In paper 1's
registered SMR table (11,224 probes) this misread 24 ``p_SMR`` values (4 as 0) and 13 ``p_HEIDI`` values (2 as 0),
and changed no gene's outcome under the BH or the Bonferroni rule. Every read here uses
``float_precision="round_trip"``, which gives the same double as Python's ``float()``.

**Missing identifiers.** Identifier columns are read as text with no automatic NA conversion, so the result does
not depend on how a pandas version treats ``NaN`` under ``astype(str)``. A value is a missing identifier if, after
stripping surrounding white space, it is empty or one of :data:`MISSING_ID_TOKENS` (case-insensitive). Rows with a
missing identifier either raise (``missing_ids="error"``, the default) or are dropped and counted
(``missing_ids="drop"``; the counts are in ``df.attrs["dropped_missing_ids"]``). Duplicated identifiers raise
unless ``unique=False``.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

__all__ = ["MISSING_ID_TOKENS", "read_table", "read_smr", "missing_id_mask"]

MISSING_ID_TOKENS = frozenset({"", "na", "n/a", "nan", "null", "none", "-", "."})
"""Identifier values treated as missing (compared after stripping white space, case-insensitively)."""


def _sep(path, sep):
    if sep is not None:
        return sep
    p = str(path).lower()
    if p.endswith(".gz"):
        p = p[:-3]
    return "\t" if p.endswith((".tsv", ".txt", ".smr")) else ","


def missing_id_mask(values) -> np.ndarray:
    """Boolean mask of missing identifiers under the rule above (``None`` and float NaN count as missing)."""
    s = pd.Series(values, dtype=object)
    out = np.zeros(len(s), bool)
    for i, v in enumerate(s):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            out[i] = True
        else:
            out[i] = str(v).strip().lower() in MISSING_ID_TOKENS
    return out


def read_table(path, id_columns=(), missing_ids="error", unique=True, sep=None, **kwargs) -> pd.DataFrame:
    """Read a CSV/TSV (optionally gzipped) or Parquet table with exact float parsing.

    Parameters
    ----------
    path : file path; ``.tsv``, ``.txt`` and ``.smr`` (optionally ``.gz``) are tab-separated, others comma-separated
        unless ``sep`` is given; ``.parquet`` needs pyarrow.
    id_columns : identifier columns, read as text and checked with the missing-identifier rule.
    missing_ids : ``"error"`` (default) or ``"drop"``.
    unique : raise if the first identifier column has duplicated values.
    kwargs : passed to :func:`pandas.read_csv` (not ``float_precision``, ``dtype`` for identifier columns or
        ``keep_default_na``, which this function sets).
    """
    if missing_ids not in ("error", "drop"):
        raise ValueError("missing_ids must be 'error' or 'drop'")
    id_columns = list(id_columns)
    for k in ("float_precision", "keep_default_na"):
        if k in kwargs:
            raise ValueError(f"read_table sets {k} itself")
    if str(path).lower().endswith(".parquet"):
        try:
            df = pd.read_parquet(path)
        except ImportError as err:
            raise ValueError("reading .parquet needs pyarrow: pip install 'nomibench[parquet]'") from err
        for c in id_columns:
            df[c] = df[c].astype(object).where(df[c].notna(), None)
            df[c] = [None if v is None else str(v) for v in df[c]]
    else:
        s = _sep(path, sep)
        df = pd.read_csv(path, sep=s, float_precision="round_trip", **kwargs)
        if id_columns:
            usecols = [c for c in id_columns if c in df.columns]
            missing_cols = sorted(set(id_columns) - set(usecols))
            if missing_cols:
                raise ValueError(f"{os.path.basename(str(path))}: identifier column(s) not found: {missing_cols}")
            raw = pd.read_csv(path, sep=s, usecols=usecols, dtype=str, keep_default_na=False,
                              **{k: v for k, v in kwargs.items() if k not in ("usecols", "dtype")})
            for c in usecols:
                df[c] = raw[c].astype(object).values
    df.attrs["dropped_missing_ids"] = {}
    for c in id_columns:
        if c not in df.columns:
            raise ValueError(f"identifier column {c!r} not found")
        m = missing_id_mask(df[c].values)
        if m.any():
            if missing_ids == "error":
                raise ValueError(f"{int(m.sum())} row(s) have a missing identifier in {c!r} "
                                 f"(empty or one of {sorted(MISSING_ID_TOKENS - {''})}); "
                                 "pass missing_ids='drop' to drop them")
            df.attrs["dropped_missing_ids"][c] = int(m.sum())
            df = df.loc[~m].reset_index(drop=True)
        df[c] = [str(v).strip() for v in df[c]]
    if id_columns and unique:
        c = id_columns[0]
        dup = df[c].duplicated(keep=False)
        if dup.any():
            raise ValueError(f"{int(dup.sum())} row(s) share an identifier in {c!r}, e.g. "
                             f"{df.loc[dup, c].iloc[0]!r}; pass unique=False if that is intended")
    return df


def read_smr(path, id_column="probeID", missing_ids="error") -> pd.DataFrame:
    """Read SMR output (one or many chromosomes concatenated) exactly.

    The probe identifier is checked with the missing-identifier rule and must be unique; ``p_SMR`` and
    ``p_HEIDI`` (when present) must lie in [0, 1] or be missing.
    """
    df = read_table(path, id_columns=(id_column,), missing_ids=missing_ids)
    for c in ("p_SMR", "p_HEIDI"):
        if c in df.columns:
            v = pd.to_numeric(df[c], errors="raise").astype(float).values
            bad = ~np.isnan(v) & ((v < 0) | (v > 1))
            if bad.any():
                raise ValueError(f"{int(bad.sum())} value(s) of {c} lie outside [0, 1]")
            df[c] = v
    return df
