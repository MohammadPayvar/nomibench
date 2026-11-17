"""The deposited frozen outputs of paper 1 (repro/paper1/frozen): checksums and the no-gene-level rule (1.0.1).

The deposit is in the repository and its archived releases, not in the source distribution, so these tests skip
when the files are absent.
"""
import glob
import hashlib
import os
import re

import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FROZEN = os.path.join(ROOT, "repro", "paper1", "frozen")
SUMS = os.path.join(FROZEN, "SHA256SUMS")
present = pytest.mark.skipif(not os.path.isfile(SUMS), reason="frozen outputs are not in the source distribution")

# column names that would identify genes, variants or probes
ID_COLUMN = re.compile(r"^(ensg|ensembl(_id)?|gene(_id|_name)?|symbol|probe(id|_id)?|snp|rsid|variant(_id)?|"
                       r"transcript(_id)?|hgnc(_id|_symbol)?|entrez(_id)?)$", re.I)
ID_VALUE = re.compile(r"ENS[GT]\d{6,}|\brs\d{3,}\b")


def _listed():
    with open(SUMS, encoding="utf-8") as fh:
        return [ln.split(None, 1) for ln in fh.read().splitlines() if ln.strip()]


@present
def test_checksums_cover_every_file_and_match():
    listed = _listed()
    names = {n.strip() for _, n in listed}
    on_disk = {os.path.relpath(p, FROZEN).replace(os.sep, "/")
               for p in glob.glob(os.path.join(FROZEN, "**", "*"), recursive=True) if os.path.isfile(p)}
    assert on_disk - {"README.md", "SHA256SUMS"} == names
    for h, n in listed:
        with open(os.path.join(FROZEN, n.strip()), "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() == h, n


@present
def test_no_gene_level_columns_or_identifiers():
    for _, n in _listed():
        p = os.path.join(FROZEN, n.strip())
        with open(p, encoding="utf-8", errors="ignore") as fh:
            assert not ID_VALUE.search(fh.read()), n
        if n.endswith(".csv"):
            cols = pd.read_csv(p, nrows=0).columns
            assert not [c for c in cols if ID_COLUMN.match(str(c).strip())], n
