"""1.0.1: exact float parsing of SMR tables and the explicit missing-identifier rule (nomibench.io)."""
import pandas as pd
import pytest

from nomibench.cli import _read
from nomibench.io import MISSING_ID_TOKENS, missing_id_mask, read_smr, read_table

# Fixed-notation P values as SMR writes them. pandas' default parser keeps only about 17 significant characters
# of such numbers, so several of these come back rounded or as 0; the first is a value from paper 1's table.
LONG = ["0.0000000000000000000000000000009043744", "0.000001234567890123456789", "0.0000000999999999999999999999",
        "0.00000000000000000000000000000000000000000000000000000000012345", "0.1541948", "1"]


def _smr(tmp_path, rows):
    p = tmp_path / "smr.csv"
    lines = ["probeID,Gene,p_SMR,p_HEIDI"] + [",".join(r) for r in rows]
    p.write_text("\n".join(lines) + "\n")
    return p


def test_exact_float_parsing(tmp_path):
    p = _smr(tmp_path, [(f"P{i:04d}", f"G{i}", v, "0.5") for i, v in enumerate(LONG)])
    df = read_smr(p)
    assert list(df.p_SMR) == [float(v) for v in LONG]          # bit for bit, as Python's float() reads them
    assert (df.p_SMR > 0).all()
    t = read_table(p)
    assert list(t.p_SMR) == [float(v) for v in LONG]


def test_tab_separated_smr_output(tmp_path):
    p = tmp_path / "chr1.smr"
    p.write_text("probeID\tp_SMR\tp_HEIDI\nGENE0001\t" + LONG[0] + "\tNA\n")
    df = read_smr(p)
    assert df.p_SMR.iloc[0] == float(LONG[0])
    assert pd.isna(df.p_HEIDI.iloc[0])                          # a missing HEIDI stays missing


def test_p_values_outside_unit_interval_raise(tmp_path):
    p = _smr(tmp_path, [("GENE0001", "A", "1.2", "0.5")])
    with pytest.raises(ValueError, match="outside"):
        read_smr(p)


@pytest.mark.parametrize("token", ["", "NA", "na", " N/A ", "NaN", "null", "NULL", "None", "-", "."])
def test_missing_identifier_tokens(tmp_path, token):
    p = tmp_path / "genes.csv"
    p.write_text(f"gene,Y\nGENE0001,1\n{token},0\nGENE0003,0\n")
    with pytest.raises(ValueError, match="missing identifier"):
        read_table(p, id_columns=("gene",))
    df = read_table(p, id_columns=("gene",), missing_ids="drop")
    assert list(df.gene) == ["GENE0001", "GENE0003"]
    assert df.attrs["dropped_missing_ids"] == {"gene": 1}


def test_identifiers_are_text_whatever_the_pandas_version(tmp_path):
    p = tmp_path / "genes.tsv"
    p.write_text("gene\tY\n0001\t1\n  NAA10  \t0\nNANOG\t0\n")
    df = read_table(p, id_columns=("gene",))
    assert list(df.gene) == ["0001", "NAA10", "NANOG"]          # leading zeros kept; symbols starting NA are genes
    assert all(type(v) is str for v in df.gene)


def test_duplicated_identifiers_raise(tmp_path):
    p = tmp_path / "genes.csv"
    p.write_text("gene,Y\nA,1\nA,0\n")
    with pytest.raises(ValueError, match="share an identifier"):
        read_table(p, id_columns=("gene",))
    assert len(read_table(p, id_columns=("gene",), unique=False)) == 2


def test_mask_and_tokens():
    assert missing_id_mask(["A", None, float("nan"), " na ", "N/A", "NAB2"]).tolist() == \
        [False, True, True, True, True, False]
    assert "" in MISSING_ID_TOKENS


def test_cli_reader(tmp_path):
    p = tmp_path / "genes.csv"
    p.write_text("gene,Y\nGENE0001,1\nNA,0\n")
    with pytest.raises(ValueError):
        _read(str(p), "gene")
    assert len(_read(str(p), "gene", "drop")) == 1
    assert len(_read(str(p))) == 2                              # no identifier column named: no rule applied
