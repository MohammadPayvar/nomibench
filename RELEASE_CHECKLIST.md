# Release checklist (v1.0.1)

These steps publish the package under your name and accounts, so only you can do them. Version 1.0.0 was never
tagged or archived; 1.0.1 is the first archived release.

## 1. Metadata (done for 1.0.1)

| File | State |
|---|---|
| `pyproject.toml` | version 1.0.1; author; project URLs |
| `CITATION.cff` | version 1.0.1; author with ORCID and affiliation; `repository-code`; `date-released` 2026-11-17 (change it if you tag on another day); `doi` to add once Zenodo has minted it (step 3) |
| `.zenodo.json` | version 1.0.1; creator with ORCID and affiliation; add the paper as a `related_identifiers` entry once it has a DOI |
| `README.md` | the paper (as a manuscript); the Zenodo DOI goes under "Citation" once minted |
| `LICENSE` | "the nomibench authors"; change it if the copyright holder should be named |

Check that no placeholder is left:

```bash
grep -rn "TO FILL\|OWNER\|FAMILY NAME\|GIVEN NAMES\|to be added" --exclude=RELEASE_CHECKLIST.md --exclude-dir=.git .
```

## 2. GitHub

The `tests` workflow runs the fast suite on Linux, Windows and macOS with Python 3.10–3.13, the slow simulation
checks on Linux, and the linter. All three jobs must pass on `main` before you tag the release.

## 3. Get a DOI from Zenodo

Choose one route. Route B archives the whole repository, including `repro/paper1/frozen/` (the frozen aggregate
outputs of paper 1); the Python distributions in route A do not contain those files.

**A. Reserved DOI (use this if the paper must cite the exact version before it is published).**
1. On Zenodo, start a new upload. For "Do you already have a DOI?" answer No and click "Get a DOI now!".
   The DOI is reserved but not yet registered.
2. Put the DOI into `CITATION.cff` (`doi:`), the README's "Citation" section and the manuscript's code-availability
   statement. Commit, push and tag `v1.0.1`.
3. Build the distributions (step 4) and a full archive of the tag
   (`git archive --format=zip -o nomibench-1.0.1-repository.zip v1.0.1`). Upload all three to the Zenodo draft,
   copy the metadata from `.zenodo.json`, and publish.

**B. GitHub integration.**
1. In Zenodo's GitHub settings, switch the repository on.
2. Create the GitHub release `v1.0.1` from `main`. Zenodo archives the tag and mints a version DOI and a concept
   DOI. It reads its metadata from `.zenodo.json`; when that file exists, Zenodo ignores `CITATION.cff`, although
   GitHub still shows it in the citation box.
3. The archived v1.0.1 cannot contain its own DOI. Add it to the README and `CITATION.cff` in a follow-up commit,
   and put the version DOI in the manuscript.

## 4. Build and check the distributions

```bash
python -m pip install --upgrade build twine
python -m build
twine check dist/*
python -m venv /tmp/nb && /tmp/nb/bin/pip install "dist/nomibench-1.0.1.tar.gz[test]"
tar -xzf dist/nomibench-1.0.1.tar.gz -C /tmp && cd /tmp/nomibench-1.0.1 && /tmp/nb/bin/python -m pytest -q -m "not slow"
```

Confirm that the code is the validated code: the source-tree hash in `docs/VALIDATION.md` (§4, 1.0.1) must not
change when only metadata and documentation were edited (the command is given there).

## 5. PyPI (optional)

`twine upload dist/*` with your own PyPI API token. The version on PyPI must be the same files as on Zenodo.

## 6. Data governance

- The repository contains no MetaBrain-derived gene-level values. `tests/test_governance_cli.py` scans for
  Ensembl gene IDs and private paths, and `tests/test_frozen_deposit.py` checks that the deposited frozen outputs
  match their checksums and carry no gene, variant or probe identifiers.
- `repro/paper1/frozen/` holds paper 1's frozen aggregate outputs only: result tables, simulation replicates,
  control sets and bootstrap replicates. Never add the three inputs (the U_gw covariate table, the SMR table, the
  E5 exposure-set table), a reproduction's own output folder, any `private/` folder or a freeze manifest with local
  paths; `.gitignore` blocks the inputs and `private/`.
