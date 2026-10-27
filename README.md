# GPCRome basal cAMP screen

Code, processed data and figures for the basal cAMP GPCRome screen in **Howard et al.
2026** — a pooled measurement of **ligand-independent (basal) Gs/cAMP signalling across
946 barcoded human GPCR constructs**.

A barcoded GPCRome library is integrated into a HEK293T CRE–eGFP reporter line, grown
without added agonist, and sorted into four bins of reporter fluorescence. Each construct's
basal activity is the **weighted mean bin** of its barcode across those bins, rescaled to
0–1. Higher score = more cAMP with no agonist present.

**Correspondence:** Matthew K. Howard, <matthewkarlhoward@gmail.com>

**This repository is scoped to Figure 1 and Supplemental Figure 1.** It contains the
pipeline, the scoring, the annotation tables those two figures read, and the scripts that
drew them — and nothing else. The exploratory analyses around the screen, and the deep
mutational scan of GPR161 that follows from it, are not here.

```
pipeline/     FASTQ -> per-barcode counts (fastp, cutadapt, count_barcodes.py)
scripts/      five table builders, three figure scripts, one style module
data/         counts, sort QC, six derived tables, library and reference annotation
figures/      the panels of Figure 1 and Supplemental Figure 1, as .pdf and .png
docs/         methods and caveats, sequencing QC, panel -> script map
```

- **[`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md)** — which script drew which panel.
- **[`docs/METHODS.md`](docs/METHODS.md)** — how the scores are computed, which sort day
  was excluded and why, and the three confounds to know about before reusing them. Read it
  before quoting any number from here.

Cartoons (Figure 1a) were drawn by hand and the FACS gating plots (Supplemental Figure 1a)
were exported from FlowJo; neither is produced here.

---

## Reproducing the figures

Raw FASTQs are deposited separately (see the manuscript's data-availability statement);
`data/counts/` is the output of `pipeline/` and is the entry point for everything else.

```bash
pip install -r requirements.txt
bash run_all.sh          # or: PY=python3 bash run_all.sh
```

The scripts are run through `$PY`, which defaults to `python`; set it if your
interpreter is `python3`.

The **five derived tables** in `data/processed/` regenerate byte-identically from
`data/counts/` plus the reference data under `data/external/`, `data/inhouse/` and
`data/library/`. The sixth file there, `literature_constitutive_activity.csv`, is a
hand-curated **input** that nothing regenerates — do not delete it.

| Step | Script | Writes |
|---|---|---|
| 1 | `scripts/score_basal.py` | `basal_scores_per_barcode.csv` |
| 2 | `scripts/build_coupling.py` | `coupling_annotation.csv` |
| 3 | `scripts/build_ranked_tables.py` | `canonical_scores_ranked.csv` |
| 4 | `scripts/build_construct_features.py` | `receptor_construct_features.csv` |
| 5 | `scripts/build_surface_join.py` | `basal_vs_tedman_pme.csv` |
| 6 | `scripts/plot_figure1.py` | **Figure 1b–e** |
| 7 | `scripts/plot_qc_panels.py` | **Supplemental Figure 1b–e** |
| 8 | `scripts/plot_surface_vs_basal.py` | **Supplemental Figure 1f** |

Steps 1–5 build tables and draw nothing; steps 6–8 draw and write nothing. Besides the
published panels, step 6 also emits `F1_main_canonical` (a canonical-only variant of panel
b that was not used) and step 7 `Q_qc_strip` (panels b–e pre-arranged as one strip).

### From raw reads

```bash
bash pipeline/process_one.sh <sample>_R1_001.fastq.gz   # one sample
bash pipeline/run_paired.sh <dir-of-per-sample-folders> # the whole run
```

fastp merges R1/R2 into a consensus read, cutadapt pulls the 10 bp barcode out with linked
primer anchors, and `count_barcodes.py` matches it to the 946-barcode reference allowing
≤1 mismatch. The closest pair of reference barcodes is 4 mismatches apart, so 1 mismatch is
unambiguous. See [`docs/SEQUENCING_QC.md`](docs/SEQUENCING_QC.md).

---

## Citing this repository

`CITATION.cff` carries the citation metadata; GitHub renders it as a **Cite this
repository** button. `.zenodo.json` supplies the metadata Zenodo uses when it archives a
release.

Tagged releases are archived on [Zenodo](https://zenodo.org) through its GitHub
integration, which mints a DOI per release plus a concept DOI that always resolves to the
latest one. **Cite the concept DOI** unless you need to pin an exact version.

> Not yet archived. To set it up: enable this repository under Zenodo → *GitHub*, then
> publish a GitHub release — Zenodo archives the tag and mints the DOI. Afterwards, add the
> DOI badge here and the `doi:` field to `CITATION.cff`, and note that Zenodo picks up
> `.zenodo.json` at the moment the release is published, so any metadata change has to be
> committed before the tag.

## Licence

Copyright (c) 2026 Matthew K. Howard.

Code in this repository is released under the MIT licence (`LICENSE`); the data generated
here is released under CC BY 4.0. Third-party reference data under `data/external/` is not
covered by either and keeps its own licence and attribution — see
[`data/external/README.md`](data/external/README.md).
