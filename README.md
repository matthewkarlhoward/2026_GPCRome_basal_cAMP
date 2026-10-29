# GPCRome basal cAMP screen

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23682602.svg)](https://doi.org/10.5281/zenodo.23682602)

Code, processed data and figures for the basal cAMP GPCRome screen in **Howard et al.
2026** — a pooled measurement of **ligand-independent (basal) Gs/cAMP signalling across
946 barcoded human GPCR constructs**.

A barcoded GPCRome library is integrated into a HEK293T CRE–eGFP reporter line, grown
without added agonist, and sorted into four bins of reporter fluorescence. Each construct's
basal activity is the **weighted mean bin** of its barcode across those bins, rescaled to
0–1. Higher score = more cAMP with no agonist present.

**Correspondence:** Matthew K. Howard, <matthewkarlhoward@gmail.com>

```
pipeline/     FASTQ -> per-barcode counts (fastp, cutadapt, count_barcodes.py)
scripts/      six table builders, four figure scripts, a supplemental-table
              writer, one style module
data/         counts, sort QC, seven processed tables, library and reference annotation
figures/      the panels of Figure 1 and Supplemental Figures 1-2, as .pdf and .png
tables/       Supplemental Table 1
docs/         methods and caveats, sequencing QC, panel -> script map
```

- **[`docs/FIGURE_MAP.md`](docs/FIGURE_MAP.md)** — which script drew which panel.
- **[`docs/METHODS.md`](docs/METHODS.md)** — how the scores are computed, where the
  coupling annotation comes from, and what to know before reusing them. Read it before
  quoting any number from here.

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

The **six derived tables** in `data/processed/` regenerate byte-identically from
`data/counts/` plus the reference data under `data/external/`, `data/inhouse/` and
`data/library/`. The seventh file there, `literature_constitutive_activity.csv`, is a
hand-curated **input** that nothing regenerates — do not delete it. Nor delete
`data/inhouse/legacy_coupling.csv`: it is superseded for most receptors but is still the
last-resort coupling source for 42 of the 901 scored barcodes, so removing it changes
Figure 1b, Figure 1d and Supplemental Table 1 (see `docs/METHODS.md`).

| Step | Script | Writes |
|---|---|---|
| 1 | `scripts/score_basal.py` | `basal_scores_per_barcode.csv` |
| 2 | `scripts/build_coupling.py` | `coupling_annotation.csv` |
| 3 | `scripts/build_ranked_tables.py` | `canonical_scores_ranked.csv` |
| 4 | `scripts/build_construct_features.py` | `receptor_construct_features.csv` |
| 5 | `scripts/build_surface_join.py` | `basal_vs_tedman_pme.csv` |
| 6 | `scripts/build_isoform_table.py` | `alt_vs_canonical_isoform.csv` |
| 7 | `scripts/plot_figure1.py` | **Figure 1b–e** |
| 8 | `scripts/plot_qc_panels.py` | **Supplemental Figure 1b–e** |
| 9 | `scripts/plot_surface_vs_basal.py` | **Supplemental Figure 1f** |
| 10 | `scripts/plot_isoform_vignettes.py` | **Supplemental Figure 2** |
| 11 | `scripts/make_supp_table.py` | **Supplemental Table 1** |

Steps 1–6 build tables and draw nothing; steps 7–10 draw and write nothing; step 11
writes `tables/supp_table1_GPCRome_basal.csv` and draws nothing.

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
latest one.

| | DOI |
|---|---|
| **Cite this** — all versions, always resolves to the latest | [10.5281/zenodo.23682602](https://doi.org/10.5281/zenodo.23682602) |
| latest version only, pinned to that exact snapshot | listed on the Zenodo record under *Versions* |

Every earlier version keeps its own DOI and stays resolvable; they are listed on the
Zenodo record under *Versions*.

Cite the concept DOI unless you need to pin an exact version: a later release moves it to
the newer snapshot, so a DOI printed in the manuscript keeps resolving to current code.

To cut a new version, publish another GitHub release. Zenodo reads `.zenodo.json` **from
the tagged commit at the moment the release is published**, so any metadata change has to
be committed before the tag — editing `main` afterwards does not alter an existing
record.

## References

### Data

The coupling annotation and the surface-expression comparison are not ours; cite them in
their own right. Redistribution terms are in [`data/external/README.md`](data/external/README.md).

- **Surface expression and the library** — Tedman A, Goel M, Shah S, Howard MK, Chamness LM,
  Bonifasi A, Adams I, Gallagher JM, Penn WD, Nemec K, McDonald EF, Corman BN, Robinson JP,
  Post CB, Clark PL, Babu MM, Manglik A, Kuntz CP, Coyote-Maestas W, Schlebach JP. Efficient
  experimental characterization of the GPCRome via deep receptor scanning. *Nature
  Communications*. 2026;17(1):9796. doi:[10.1038/s41467-026-76564-7](https://doi.org/10.1038/s41467-026-76564-7)
- **G-protein coupling** — Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH,
  Davenport AP, Spedding M, Davies JA. The IUPHAR/BPS Guide to PHARMACOLOGY in 2026.
  *Nucleic Acids Research*. 2026;54(D1):D1446–D1456.
  doi:[10.1093/nar/gkaf1067](https://doi.org/10.1093/nar/gkaf1067)
  <br>`build_coupling.py` reads only the `GtoPdb` sheet, so this is the provenance of every
  coupling call.
- **The workbook that sheet was downloaded from** — GproteinDb asks that both its own
  reference and the main GPCRdb reference be cited:
  - Pándy-Szekeres G, Taracena Herrera LP, Caroli J, Kermani AA, Kulkarni Y, Keserű GM,
    Gloriam DE. GproteinDb in 2024: new G protein-GPCR couplings, AlphaFold2-multimer models
    and interface interactions. *Nucleic Acids Research*. 2024;52(D1):D466–D475.
    doi:[10.1093/nar/gkad1089](https://doi.org/10.1093/nar/gkad1089)
  - Taracena Herrera LP, Andreassen SN, Caroli J, Rodríguez-Espigares I, Kermani AA,
    Keserű GM, Kooistra AJ, Pándy-Szekeres G, Gloriam DE. GPCRdb in 2025: adding odorant
    receptors, data mapper, structure similarity search and models of physiological ligand
    complexes. *Nucleic Acids Research*. 2025;53(D1):D425–D435.
    doi:[10.1093/nar/gkae1065](https://doi.org/10.1093/nar/gkae1065)

### Software

Versions are the ones this screen was actually run with, read from the run's own logs and
from `requirements.txt`.

| Tool | Version | Citation |
|---|---|---|
| fastp | 1.3.2 | Chen S. fastp 1.0: an ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta*. 2025;4(5):e70078. doi:[10.1002/imt2.70078](https://doi.org/10.1002/imt2.70078) |
| cutadapt | 5.0 | Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*. 2011;17(1):10–12. doi:[10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200) |
| Biopython | 1.83 | Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25(11):1422–1423. doi:[10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) |
| pandas | 3.0.1 | The pandas development team. pandas-dev/pandas: Pandas. Zenodo. doi:[10.5281/zenodo.10697587](https://doi.org/10.5281/zenodo.10697587) · McKinney W. Data structures for statistical computing in Python. *Proc. 9th Python in Science Conf.* 2010:56–61. doi:[10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a) |
| NumPy | 2.4.3 | Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585(7825):357–362. doi:[10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| SciPy | 1.17.1 | Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*. 2020;17(3):261–272. doi:[10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) |
| matplotlib | 3.10.8 | Hunter JD. Matplotlib: a 2D graphics environment. *Computing in Science & Engineering*. 2007;9(3):90–95. doi:[10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55) |

fastp is cited as its 1.x paper because version 1.3.2 was used; the 2018 *Bioinformatics*
paper covers fastp 0.x. openpyxl and pillow have no citable publication.

### Identifiers

Used as join keys and labels only — no data was taken from these resources, so cite them
only if you need identifier provenance. Gene symbols follow **HGNC**
(doi:[10.1093/nar/gkaf1229](https://doi.org/10.1093/nar/gkaf1229)), transcript IDs follow
**Ensembl** (doi:[10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239)), and protein
accessions follow **UniProt** (doi:[10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)).

## Licence

Copyright (c) 2026 Matthew K. Howard.

Code in this repository is released under the MIT licence (`LICENSE`); the data generated
here is released under CC BY 4.0. Third-party reference data under `data/external/` is not
covered by either and keeps its own licence and attribution — see
[`data/external/README.md`](data/external/README.md).
