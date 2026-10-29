# Third-party reference data

Redistributed here so the analysis is self-contained. Each dataset keeps its own licence
and should be cited in its own right, not as part of this work.

## `tedman2026/`

Tedman et al. 2026, *Nature Communications*, DOI
[10.1038/s41467-026-76564-7](https://doi.org/10.1038/s41467-026-76564-7) (PMC13575135).
Open access, CC BY 4.0. The one supplementary file used here is their Supplementary
Data 1:

Plasma-membrane expression for **this exact library** — 766 canonical GPCRs and 174 splice
variants, same HEK293T system, anti-HA immunostaining across five sorted quintiles. Joined
to our scores on Ensembl transcript ID. Several of this manuscript's authors are authors
of that paper.

One caution, documented in `../../docs/METHODS.md`:

- **The intensity estimate saturates at the top.** Cells were sorted into quintiles, so the
  top gate is open-ended and a per-receptor estimate reconstructed from bin occupancy
  cannot exceed the top bin's anchor. Reported CV collapses 11-fold above 3000. Do not read
  fine structure above ~4000.

| File | Sheet used for |
|---|---|
| `41467_2026_76564_MOESM3_ESM.xlsx` | per-construct plasma-membrane expression (Canonical, Isoforms) |

## `gproteindb/`

GproteinDb / GPCRdb bulk G-protein coupling table, <https://gproteindb.org/signprot/couplings>
(downloaded 2026-10-20; see `SOURCE.txt`). The `GtoPdb` sheet — IUPHAR/Guide to
Pharmacology curation — is what `build_coupling.py` uses for primary and secondary
transducers, and it resolves most of the library. It is **external to this screen** and is
not derived from our scores. It is not, however, the only layer: where it is silent
`build_coupling.py` falls back to Tedman's structural class and then to an in-house
annotation, which between them decide 42 barcodes. See `../../docs/METHODS.md`.
