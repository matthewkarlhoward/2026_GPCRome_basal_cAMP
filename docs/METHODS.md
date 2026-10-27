# Methods and caveats

What the two figures measure, how the numbers were produced, and what not to read into
them. Distilled from the working analysis record; the full exploratory analysis is not
part of this release.

## Design

A barcoded GPCRome library of 946 constructs is integrated into a HEK293T CRE–eGFP
reporter line, grown without added agonist, and sorted into four bins (`a`–`d`) of
increasing reporter fluorescence. Each bin is sequenced; a construct's basal activity is
where its barcode sits across those bins.

Three sort days of four bins each were run, plus an unsorted baseline:

| Day | Samples | Used |
|---|---|---|
| 20260922 | s1–s4 | **no** |
| 20260924 | s5–s8 | yes (replicate 1) |
| 20260926 | s9–s12 | yes (replicate 2) |
| baseline | s13 | **no** — unsorted |

**Day 1 was excluded** because its four bins had near-identical barcode composition — the
sort did not separate cells. Its weighted mean bin correlated at r = −0.2 with the other
two days and no relabelling of bins rescued it. Day 1 and the baseline also carry an
off-library barcode (`TTTTGCGCTG`) absent from the 946-barcode reference and from every
earlier run, at 9–45% of reads against index-hopping levels in s5–s12. Only the two
retained replicates are included in this release; `docs/SEQUENCING_QC.md` records what was
seen in the excluded samples.

The bench sheet lists day 3 as 20240926; 2026 is assumed.

## Scoring

The unit is the **barcode**, never the gene: canonical and alternative isoforms of one
receptor are separate constructs and are never collapsed.

Within each sort day, each bin's counts are converted to a frequency within that bin's
library-matched reads, so sequencing depth cancels:

```
WMB_day = Σ w·f / Σ f        w = a:1, b:2, c:3, d:4
```

This assumes an equal number of cells sorted into each bin; set `CELLS` in
`score_basal.py` if bin cell counts are known. Weighting by raw counts instead agrees at
r = 0.999. A day counts for a barcode only if that barcode has **≥50 reads** summed over
its four bins.

`score` is the mean over the two retained days, rescaled as `(WMB − 1) / 3`; `score_sd` is
the SD between days; `wmb` keeps the raw 1–4 value. **901 of 946 barcodes are scored on
both days** (727 canonical + 174 alternative isoforms) and 930 on at least one.

The 0–1 scale is a change of units only. It is a **relative activity index, not a cAMP
concentration**, and its endpoints are the bin bounds, not a biological zero or maximum.

## Data quality

- **Match rate** 93–96% of trimmed reads in every retained sample.
- **Coverage** 924 (day 2) and 907 (day 3) barcodes pass the read threshold.
- **Library skew** the 90/10 abundance ratio is 205–285 — uneven but usable.
- **Replicate agreement** r = 0.90, ρ = 0.88 (Supplemental Figure 1e). Bland–Altman bias
  +0.006 with 95% limits of ±0.16 **on the 0–1 score** (+0.02 and ±0.48 on the raw 1–4 WMB).
- **Depth and precision** the median `score_sd` between days falls from 0.054 below 1,000
  reads (n = 261) to 0.017 above 10⁴ (n = 444), on the 0–1 score.
- **Bin depth is even** the eight bin libraries are sequenced to within 1.57× of each
  other, and **bin *d* is the deepest on neither day** (the deepest is bin *a* on
  20260924 and bin *c* on 20260926), so a high-scoring receptor cannot gain reads by
  sitting in a deeper bin.

## Coupling annotation (Figure 1b, 1d)

`coupling_display` is assigned in `build_coupling.py` from the **IUPHAR/Guide to
Pharmacology** curation in the GproteinDb bulk table, which is external to this screen and
not derived from our scores. Receptors are labelled `Gs only`, `Gs + other`,
`Gs secondary`, `Gq`, `Gi`, `multi (no Gs)`, `Golf`, `gustducin` or `unassigned`.

Two Gs definitions are carried separately and **an AUC must always state which one it
refers to**, because the separation weakens as the definition loosens. The comparison set
is the 338 **non-sensory** canonical receptors: olfactory (362) and taste (27) barcodes are
54% of the canonical library and cannot signal through a Gs/cAMP readout, so leaving them
in the "rest" group inflates every value below by 0.05–0.08.

| Definition | n | AUC vs non-sensory rest | p (one-sided Mann–Whitney) |
|---|---|---|---|
| Gs sole primary (`coupling_display == "Gs only"`) | 54 | 0.853 | 1e-16 |
| Gs among primaries (`gs_primary`) | 73 | 0.786 | 4e-14 |
| Gs primary or secondary (`gs_any`) | 91 | 0.750 | 9e-13 |

## Prior constitutive-activity evidence (Figure 1c)

`literature_constitutive_activity.csv` is a **hand-curated** table merged into
`canonical_scores_ranked.csv` as `ca_evidence`, `ca_mechanism`, `ca_caveat` and
`ca_source`. **Only the top 40 canonical receptors were reviewed.** A blank or
`not reviewed` means nobody looked, not that no evidence exists. Grades are `established`
(primary papers, structures or textbook-level), `reported`, `weak`, `none found` and
`not assessed`.

Of the top 40: 22 established, 7 reported, 2 weak/none, 9 not assessed; 23 carry a caveat.
All 11 well-documented constitutively active orphan receptors land in the top 40: GPR6 (rank 1), GPR52 (7), GPR161 (9), GPR12 (10), GPR62 (12), GPR101 (22), GPR26 (28), ADGRG5 (29), GPR3 (30), GPR21 (33) and GPR61 (39). This is a different, larger set than the `ORPHAN` flag drawn in Figure 1c, which marks orphans among the top 30 only.

The two flags drawn in the panel mark receptors that need a caveat rather than hiding them:

- **agonist in situ** — prostanoid and adenosine receptors, whose agonist these cells can
  supply (medium is not exchanged for 24 h and HEK293T make PGE2/PGD2). PTGIR is the
  internal control: same family, same transducer, comparable surface expression, rank 106,
  and its synthase is absent from HEK293T. The effect is receptor-specific, not a blanket
  artifact.
- **orphan** — receptors with no known endogenous ligand.

Secretin-family class B receptors (GIPR, VIPR2, PTH1R, PTH2R, SCTR, GLP1R, GHRHR) show
ligand-independent cAMP in proportion to expression level and should be read with that in
mind.

## Three things to know before reusing these scores

1. **The readout is Gs/cAMP.** A Gi- or Gq-coupled receptor with genuinely high
   constitutive activity does not score high here. "The top of the ranking is almost pure
   Gs" is partly true by construction, and the Gs AUC measures assay performance as much
   as biology.
2. **Score correlates with library abundance** (Spearman ρ ≈ 0.57) — the largest single
   confound in the screen. It is not a sequencing artifact (recomputing abundance from
   within-bin frequencies agrees at ρ = 0.9999) and not counting noise (it holds inside the
   high-abundance mode, ρ = +0.29). Composition explains roughly half: the library is
   bimodal, 74% of the low mode is olfactory, and 67% of all olfactory barcodes sit in it.
   Alternative isoforms are systematically better represented than canonical constructs
   (median 10³·¹⁷ vs 10²·⁴⁹ rpm), which matters for any canonical-vs-isoform comparison.
3. **Olfactory receptors drive a Simpson's paradox in anything involving surface
   expression.** They are ~40% of the scored library, ~42% of the constructs that join to Tedman, and carry
   their own internal negative basal~surface correlation (ρ = −0.44). Any filter that
   changes their share changes the apparent sign of the whole-library correlation.

## Surface expression (Supplemental Figure 1f)

Compared against **Tedman et al. 2026** (*Nat Commun*, DOI 10.1038/s41467-026-76564-7),
which measured plasma-membrane expression for this exact library in the same HEK293T
system. Joined on Ensembl transcript ID; 859 of the 901 scored barcodes have a value.

The panel's point is that the two measurements are close to independent: a receptor sitting
at the plasma membrane is not thereby signalling. Across all constructs ρ = 0.05
(n.s.); within the Gs-inclusive set ρ = 0.27.

Two cautions:

- **The Tedman readout saturates at the top.** Cells were sorted into *quintiles* of
  staining intensity, so the top gate is open-ended and a per-receptor estimate
  reconstructed from bin occupancy cannot exceed the top bin's anchor. The reported CV
  collapses 11-fold above 3000, and 9.9% of receptors sit within 10% of the 6515 maximum.
  **Do not read fine structure above ~4000.**
- **Do not substitute an in-house surface screen for an independent one.** On identical
  canonical non-olfactory barcodes (n = 265), an earlier same-lab surface screen correlates
  with basal score at ρ = +0.23 while Tedman gives ρ = +0.003 — yet the two surface screens
  agree with each other at ρ = +0.50. The residual is shared technical variance between two
  in-house screens, plausibly the same recombinant cell pool, not biology. That in-house
  screen is not part of this release; use Tedman for anything independent.
