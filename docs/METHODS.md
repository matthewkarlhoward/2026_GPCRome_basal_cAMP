# Methods and caveats

What the panels measure, how the numbers were produced, and what not to read into them.
Distilled from the working analysis record; the full exploratory analysis is not part of
this release.

## Design

A barcoded GPCRome library of 946 constructs is integrated into a HEK293T CRE–eGFP
reporter line, grown without added agonist, and sorted into four bins (`a`–`d`) of
increasing reporter fluorescence. Each bin is sequenced; a construct's basal activity is
where its barcode sits across those bins.

Two sort days of four bins each were run:

| Day | Samples | Used |
|---|---|---|
| 20260924 | s5–s8 | replicate 1 |
| 20260926 | s9–s12 | replicate 2 |

The screen has **two replicate sorts**, each four bins. A construct is scored only when it
clears the read threshold in both, and `scripts/score_basal.py` asserts that the two days
it reads are the two in the release.



## Scoring

The unit is the **barcode**, never the gene: canonical and alternative isoforms of one
receptor are separate constructs and are never collapsed.

Within each sort day, each bin's counts are converted to a frequency within that bin's
library-matched reads, so sequencing depth cancels:

```
WMB_day = Σ w·f / Σ f        w = a:1, b:2, c:3, d:4
```

This assumes an equal number of cells sorted into each bin. A day counts for a barcode
only if that barcode has **≥50 reads** summed over its four bins.

`score` is the mean over the two retained days, rescaled as `(WMB − 1) / 3`, and
`score_sd` is the SD between days. **901 of 946 barcodes are scored on both days**
(727 canonical + 174 alternative isoforms) and 930 on at least one.

The 0–1 scale is a change of units only. It is a **relative activity index, not a cAMP
concentration**, and its endpoints are the bin bounds, not a biological zero or maximum.

## Data quality

- **Match rate** 93–96% of trimmed reads in every retained sample.
- **Coverage** 924 (20260924) and 907 (20260926) barcodes pass the read threshold.
- **Replicate agreement** r = 0.90, ρ = 0.88 (Supplemental Figure 1e).

## Coupling annotation (Figure 1b, 1d)

`coupling_display` is assigned in `build_coupling.py`, principally from the **IUPHAR/Guide
to Pharmacology** curation in the GproteinDb bulk table, which is external to this screen.
Receptors are labelled `Gs only`, `Gs + other`, `Gs secondary`, `Gq`, `Gi`,
`multi (no Gs)`, `Golf`, `gustducin` or `unassigned`. No layer is derived from our own
scores, but not every layer is external to this lab, and `coupling_source` records which
one resolved each barcode.

Where IUPHAR is silent the structural class from **Tedman et al. 2026** supplies `Golf` for
olfactory and `gustducin` for taste receptors, and a short literature set supplies
`Gs only` for receptors the IUPHAR transducer table omits entirely — GPR161 among them.
Where all of those are silent the last resort is an **in-house curated annotation**
(`data/inhouse/legacy_coupling.csv`, carried over from an earlier screen in this lab).

That in-house annotation decides the label for **42 of the 901 scored barcodes**: 40 through
the fallback layer itself, which stamps them `coupling_source = "legacy curated dict
(unverified)"`, and 2 more (TAAR6, TAAR9) that it calls olfactory although Tedman classes
them `A`. It is hand-written and unverified, and it is not independent of this lab, so
those 42 calls should be treated as weaker than the rest. It is nonetheless an input the
figures depend on: deleting it recolours Figure 1b, changes Figure 1d's groups and printed
`n`, and changes Supplemental Table 1's `coupling` column.

## Prior constitutive-activity evidence (Figure 1c)

`literature_constitutive_activity.csv` is a **hand-curated** table merged into
`canonical_scores_ranked.csv` as `ca_evidence`, `ca_mechanism`, `ca_caveat` and
`ca_source`. **Only the top 40 canonical receptors were reviewed.** A blank or
`not reviewed` means nobody looked, not that no evidence exists. Grades are `established`
(primary papers, structures or textbook-level), `reported`, `weak`, `none found` and
`not assessed`.

Of the top 40: 22 established, 6 reported, 3 weak/none, 9 not assessed; 23 carry a caveat.
Eleven receptors that are constitutively active and have no well-established endogenous
ligand land in the top 40: GPR6 (rank 1), GPR52 (7), GPR161 (9), GPR12 (10), GPR62 (12),
GPR101 (22), GPR26 (28), ADGRG5 (29), GPR3 (30), GPR21 (33) and GPR61 (39). This is a
looser set than the `ORPHAN` flag drawn in Figure 1c, and the two must not be conflated —
see the reconciliation below.

### Reconciling the Figure 1c keys with Supp. Table 1

Three different counts of the same 30 receptors were in circulation; these are the ones to
use. All are now printed by `scripts/make_supp_table.py`, which asserts them on every run.

| Key in Figure 1c | Count | Definition |
|---|---|---|
| filled dot, "prior report" | 23 | `ca_evidence` is `established` (18) or `reported` (5) |
| open dot, "weak or not assessed" | 7 | everything else: `not assessed` (6) plus GPR26, graded `weak` |
| square, "orphan" | 8 | the curated `ORPHAN` set, mirrored into Supp. Table 1 as `orphan_curated` |
| triangle, "agonist in situ" | 5 | the curated `AUTOCRINE` set, mirrored as `agonist_in_situ` |

Two traps produced the earlier mismatched numbers:

- **23 vs 24.** Counting "receptors with some prior claim" gives 24, because it folds in
  GPR26, whose evidence is graded `weak`. The figure's filled dot is `established` +
  `reported` only, so it is 23. The open dot was previously labelled "not assessed", which
  was wrong for GPR26; it now reads "weak or not assessed".
- **8 vs 9 orphans.** Two different routes both give 9, and neither matches the figure.
  Counting `ligand_family == "orphan"` adds GPR4, which is annotated an orphan but is a
  deorphanised proton sensor. Counting off the eleven-receptor list above adds ADGRG5,
  an adhesion receptor activated by its own tethered Stachel agonist — constitutively
  active, but not ligand-free. The figure excludes both, giving 8: GPR6, GPR52, GPR161,
  GPR12, GPR62, GPR101, GPR26, GPR3.

GPR4 therefore carries no flag in Figure 1c, although its pH caveat is real (its pH50 for
cAMP is close to 8.0, so it is more than half maximally active at the assay pH). That
caveat is carried in the text rather than in the panel.

The two flags drawn in the panel mark receptors that need a caveat rather than hiding them:

- **agonist in situ** — prostanoid and adenosine receptors, flagged because the medium is
  not exchanged for 24 h and an accumulating endogenous agonist would be read as basal
  activity. Treat this as a hypothesis we have not tested, not an established mechanism:
  no primary study shows HEK293 releasing prostaglandin at a concentration that would
  activate a transfected prostanoid receptor, and adenosine concentration in conditioned
  medium has not been measured here.
- **orphan** — receptors with no known endogenous ligand, by the curated `ORPHAN` set
  rather than the `ligand_family` column; 8 of the top 30.

## What to know before reusing these scores

**The readout is Gs/cAMP.** A Gi- or Gq-coupled receptor with genuinely high constitutive
activity does not score high here. "The top of the ranking is almost pure Gs" is partly
true by construction: the Gs enrichment measures assay performance as much as biology.

## Surface expression (Supplemental Figure 1f)

Compared against **Tedman et al. 2026** (*Nat Commun*, DOI 10.1038/s41467-026-76564-7),
which measured plasma-membrane expression for this exact library in the same HEK293T
system. Joined on Ensembl transcript ID; 859 of the 901 scored barcodes have a value.

The panel's point is that the two measurements are close to independent: a receptor sitting
at the plasma membrane is not thereby signalling. Across all constructs ρ = 0.05
(n.s.); within the Gs-inclusive set ρ = 0.27.

One caution:

- **The Tedman readout saturates at the top.** Cells were sorted into *quintiles* of
  staining intensity, so the top gate is open-ended and a per-receptor estimate
  reconstructed from bin occupancy cannot exceed the top bin's anchor. The reported CV
  collapses 11-fold above 3000, and 9.9% of receptors sit within 10% of the 6515 maximum.
  **Do not read fine structure above ~4000.**
