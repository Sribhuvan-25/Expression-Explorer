# Methods and analysis decisions

Every analytical choice this tool makes, why it was made, and where it
lives in the code. Anything that changes a number a user reads belongs
here.

Each entry is dated and marked with its status:
**Decided** (settled, implemented) · **Provisional** (in use, open to
change) · **Open** (needs a decision).

When a decision changes, update the entry in place and add a line to
[Change log](#change-log) at the bottom — don't silently rewrite history.

---

## 1. Statistics

### 1.1 Group comparison — pairwise test
**Decided.** Mann–Whitney U, two-sided.
`backend/app/analysis/compare.py` → `pairwise_tests()`

Non-parametric, so it makes no normality assumption. Expression data is
routinely skewed and group sizes here are often small and unequal
(TARGET ETP n=19 vs non-ETP n=146), which is exactly where a t-test's
assumptions are least safe.

### 1.2 Group comparison — omnibus test
**Decided.** Kruskal–Wallis across all groups at once.
`backend/app/analysis/compare.py` → `kruskal_wallis()`

The non-parametric counterpart to one-way ANOVA, consistent with 1.1.
Reported alongside pairwise tests so a reader can see whether *any*
difference exists before reading individual pairs.

### 1.3 Multiple-testing correction
**Decided.** Benjamini–Hochberg FDR, applied across the pairwise tests
shown together. Reported as `q_value` next to each raw `p_value`.
`backend/app/analysis/compare.py:50-55`

Comparing 35 subtypes means 276 pairwise tests; uncorrected p-values
would be badly misleading. BH controls false-discovery rate rather than
family-wise error, which is the right trade-off for exploratory
screening. Verified against an independent implementation (max abs
difference 0.0 across all 276 tests) with q ≥ p and monotonicity holding.

Both raw p and q are shown — the correction is visible, not silently
applied.

### 1.4 Minimum group size for testing
**Decided.** A group needs **n ≥ 2** to enter a pairwise test or the
Kruskal–Wallis omnibus. Smaller groups are skipped, not dropped from
the chart.
`compare.py:45-46` (pairwise), `compare.py:62` (omnibus)

A single sample has no within-group variance, so a test against it is
meaningless. The sample still appears in the plot (as a point — see 3.3)
because it is real data worth seeing; it just isn't given a p-value.

---

## 2. Signature scoring

### 2.1 Default method
**Decided.** AUCell-equivalent rank-based scoring.
`backend/app/analysis/signature.py` → `auc_signature_score()`

Rank-based, so it's insensitive to absolute scale and to between-sample
normalisation differences. This matters more once multiple datasets are
in play.

### 2.2 aucMaxRank
**Decided.** `top_fraction = 0.25` — the top 25% of expressed genes per
sample. `signature.py:27`

Matches the reference method's stated setting for bulk RNA-seq.

### 2.3 Score normalisation
**Decided.** AUC is divided by `max_rank × n_signature_genes`, bounding
scores to **[0, 1]**. `signature.py:63-64`

Verified empirically: observed range 0.0–0.622 (DepMap ETP-TF5),
0.0–0.917 (TARGET). A score outside [0,1] would indicate a bug.

A score of exactly 0 is meaningful, not missing: it means none of the
signature's genes reached that sample's top 25%. 46 of 186 DepMap lines
score 0 on ETP-TF5.

### 2.4 Alternative method
**Decided.** `log2_mean` offered as a toggle — mean of log2(x+1) across
signature genes. Unbounded, unlike AUC. `signature.py` →
`log2_mean_signature_score()`

Kept because some published figures use a simple mean; the toggle makes
the choice explicit rather than hidden.

### 2.5 Unknown genes
**Decided.** Reject the whole request with a clear 404 naming the gene,
rather than silently dropping it. `backend/app/api/main.py` →
`_resolve_gene()`

Silently scoring a 5-gene signature as a 4-gene signature would change
the result with no indication. The library layer also filters missing
genes defensively, but the route stops it first.

---

## 3. Survival analysis

### 3.1 Group split
**Decided.** Binarise at the **median** signature score → LOW / HIGH.
`backend/app/analysis/survival.py` → `binarize_by_median()`

Matches the reference method's binarisation. Gives balanced arms
(233/233 on TARGET) and avoids a threshold chosen to maximise
separation, which would be circular.

Note `> median` → HIGH, so exact-median samples fall in LOW.

### 3.2 Event and duration
**Decided.**
- `event = 1` if `vital_status == "Dead"`, else 0 (right-censored)
- `duration = days_to_death` if the patient died, else
  `days_to_last_follow_up`

`survival.py:26-27`

### 3.3 Excluded patients
**Decided.** A patient with neither `days_to_death` nor
`days_to_last_follow_up` has no time axis and is dropped.
`survival.py:38`

On TARGET this is **3 of 469** (n = 466). This was previously invisible
and read as a mystery third dataset during expert review — the UI now
states the exclusion explicitly. See §6.1.

### 3.4 Tests reported
**Decided.** Log-rank test between arms, plus a Cox proportional-hazards
model. Cox reports coefficient, hazard ratio, **95% CI**, and p.
`survival.py` → `kaplan_meier_curves()`, `cox_model()`

The CI was previously computed by lifelines but dropped before display;
an HR without one gives no sense of precision.

### 3.5 Cox failures
**Decided.** If the Cox model cannot be fit (too few complete cases, an
unknown covariate, a degenerate/collinear design), return `cox: null`
and still show the KM curves and log-rank result.
`backend/app/api/main.py`, catching `ValueError` /
`ConvergenceError` / `StatError`

A failed model shouldn't destroy a valid curve.

### 3.6 Categorical covariates
**Open.** `cox_model()` does not one-hot encode covariates, so a raw
categorical column (e.g. `etp_status`) fails the fit and degrades to
`cox: null`. No UI currently exposes covariate selection, so this is
unreachable in the app today. Needs a decision before covariates are
surfaced.

---

## 4. Data handling

### 4.1 Grouping columns with partial coverage
**Decided.** Samples with no value for the grouping column are excluded
from that comparison, and the exclusion is reported.
`backend/app/analysis/compare.py` → `expression_by_group()`

This is substantial and easy to miss: grouping TARGET by `etp_status`
runs on **190 of 469** samples (279 excluded), because ETP status is
only classified for part of the cohort. `mrd_status` covers 265 of 469.

### 4.2 Missing metadata in output
**Decided.** Missing values serialise as JSON `null` (rendered "—"), not
`NaN`. `backend/app/analysis/ranking.py`

`NaN` is not valid JSON; Starlette rejects it with `allow_nan=False`,
which previously turned any ranking touching an unclassified sample into
an unhandled 500.

### 4.3 Which metadata appears in rankings
**Decided.** Every declared group column, whether it lives as a
top-level sample field or inside the nested `group_columns` dict.
Survival inputs (`days_to_death`, `days_to_last_follow_up`) are
deliberately **excluded**. `ranking.py` → `_NON_META_COLUMNS`

Survival fields are model inputs, not annotation someone scans while
picking candidates.

---

## 5. Visualisation

### 5.1 Log vs linear axis
**Decided.** Log axis auto-engages when `max / median > 10`; user can
override with a toggle, and the auto-choice is labelled "auto".
`frontend/src/components/BoxPlot.tsx:64`

Expression data is routinely right-skewed — on a linear axis a single
high outlier compresses every other group into a flat line. Defaulting
to linear meant users saw the unreadable view first.

Uses `log1p` / `expm1`, not bare log: true zeros are common and
`log(0)` is undefined, while `log1p(0) = 0`.

### 5.2 Axis headroom
**Decided.** 4% headroom added in **raw** space, before the log
transform. `BoxPlot.tsx:134-135`

Scaling the transformed max instead inflates exponentially once
untransformed for tick labels: a real max of 6.402 printed as 7.02, an
axis claiming an expression level that does not occur in the data.

### 5.3 Orientation
**Decided.** Vertical boxes up to 7 groups; **horizontal rows** beyond
that. `BoxPlot.tsx:12`

35 vertical bands need ~3000px to keep labels legible, forcing
horizontal scrolling that hides most of the chart. Horizontal rows put
names in a readable left column and let the chart grow downward.

### 5.4 Single-sample groups
**Decided.** A group with n = 1 renders as a **point**, not a box.

One sample has no quartile range; a 1px "box" reads as a rendering bug
and implies a distribution that isn't there.

### 5.5 What the statistics panel covers
**Decided.** Statistics are computed across **all** groups, independent
of the chart's group filter, and the panel says so when filtering is
active.

Filtering the chart is a display choice; silently recomputing p-values
to match would let a user change significance by hiding groups.

---

## 6. Reporting and transparency

### 6.1 Sample accounting
**Decided.** Any analysis running on fewer samples than the dataset
holds reports `n_dataset_total`, `n_excluded`, and a plain-language
reason, shown in the UI.

Added after expert review: "n = 466" against a listed 469 was read as an
undisclosed third dataset. Applies to both survival (§3.3) and grouped
comparison (§4.1).

### 6.2 Source attribution
**Decided.** Dataset and paper citations live **only** in the Sources
panel. Feature pages describe the operation, never a specific study.

The tool is general-purpose; page copy naming one study implies it isn't.

---

## 7. Planned — not yet implemented

### 7.1 GDS4299 / GSE28703 probe collapse
**Decided and implemented (2026-08-25).** `app/ingest/gds4299.py`.
For genes with multiple probes,
select the probe with the **highest mean expression** across all samples.

The selection rule must be **independent of the ETP / non-ETP grouping**.
Prior exploratory work (`Protein-Analysis/scripts/etp_analysis_optimal_probes.py`)
selected, per gene, the probe with the *lowest p-value between ETP and
non-ETP*. That is circular for inference — choosing the probe that best
separates the groups and then testing whether the groups separate
inflates significance. Acceptable for generating exploratory figures,
not for a tool that reports p-values as results.

Highest-mean-expression is outcome-independent and favours probes with
real signal over background. Alternatives considered: max variance
(outcome-independent but noise-prone), mean of all probes (unbiased but
diluted by poor probes).

**Empirically checked before adopting** (2026-08-25), across the MYCN /
PP2A / ETP-signature gene panel on GSE28703:

- The two rules **pick the same probe for 8 of 12** multi-probe genes,
  including MYCN.
- MYCN has **10 probes**; the highest-mean probe (`209757_PM_s_at`,
  mean 5.36) is unambiguous — the other nine sit at ~2.7, i.e.
  background, with scattered p-values (0.055–0.88) consistent with
  noise. Both rules select it.
- Where the rules differ (MYB, PPP2R5C, PPP2R1B, KIAA1524), the
  unbiased rule **still finds the comparison significant**
  (p = 0.008, 0.004, <0.001, 0.003).

So the unbiased rule costs essentially nothing in sensitivity here while
removing the circularity. Worth raising with the expert since it differs
from the earlier analysis, but it does not overturn its conclusions.

**Aliasing note:** CIP2A is not on GPL13158 under that symbol — it is
annotated as `KIAA1524`. Gene-symbol aliasing has to be handled at
lookup or genes will silently appear absent. See §7.4.

**Probe-to-gene mapping details** (carried over from the prior verified
analysis in `Protein-Analysis/scripts/etp_analysis_optimal_probes.py`,
which had already solved these):
- A probe may list several symbols (`A /// B`); take the first.
- Sentinel values `---` and `NA` are not gene symbols and must be
  dropped, or ~10k probes map to a junk "gene".
- Legacy symbols need mapping: `CIP2A → KIAA1524`, `PME1 → PPME1`.
- Of 54,715 probes on GPL13158, **44,238 carry a usable gene symbol**.

### 7.1a GDS4299 ETP labels — provenance
**Decided.** ETP / non-ETP status comes from a curated label table
(`Protein-Analysis/data/sample_labels.csv`), **not** from GEO metadata.

This needs stating because GEO's own `!Sample_characteristics_ch1` for
GSE28703 says only `"cell type: tumor cells"` for every sample — the
series carries **no ETP annotation at all**. Grouping is impossible
without an external label source, exactly as with TARGET (§ Liu 2017
enrichment, `app/ingest/liu2017_etp.py`).

**Verified before use:** the label table splits **12 ETP / 40 non-ETP**,
which matches the published composition of GSE28703 (Zhang et al.,
Nature 2012, 481:157-163; PMID 22237106) confirmed against an
independent source. All 52 labels join cleanly to the matrix with no
missing or extra samples.

### 7.2 Cross-dataset comparison
**Decided (2026-08-25), not yet built.** Datasets are shown
**side by side, each with its own axis and its own statistics.** No
pooling of values across datasets by default.

TARGET is RNA-seq TPM; GDS4299 is Affymetrix log2 intensity. Different
scales, dynamic ranges, and zero-behaviour — pooling raw values across
platforms is not valid, and platform/batch effect would typically
dominate the biological ETP signal.

Faceting answers "is the ETP effect reproducible across independent
cohorts?", which is a stronger claim than a single pooled p-value from
mismatched sources.

**Deferred, not rejected:** pooling after per-dataset z-score or rank
normalisation, and formal batch correction (ComBat). Either could be
added later as an explicit, clearly-labelled opt-in. Group sizes are the
motivation — pooling GDS4299 would take ETP from 19 → 31 (+63%) — so
this is worth revisiting with the domain expert.

### 7.3 AALL0434
**Open — blocked on data access.** Controlled-access via dbGaP
`phs002276`; needs a Data Use Certification and institutional sign-off.
Nothing to build until data is in hand — this is paperwork, not
engineering, and no public mirror substitutes for it.

Parked deliberately rather than chased: everything else proceeds on
public data. Picked up if/when access is confirmed.

### 7.4 Gene-symbol aliases
**Open.** `_resolve_gene()` (`backend/app/api/main.py`) matches an exact
symbol or feature_id only. `FeatureMetadata` already carries an
`aliases` field, but nothing populates or searches it.

This bites as soon as a second platform is added: **CIP2A** is annotated
as **KIAA1524** on GPL13158, so a user searching "CIP2A" gets "not found"
even though the gene is present. Older microarray annotations use
superseded symbols routinely.

Options: populate `aliases` per-dataset at ingest from the platform
annotation, or resolve through a shared symbol-alias table (e.g. NCBI
gene_info) at lookup. Needs a decision before the alias problem
multiplies across platforms.

---

## Change log

| Date | Change |
|---|---|
| 2026-08-25 | Document created. Existing decisions (§1–6) recorded retrospectively from the implementation and verified against source. |
| 2026-08-25 | §7.1 probe collapse → highest mean expression; §7.2 cross-dataset → side-by-side, no pooling. |
| 2026-08-25 | §6.1 sample accounting added, prompted by expert review of the n=466 survival count. |
| 2026-08-25 | §3.4 Cox hazard-ratio 95% CI now reported (was computed then dropped). |
| 2026-08-25 | GDS4299/GSE28703 ingested (`app/ingest/gds4299.py`). 21,596 genes × 52 samples (12 ETP / 40 non-ETP). Probe collapse and ETP labels verified — see §7.1, §7.1a. |
