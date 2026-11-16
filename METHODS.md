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

### 1.4a Small groups: exact Mann–Whitney, and the smallest p it can give
**Decided (2026-11-15).** `backend/app/analysis/compare.py` →
`mann_whitney_pvalues()`, `min_p_achievable()`; used by `/compare`,
`/compare-multi` and `/differential`.

**Exact test for small groups.** When the two groups admit at most
5,000 splits (`EXACT_SPLIT_LIMIT`; e.g. 2 vs 2 = 6, 5 vs 3 = 56, 7 vs 7 =
3,432) the p-value is the exact permutation distribution of U over every
split, average ranks for ties, two-sided as min(1, 2 × smaller tail) —
scipy's exact-test convention. On tie-free data it equals
`mannwhitneyu(method="exact")` to 1e-12; on tied data it equals
`scipy.stats.permutation_test` with full enumeration (both pinned in
`tests/test_mann_whitney_exact.py`). Larger groups use the normal
approximation with tie and continuity correction, now chosen
explicitly. Cost: 7 vs 7 across 65,535 genes takes 4.7 s.

**Why this replaced scipy's default (a defect found in review).** The
genome-wide scan called `mannwhitneyu(..., axis=1)`; in scipy 1.15 one
tied gene anywhere switches *every* gene to the normal approximation. At
2 vs 2 that served p = 0.245 for a cleanly separated gene (exact: 0.333)
and p = 0.2207 for every gene with one all-zero arm, so GSE324984's
"Rank genes" list was 500 genes with a silent arm, MYCL absent. At 3 vs 3
the approximation crossed 0.05 (p = 0.047) where the exact test cannot go
below 0.1. Single-gene Compare had the same fallback for tied genes.
**This changes numbers** for any comparison with ≤ 5,000 splits, including
lopsided ones such as 2 vs 30 that scipy previously approximated
(e.g. small DepMap subtypes), and only those; TARGET's and GDS4299's ETP
comparisons are far above the limit and unchanged.

**The floor is reported.** Every test carries `min_p_achievable` =
2 / C(n_a + n_b, n_a) on the exact path (null on the approximate path).
It genuinely bounds every served p now. At the November 2026 datasets'
sizes it *is* the result: 2 vs 2 cannot go below **0.333**, 5 vs 2 not
below **0.095**, 3 vs 2 not below **0.2**, 5 vs 3 not below **0.036**. The UI
states it wherever it exceeds 0.05, because "p = 0.33, not significant"
otherwise reads as "no effect" when the test could not have said
anything else. For replicated designed experiments, read the count-based
contrast in §11.2 instead.

Kruskal–Wallis (§1.2) still uses the chi-square approximation; at 4
groups of 2 its smallest attainable p is 0.083 (complete separation,
measured), so it cannot produce
a spurious p < 0.05 at these sizes, but it is approximate.

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

**Update 2026-11-15.** For datasets already stored on a log2 scale
(`log2_tpm`: DepMap; `log2_intensity`: GDS4299, GSE66586) the values are
averaged as they are; log2(x+1) is applied only to linear units (TPM).
Previously every matrix was transformed, so on three of five datasets
the score was a mean of log2(log2(x)+1) — GSE66586 samples at 9–11 log2
units scored about 3.3. `is_log2_unit` / `already_log2` in
`app/analysis/signature.py`, used by `/signature-score` and
`/rank-signature`. Verified: GSE66586 PP2A-inhibitor log2-mean scores
equal the raw series-matrix probe means exactly (highest-mean probe per
gene, §7.1).

### 2.4a Entering a gene set
**Decided (2026-11-15).** The gene-set tag input (Signature Score,
Sample Ranking, PCA, Survival) splits typed or pasted text with the gene
panel's parser (`frontend/src/lib/geneInput.ts`): spaces, commas,
newlines, "Heading:" prefixes and "(former name)" notes. Typing "CIP2A
SET ANP32A ARPP19 PABIR1" then Enter used to make one tag that no
dataset could resolve; text left in the box when the user clicked
Score was dropped and is now kept.

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
**Decided.** Three cutoff modes, chosen per query (median is still the
default): `backend/app/analysis/survival.py` → `binarize_by_cutoff()`
(added 2026-08-29; `binarize_by_median()` kept as the underlying
median-mode implementation, unchanged).

- **median** — binarise at the median signature score → LOW/HIGH,
  covering every sample. Matches the reference method's binarisation,
  gives balanced arms (233/233 on TARGET), avoids a threshold chosen to
  maximise separation (which would be circular). Note `> median` → HIGH,
  so exact-median samples fall in LOW. This is what "Survival" meant
  before cutoff modes existed, and remains the default.
- **quartile** — top 25% (HIGH) vs. bottom 25% (LOW) by score; the
  middle 50% are excluded from everything downstream (KM curves,
  log-rank, Cox), not assigned to either group. This is a materially
  different comparison from median-split (fewer samples, cleaner
  separation), modelled on GEPIA3's identical mode — see
  reference/gepia3-feature-review.md. The excluded middle band is
  additive to §3.3's existing "no follow-up time" exclusion, and both
  reasons are combined in one `exclusion_reason` string rather than
  only reporting one of the two causes. When scores tie exactly on a
  cutoff the arms are **not** a clean 25/25 split — see §3.1a, which
  governs what happens then and is user-selectable.
- **custom** — same idea as quartile but with independently chosen
  high/low percentiles (e.g. top 30% vs. bottom 10%), since a real use
  case is an asymmetric split, not just "quartile with a different
  number." Rejects (422) an invalid combination where
  `high_pct + low_pct > 100`, since the two groups would overlap.

Verified: median mode reproduces the exact prior reference numbers
(n=466/469, HIGH=233/LOW=233) through the new shared function, not just
independently by a separate code path — a regression here would have
been silent otherwise.

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

**Correction (2026-08-29):** a QA pass found the CI was computed by the
API and even typed in the frontend client, but the Survival page's
results table never actually rendered the column — "reported" was true
of the backend only. Fixed in `frontend/src/pages/SurvivalPage.tsx`; the
95% CI column is now confirmed rendering live (e.g. `9.209 – 74.688` for
the ETP-TF5 gene set on TARGET). Recorded here so this specific
computed-but-not-displayed gap doesn't quietly recur.

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

**Update 2026-11-15 — units already on a log scale.** `BoxPlot` now takes
the dataset's `expression_unit`. For `log2_*` and `delta_ct` units the
toggle is hidden (a further `log1p` would plot log of a log2 value) and
the axis is **fitted to the data with 8% padding instead of starting at
0** (`valueDomain`). Zero is meaningful on a TPM axis ("not expressed")
and stays there; on a log2 axis it is just a low reading. Found in
user-standpoint testing: GSE66586 PPP2R5B (7.5–10.0 log2 intensity) on a
0-based axis filled the top quarter of the plot and the MYCN-status
difference was barely visible. The axis title now names the unit
("EXPRESSION (LOG2 INTENSITY)"), so an exported PNG says what its numbers
are.

### 5.2 Axis headroom
**Decided.** 4% headroom added in **raw** space, before the log
transform. `BoxPlot.tsx:134-135`

Linear axes (2026-11-15) put ticks on round values — steps of 1, 2, 2.5
or 5 × 10^k (`niceTicks`) — with the domain widened to the outer ticks,
all printed at one precision; the previous evenly divided axis printed
"10.4, 8.35, 6.25, 4.18". The headroom rule above still governs the log
mode.

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

### 5.6 Histogram bin count adapts to cohort size
**Decided (2026-09-13), after QA.** The Signature Score distribution uses
Freedman–Diaconis bins with a `sqrt(n)` floor, clamped to [6, 36], rather
than a fixed 36. `frontend/src/pages/SignaturePage.tsx`

A fixed 36 bins left **39% of bins empty** on the 52-sample GDS4299
cohort — a gap-toothed comb of 1–2 sample bars, which cannot answer the
question the histogram exists to answer (whether the top scores are a
distinct population or the tail of one continuous distribution).

FD rather than Sturges because it keys off the IQR, so it stays sane on
the zero-inflated score distributions this app produces (§3.1a). The
`sqrt(n)` floor exists because FD alone is very conservative on these
tightly-clustered scores — it gave 9 bins for the 469-sample TARGET
cohort, which rendered perfectly well at 36, and coarsening a working
view is its own regression. Measured on real signature scores:

| dataset | n | bins | empty bins | was (fixed 36) |
|---|---|---|---|---|
| GDS4299 | 52 | 8 | **0%** | 39% |
| DepMap | 186 | 14 | 7% | — |
| TARGET | 469 | 22 | 5% | 6% |

### 5.7 Exported figures must identify their subject
**Decided (2026-09-13), after QA.** Every PNG export puts the analysed
gene (or gene set) in both the filename and the in-image title.

Expression Compare wrote every figure to `<dataset>-by-group.png` with
the dataset name as the only title, so exporting two genes silently
**overwrote the first file**, and the gene appeared nowhere inside the
image — a saved figure could not be identified after the fact, which
makes it useless as a record. Every other export site already named its
subject; this one and the PCA title (a bare "PCA") were the exceptions.

**2026-11-15 additions.** Box-plot PNGs print each pair's FDR q beside
its raw p (the subtitle says "with FDR correction"; the image printed
only raw p, e.g. 0.0357 for a pair whose q is 0.107). The box plot uses
the shared scatter/PCA palette: its own order put two near-identical
blues on groups 1 and 3, the PCA defect above, so GSE66586's
"non-malignant" and "MYCN-amplified" boxes looked alike. With ≤ 40
samples, sample dots are drawn larger and opaque (each one is the data
in a cell-line panel). The export stamp uses the local date (it was UTC,
so evening exports in the Americas carried tomorrow's date).

**Correlation and PCA scatter (2026-11-15).** The correlation scatter is
coloured by a grouping column (default: the dataset's first declared one;
"None" turns it off) and, whenever the points span more than one group,
says under the result and in the PNG that r pools every sample and a
between-group difference alone can produce a strong r. GSE66586 MYCN vs
PPP2R5B, r = −0.847 (recomputed from the raw series matrix: −0.84696,
p = 0.00199), is mostly the five MYCN-amplified lines sitting apart; one
colour hid that. Axes name the unit and use round ticks (§5.2). Both
endpoints attach each sample's declared group columns through one
helper, `_sample_groups` in `app/api/main.py`, which reads top-level
columns as well as the nested dict: PCA used to read only the dict, so
on TARGET (where `etp_status` is top-level) its default "Color by
etp_status" coloured nothing. Samples with no value in the colouring
column are grey, drawn underneath and counted in the legend ("no value
279"); they had fallen back to the first group's blue, so TARGET's 279
uncalled samples read as ETP. Legend and PNG key carry per-group n.

### 5.8 A disabled control must say why
**Decided (2026-09-21), found during regression testing.** The
Differential Genes tab's "Rank genes" button and both group selects go
disabled whenever the chosen grouping column has fewer than 2 distinct
values on the selected dataset — but for the single-value case
specifically, nothing on the page said so.

Root cause: `groupA`/`groupB` only get auto-populated by an effect gated
on `groupOptions.length >= 2`; with exactly one option (DepMap's own
`lineage`, filtered to Lymphoid-only upstream — see §7.2's dataset
filtering) that effect never runs, so both stay `""`. The pre-existing
"Group A and Group B must be different" hint requires `groupA` to be
truthy first, so it never fires either. The net result was a visibly
disabled control with zero explanatory text anywhere on the page —
worse than the "must be different" case, which at least says something.

Distinct from §7.10's `unavailable` flag (a column that is declared but
currently has zero values — e.g. a degraded upstream fetch): this is the
**structural** case, where the column genuinely, permanently has exactly
one value for this dataset. The API guarantees `values.length === 0`
implies `unavailable === true`, so the two conditions never overlap —
verified by a test asserting the single-value message and the
`unavailable` message are mutually exclusive for the same disabled-button
symptom.

Fixed with a third, explicit message: `"{column}" has only one value on
this dataset ({value}) — nothing to compare. Pick a different grouping
column.` Verified live: switching from DepMap/`lineage` (message shown,
button disabled) to TARGET-ALL-P2/`vital_status` (2 values, no message,
button enabled, ranking runs normally) toggles correctly with no false
positives on a working configuration.

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
  unbiased rule **still finds the comparison significant**: Mann–Whitney
  (§1.1, what Compare serves) **p = 0.010, 0.023, < 0.001, 0.0046**
  respectively (CIP2A is the current symbol for KIAA1524).
  *Revised 2026-11-15:* this bullet originally gave p = 0.008, 0.004,
  < 0.001, 0.003 without naming the test. Those are **Welch t-tests** from
  the exploratory script on the same probes (they recompute to 0.0076,
  0.0037, 0.0001, 0.0031), not the app's test, so they never matched the
  UI. The Mann–Whitney values above are recomputed from the GEO files and
  equal what the running app serves.

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

### 7.1b GSE66586 (GPL6244) — probe collapse, symbols, cell-line identity
**Decided and implemented (2026-11-15).** `app/ingest/gse66586.py`;
shared GEO helpers moved to `app/ingest/_geo.py`.

Neuroblastoma cell lines (Gu et al. 2016, PMID 26844271), Affymetrix
Human Gene 1.0 ST, RMA log2 as submitted. 10 samples: 5 MYCN-amplified,
3 MYCN non-amplified, 2 non-malignant (bone-marrow MSC, HCN-1A).

**Probe collapse.** The §7.1 rule, now shared code
(`_geo.collapse_probes`): first symbol per transcript cluster, `---`/`NA`
dropped, highest mean across all 10 samples wins. It never looks at MYCN
status. 22,195 of 33,297 clusters carry a symbol; they collapse to 20,002
genes (1,310 had more than one cluster). In the MYC/PP2A panel only
ARPP19 has two (7988963, mean 11.24, kept; 8110618, mean 9.97). Moving the
helpers out of `gds4299.py` was checked to leave GDS4299's matrix,
samples and features bit-identical (`assert_frame_equal`, exact).

**Symbols.** GPL6244's annotation (Aug 2016) carries two panel genes
under old names: `KIAA1524` (now CIP2A) and `FAM122A` (now PABIR1). Both
are renamed to the current HGNC symbol at ingest and the old name kept in
`aliases`. Neither current symbol has its own cluster on this platform,
so there is no collision; if one existed, the highest-mean rule would
pick between them.

**Labels** are parsed from the series matrix (`!Sample_source_name_ch1`
→ `mycn_status`, characteristics → `cell_type`, title → `cell_line`).
Ingest fails unless the parse reproduces the published 5 / 3 / 2 split
and every cell type agrees with its group.

**Not every line is independent.** GSM1625938 "BE(2)c" is ATCC CRL-2268
= SK-N-BE(2)-C (CVCL_0529). GSM1625941 is *titled* "SK-N-BE(2)c" but is
ATCC CRL-2271 = the parental SK-N-BE(2) (CVCL_0528), and its CEL file
name agrees. SH-SY5Y (CVCL_0019) is a subclone of SK-N-SH (CVCL_0531).
So the MYCN-amplified arm is 5 samples from 4 patients and the
non-amplified arm 3 samples from 2, and any group test here overstates
the number of independent biological samples. Pearson r on the collapsed
matrix: BE(2)-C vs SK-N-BE(2) 0.954, SH-SY5Y vs SK-N-SH 0.958 — against
0.960 for the unrelated SK-N-DZ vs IMR-32 and a 0.776–0.960 range over
all 45 pairs. Neither pair looks like two runs of one line, so all 10
samples are kept unmerged. `cell_line_identity` (Cellosaurus name,
resolved from each sample's ATCC catalogue number rather than its title)
is a declared group column so the relationship is visible. The dataset
notes disclose this, and the §1.4a floors for its group sizes.

Side observation, not acted on: HCN-1A correlates with hBM-MSC at 0.947
but at most 0.866 with any neuroblastoma line.

### 7.2 Cross-dataset comparison
**Decided and implemented (2026-08-27).** Datasets are shown
**side by side, each with its own axis and its own statistics.** No
pooling of values across datasets, ever, by default.
`GET /datasets/compare-multi` (`backend/app/api/main.py`) runs the
existing single-dataset `expression_by_group` independently per dataset
and returns one block per dataset — points, pairwise tests, Kruskal-Wallis,
assay type, expression unit, sample-exclusion accounting — with no shared
axis or merged statistic anywhere in the response.

**Unit labels must be distinguishable (2026-09-13, after QA.)** The
`ExpressionUnit` enum had no `log2_tpm`, so DepMap — RNA-seq, log2(TPM+1)
— was labelled `LOG2_INTENSITY`, the *microarray* unit. The UI then
displayed the identical string "log2 intensity" for DepMap and for
GDS4299 (genuinely microarray), which undercut the never-pool guarantee
printed right beside it: two units that are not comparable read as the
same unit. Added `ExpressionUnit.LOG2_TPM`; DepMap now reads
`rna-seq · log2 tpm` against GDS4299's `microarray · log2 intensity`.
Nothing branches on the value (it is display-only), so this changes no
computed number — but a reader can no longer mistake the two for one
scale.

TARGET is RNA-seq TPM; GDS4299 is Affymetrix log2 intensity. Different
scales, dynamic ranges, and zero-behaviour — pooling raw values across
platforms is not valid, and platform/batch effect would typically
dominate the biological ETP signal.

Faceting answers "is the ETP effect reproducible across independent
cohorts?", which is a stronger claim than a single pooled p-value from
mismatched sources.

**Compatibility rule (which datasets can appear together):** a dataset is
eligible for a given comparison if it declares the requested
`group_column` in its `DatasetDescriptor.group_columns` — nothing about
assay type or expression unit gates eligibility, since results are never
pooled across datasets anyway (that's what makes mixing RNA-seq and
microarray in one faceted view safe here but would NOT make it safe to
pool their raw values). Concretely: grouping by `etp_status` makes
`target_all_p2` and `gds4299` both eligible and `depmap` ineligible
(no `etp_status` column); grouping by `lineage` makes only `depmap`
eligible.

**Default selection:** every eligible dataset is included by default —
"use everything compatible" is the starting state, not an opt-in. The
user can narrow the selection with per-dataset checkboxes
(`frontend/src/pages/ComparePage.tsx`); an ineligible dataset is not shown
disabled, it's omitted from the checkbox list entirely and named in a
one-line note instead ("Not shown (no '<column>' column): ..."), so
absence reads as "doesn't apply" rather than "failed" or "forgot to
check it." Changing the grouping column resets the selection back to
"all eligible for the new column" rather than silently carrying forward a
selection that may no longer apply.

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
**Partially resolved (2026-08-29) — annotation done, resolution still
open.** `_resolve_gene()` (`backend/app/api/main.py`) still matches an
exact symbol or feature_id only; **this has not changed**, and gene
*lookup for analysis* (compare/rank/signature/survival) is unaffected by
what follows.

What's new: a read-only gene annotation panel
(`backend/app/services/gene_info.py`, `GET /genes/{symbol}`, rendered by
`frontend/src/components/GeneAnnotation.tsx` on the Compare page) queries
MyGene.info -- free, public, no auth (verified directly against the live
API). It resolves aliases for *display purposes*: searching "CIP2A" in
the annotation panel finds the gene even though GPL13158 only knows it as
KIAA1524, because MyGene.info's own alias index carries both names.
Verified directly: `mygene.info/v3/query?q=KIAA1524` returns the CIP2A
record.

**Update (2026-11-15).** Two narrow steps, neither of which changes
which gene an existing query matches. (1) When a query misses but the
dataset's *own* feature table records the name as an older symbol,
`_resolve_gene()`'s 404 now names the current symbol ("KIAA1524 … is an
older name for CIP2A") instead of a bare "not found", which read as a
coverage claim. (2) The design-contrast gene panel (§11.4) — and only it —
resolves such a name, and says so on that row. The open question below
is unchanged for every other endpoint. (3) Mixed-case symbols resolve:
the resolver upper-cased the query and then required an exact symbol
match, so `C1orf112`-style symbols could never be found on an
Ensembl-indexed dataset (TARGET, GSE324984). It now tries exact, then
upper-case, then a case-insensitive match that resolves only when
unambiguous (two symbols differing only by case → 422). `/compare-multi`
skip reasons now carry the legacy-symbol hint too.

**Deliberately not wired into `_resolve_gene()` yet.** Making alias
resolution affect which gene a *query* matches is a different, bigger
decision than showing alias information underneath an already-resolved
gene -- it would mean a dataset's own exact-match feature table is no
longer the sole authority on what "found" means for that dataset, which
needs its own explicit sign-off (e.g. does a GDS4299 query for "CIP2A"
silently start matching the KIAA1524 probe row, and if so, is that
disclosed the same way the probe-collapse rule in §7.1 is). Revisit this
half of the decision separately; don't fold it in silently as a side
effect of the annotation panel shipping.

### 7.5 DepMap release staleness
**Decided and implemented (2026-08-29).** The DepMap loader
(`backend/app/ingest/depmap.py`) can only ever find **DepMap 24Q4 Public**
(Dec 2024) going forward, and there is currently no fix available — this
is a bounded, monitored limitation, not a bug to chase further right now.

**Verified directly, not assumed:**
- `depmap.org/portal/api/download/files` — the endpoint DepMap staff
  themselves recommend for bulk/programmatic access (per the DepMap
  community forum) — returns a Cloudflare Turnstile bot-verification page
  (confirmed via a direct request: `content-type: text/html`, a
  "DepMap — Verification" challenge page, not CSV data). This is not new;
  it's the same block the original loader docstring already described.
- Figshare's public search API was queried directly for every article
  titled `DepMap <YY>Q<N> Public`: the newest match is **24Q4**
  (published 2024-12-10). Nothing for 25Q1 onward exists under that
  title pattern — DepMap's own forum confirms they stopped publishing
  full quarterly bundles to Figshare after 24Q4, moving newer releases
  portal-only (25Q2, 25Q3, ...), which is exactly the interface that's
  bot-gated above.
- So there is currently **no scriptable path** to anything newer than
  24Q4. Fetching a fresh Figshare snapshot is not the fix — the loader
  was already fetching the newest thing Figshare has; the actual problem
  was that this fact was invisible.

**What was actually wrong, and the fix:** not that the app was silently
lying about its data (the loaded `DatasetSource.accession` already
correctly says "DepMap 24Q4 Public", verified against the real title),
but that nothing measured or surfaced *how far behind* that release was,
and nothing bounded how much staler it could silently get if this goes
unreviewed for years. `_release_age_warning()` parses the release title,
computes how many quarters behind the current calendar quarter it is,
and past a 2-quarter threshold appends a plain-language note to
`DatasetSource.notes` (shown wherever the dataset's provenance is
surfaced) and logs a warning at load time. Below the threshold, nothing
is appended — a release that's merely a quarter or two behind Figshare's
own publishing lag isn't treated as suspicious.

**Deliberately not attempted:** scripting around the Cloudflare
Turnstile challenge on DepMap's own portal API. That would mean solving
a bot-detection challenge designed specifically to stop this kind of
automated access — fragile even if it worked once, likely to break
silently on DepMap's next Cloudflare config change, and arguably against
the spirit of why that gate exists. If DepMap re-enables Figshare
publishing, or a legitimate authenticated API appears, this loader picks
it up automatically (same Figshare-search mechanism, no code change
needed) — revisit then, not by working around the gate now.

### 7.6 Gene-vs-gene correlation and co-expression network
**Decided and implemented (2026-08-29).**
`backend/app/analysis/correlation.py` — `GET
/datasets/{id}/correlation?gene_a=X&gene_b=Y&method=...` and `GET
/datasets/{id}/co-expression?gene=X&top_n=...&direction=...`. Frontend at
`frontend/src/pages/CorrelationPage.tsx` (new "Correlation" pane).
Identified as a cheap win against GEPIA3's Correlation Analysis /
Expression Network features (reference/gepia3-feature-review.md) —
"cheap" specifically because both are computable from a dataset's
existing expression matrix with zero new ingestion.

**Correlation:** Pearson / Spearman / Kendall between two genes across a
dataset's samples, scatter plot + coefficient + p-value. Samples missing
a value for either gene are dropped before computing (not zero-filled or
otherwise imputed); fewer than 3 paired samples raises rather than
returning a degenerate/undefined correlation.

**Co-expression network:** the top-N genes most (or least) correlated
with a query gene, via a single vectorised `DataFrame.corrwith` sweep
across the whole matrix rather than a per-gene Python loop — verified
this matters at this matrix's real size (tens of thousands of rows for
TARGET), not just in principle. Restricted to Pearson/Spearman for this
sweep specifically (excludes Kendall): pandas has no vectorised
Kendall-tau, so an O(genes) sweep would mean one O(n log n) kendall call
per gene, prohibitively slow, unlike the single-call pearson/spearman
path. `direction=negative` returns the most *negative* coefficients
specifically (true anti-correlation), not the largest `|r|` regardless of
sign — a user asking for anti-correlation wants the negative end, not
whichever sign happens to dominate.

**Never pools across datasets**, matching §7.2's rule: each dataset's
correlation is computed and shown independently; a caller wanting the
same gene pair across cohorts runs this once per dataset.

**Verified against real biology, not just unit tests:** on GDS4299
(ETP-ALL), MYCN's top anti-correlates are RAG1/CD1B/CD1E — T-cell
differentiation markers consistent with MYCN marking the immature ETP
state — and its top co-expression network hit on TARGET-ALL-P2 is
MYCNOS (the MYCN antisense RNA transcript), which is exactly the
correlation this method should surface if it's working correctly. Symbol
mapping confirmed correct on both an Ensembl-ID-keyed matrix (TARGET) and
a symbol-keyed matrix (GDS4299) — the network result always reports
`symbol`, regardless of which id space the underlying dataset's matrix
index uses.

**QA pass (2026-08-30, 2 parallel agents) confirmed correct, no changes
needed:**
- Every reference number above independently reproduced via a *second*,
  from-scratch recomputation (loading the raw matrix directly, not
  calling this app's own code at all) — bit-identical to full float
  precision, not just "close."
- `n` always equals `len(points)`; the query gene is structurally
  excluded from its own network (via `matrix.drop(index=gene)`, not an
  incidental filter); `direction=negative` genuinely sorts most-negative
  first (checked strictly monotonic, not just "some negative genes
  present"); the kendall/pearson/spearman asymmetry between
  `/correlation` and `/co-expression` is a deliberate `Literal` type
  difference, not an accidental omission.
- Noted, not a bug: `/correlation` (scipy) and `/co-expression` (pandas
  `.corrwith`) can differ in the last float digit (~1e-16 relative) for
  the same gene pair — two correct numeric implementations of the same
  statistic, not a logic error.
- Noted, a real design choice worth being explicit about: `top_n` is
  hard-capped at 200 (422 error above that), not silently clamped to
  "everything available." Deliberate, not accidental — 200 is already a
  generous network size — but flagging here so a future reader doesn't
  mistake the cap for an oversight.

**Two UI bugs found and fixed (2026-08-30):**
1. Switching from the Correlation tab to Co-expression network and back
   reset the Correlation tab's results to empty, even though the query
   inputs looked unchanged (they'd just reset to the same defaults on
   remount). Root cause: `CorrelationPage.tsx` conditionally rendered
   `{tab === "correlation" && <CorrelationTab/>}`, which unmounts the
   inactive tab and destroys its local state entirely. Fixed by keeping
   both tabs permanently mounted and toggling visibility with CSS
   (`hidden`/`contents`) instead of conditional rendering.
2. At split-pane widths around 700-900px, the Correlation tab's
   Pearson/Spearman/Kendall method toggle clipped "Kendall" to "Ken"
   (confirmed via `getBoundingClientRect`: the button extended ~35px past
   its `overflow-hidden` container). The Co-expression tab's own method
   toggle (Pearson/Spearman only, one fewer/shorter option) didn't hit
   this because it never needed as much width. Fixed by giving the method
   group a wider minimum width and letting it wrap onto a second row if
   still squeezed, instead of silently clipping.

### 7.7 PCA / dimensionality reduction
**Decided and implemented (2026-08-30).** `backend/app/analysis/
dimensionality.py`, `POST /datasets/{id}/pca`, new "PCA" tab on the
Correlation pane (`frontend/src/pages/CorrelationPage.tsx`). Identified
against GEPIA3's Dimensionality Reduction feature
(reference/gepia3-feature-review.md) — an unbiased "do these samples
actually separate structurally" check, e.g. does the ETP-TF5 signature
cluster ETP vs non-ETP on more than the single-gene view Compare shows.

**No default gene set — user-supplied only, matching GEPIA3 exactly.**
Deliberate: PCA on all ~20-60k genes in a dataset is noisy and slow, but
an implicit "top N most-variable genes" default was considered and
rejected in favour of always requiring the caller to name the genes —
same reasoning as not silently inferring intent elsewhere in this app.

**Implemented via SVD directly (numpy), not scikit-learn** — PCA-via-SVD
is a handful of lines (center, standardize, `np.linalg.svd`), not a
fragile reimplementation of something genuinely complex, and adding a
dependency as heavy as scikit-learn for one function isn't proportionate.
Cross-checked at test time against sklearn's own `PCA` (available in the
dev environment, not a runtime dependency) on identical standardized
input — variance-explained ratios matched to 1e-8, confirming the
hand-rolled implementation is mathematically correct, not just internally
self-consistent.

**Standardises before the SVD** (zero mean, unit variance per gene —
i.e. PCA on the correlation matrix, not covariance): without it, a
handful of high-expression genes would dominate every component
regardless of whether they actually separate samples, defeating the
point of using a curated gene set. A gene with zero variance across every
sample is dropped rather than propagated as NaN (division by zero under
standardisation), and reported back in `genes_zero_variance` so the
caller can see what was silently excluded.

**Verified against real biology:** on GDS4299, the ETP-TF5 signature
(MEF2C/LYL1/HHEX/LMO2/MYCN) gives PC1 explaining 73.2% of variance, and
ETP samples cluster visibly apart from non-ETP on PC1 (ETP range
-4.08 to -0.94; non-ETP range -2.24 to 2.53) — confirmed both via the
raw API response and by looking at the actual rendered, color-coded
scatter plot in the browser, not just the numbers.

### 7.8 Genome-wide top differential genes
**Decided and implemented (2026-08-30).** `backend/app/analysis/
differential.py`, `GET /datasets/{id}/differential`, new "Top
differential genes" tab on the Expression Compare pane. The inverse of
Compare's existing single-gene view: instead of "how does this gene
differ across groups", this answers "which genes differ most between
these two groups" — modelled on GEPIA3's Differential Genes feature.

**Two-group only, deliberately.** An omnibus (Kruskal-Wallis) scan across
more than two groups has no single natural "ranking by effect" (which
group is "up" vs "down"?) the way a two-group test does. Compare's
existing `/compare` endpoint already covers the >2-group case for one
gene; this endpoint's job is specifically "rank every gene for A vs B."

**Vectorised across every gene in one scipy call
(`mannwhitneyu(..., axis=1)`), not a per-gene Python loop — measured, not
assumed.** A naive loop over TARGET's real 60,660 genes took ~32s; the
vectorised call takes ~2s. Verified the two approaches produce
bit-identical p-values before committing to the vectorised path, so this
is a performance change with no numeric difference, not a different
(faster but less correct) statistic.

**Tie-break fixed (2026-11-15).** The code sorted by p only, although
this entry and the docstring said ties are broken by |median difference|.
At 2 vs 2 thousands of genes share the floor p, so the order among them
was arbitrary (MYCL, +230-fold in GSE324984, fell outside the top 500).
Now sorted by p, then |median difference| descending; tied p keep tied q.

**BH-FDR computed over every gene tested, not just the returned top N**
— an FDR computed only across the already-most-significant genes shown
would be a biased correction (same reasoning as the pairwise-test FDR in
`compare.py`, applied genome-wide here).

**A new small endpoint, `GET /datasets/{id}/group-values`, was added
alongside this** — not because the differential-genes feature strictly
needed it, but because the frontend's "pick Group A / Group B" dropdown
needs to know what values a grouping column actually has (e.g. "ETP" /
"non-ETP" / "near-ETP" for `etp_status`), and that information didn't
exist anywhere in the existing API surface. An earlier draft probed this
by calling `/compare` with a hardcoded gene name ("MYCN") and reading off
the distinct groups in the response — fragile, since it assumed that
gene existed in every dataset. `/group-values` reads distinct values
directly off `ds.samples`, independent of any gene, and is now the
correct building block for any future "pick a group value" UI.

**Observation, not a bug:** on TARGET-ALL-P2 (n_a=19, n_b=146), the very
top-ranked genes by p-value tend to be near-zero-in-one-group genes
(e.g. a gene expressed only in the smaller group) rather than the
well-known ETP markers (HHEX, LYL1) — this is a real, correct
consequence of Mann-Whitney on an imbalanced small group with many
zero/near-zero genes, not a statistics bug. HHEX and LYL1 (both in the
ETP-TF5 signature used elsewhere in this app) do rank highly (65th and
458th of 60,660) once checked directly, confirming the method is finding
real biology — just competing against genuinely low p-values from a
different, also-legitimate source at the very top of the list.

### 7.9 Dataset warm-up must not gate the deploy healthcheck
**Decided (2026-09-20), after a failed Railway deploy.**
`app/services/warmup.py` loads every dataset on a background thread
started by the app's lifespan hook. `/health` answers as soon as uvicorn
binds; `/readiness` reports per-dataset state.

The container previously ran `scripts/prewarm_cache.py` to completion
*before* uvicorn bound, which made deploy success a race against a
download. That race was lost twice:

| Timeout | Outcome |
|---|---|
| 30s | Every deploy failed. Build and Deploy succeeded; Healthcheck failed at exactly 30s. |
| 900s | Bought time instead of removing the race. Then §8's aux layers added a 739-file GDC MAF pull and a 198MB CRISPR matrix to the same blocking path. |
| 120s | The healthcheck no longer waits on data at all. |

Measured cold start after the aux layers, on a home connection:

| Stage | Seconds |
|---|---|
| target_all_p2 RNA-seq (530-file GDC pull) | ~300 |
| target_all_p2 MAF mutations (739 files @ 0.47s, measured) | ~348 |
| DepMap expression (124MB) | ~60 |
| DepMap CRISPR (198MB, added §8.3) | ~95 |
| DepMap PRISM (40MB, added §8.4) | ~25 |
| GDS4299 | ~20 |
| **Total** | **~848s of a 900s budget (94%)** |

~52s of headroom, from a home connection. GDC is typically slower to
datacenter IPs and the MAF loader retries up to 4x per flaky file, so
the real margin was plausibly negative — which matches the reported
failure.

**Raising the timeout again was rejected.** It is a treadmill: every
dataset added restarts the same race, and the failure mode is a dead
deploy rather than a slow one. Binding first makes dataset size cost
warm-up time instead of deploy success.

Two consequences that had to be handled rather than hidden:

- **`/datasets` must never call a loader for an unwarmed dataset.** It is
  the frontend's first call; blocking there would hang the entire UI
  behind a ~14-minute pull — the same failure wearing different clothes.
  Unwarmed datasets are returned with `warming: true` and no provenance.
- **The UI says "preparing — downloading source data…"** rather than
  rendering a name above a blank line, which reads as broken rather than
  pending, and polls every 5s so it flips to ready on its own.

Verified end to end against an empty cache: `/health` returned 200 in
**2s** (was ~848s), `/datasets` answered in **16ms** while a 300MB+
download ran, DepMap went ready at ~90s and immediately served its
reference value (MYB r=-0.652361, n=91) while the other two were still
downloading.

**Disk note:** the warm cache is now ~495MB (DepMap CRISPR alone is a
198MB parquet). The Railway volume must exceed that or the deploy fails
on disk instead of time.

### 7.10 Single-file clinical-supplement fetches need the same defense as the MAF/STAR-Counts fetches
**Decided (2026-09-21), after a production incident.** After §7.9's fix
landed and the deploy succeeded, `/readiness` showed `target_all_p2`
itself failing: `"Excel file format cannot be determined, you must
specify an engine manually."` after 502s. `depmap` and `gds4299` warmed
fine (75.7s and 3.0s) — this was isolated to one dataset's clinical
enrichment, and it took the whole 469-sample dataset down with it.

Root cause: `app/ingest/target_mrd.py`'s `load_mrd_status()` did a
single-shot `httpx.get(...)` for TARGET's clinical supplement `.xlsx`
with no retry and no check that the response body was actually an Excel
file before handing it to `pd.read_excel`. A 200 response whose body
isn't the expected file — an HTML error/redirect page, a truncated
transfer — passes `raise_for_status()` (it's a 200) and isn't an
`httpx.TransportError` (the connection didn't drop), so nothing in the
loader caught it before pandas' own format sniff failed with the
"cannot be determined" message. Fetched directly outside the app, the
same URL returned a valid file — the failure is intermittent, upstream,
and not reproducible on demand, exactly like the GDC connection drops
`gdc_target.py`'s `_download_star_counts()` already retries with the
comment *"observed repeatedly from Railway, not reproducible locally"*.
`liu2017_etp.py`'s `load_etp_status()` has the identical single-shot
shape (fetches a different third-party host, Springer Nature) and while
investigating this incident its endpoint was independently observed
timing out for 277s before eventually succeeding — the same failure
class, not yet triggered in production only because it hadn't happened
to fail on a warmup run yet.

**Two independent layers of defense, not one:**

1. **`app/ingest/_http.py`** — `fetch_with_retry()` generalises the
   existing MAF-download retry pattern to the single-file case, and adds
   the piece that pattern didn't need: `validate_excel_bytes()` checks
   the response for an XLSX/XLS magic number before parsing, so a bad
   body is retried (like a transport error) instead of reaching pandas
   as garbage. A real 4xx/5xx still raises immediately — that's not
   transient, retrying only delays a failure that won't change. Both
   `target_mrd.py` and `liu2017_etp.py` now go through it.
2. **`gdc_target.load()` no longer lets either fetch be fatal.**
   `etp_status`/`mrd_status` are supplementary grouping columns on an
   existing cohort, not the cohort itself — the same "additive, never
   fatal" contract §8.1's aux layers already follow. Each is wrapped in
   its own try/except; a failure after retries is logged and the dataset
   loads with that one column simply absent, rather than propagating up
   and killing the 60,660×469 expression matrix and survival support
   along with it.

Verified by mocking both sources to fail unconditionally: the dataset
still loads with the full matrix and correct `n_samples`, `vital_status`
(sourced from GDC's own clinical file, not either flaky third-party
endpoint) is unaffected, and `etp_status`/`mrd_status` are absent from
every sample's `group_columns` rather than present-but-wrong. A second
test confirms partial degradation (only `mrd_status` failing leaves
`etp_status` intact) — the two sources are independent and one failing
must not take the other down too.

**A related, narrower gap closed at the same time:** `GET
/datasets/{id}/group-values` could return `{"values": []}` for a column
that is genuinely declared on the dataset but has zero non-null values
on the current load (exactly what a degraded `etp_status` fetch would
produce) — indistinguishable from the "unknown column" case §8.6-era QA
already fixed once for a different reason. Added an explicit
`"unavailable": true` flag so the UI states why a picker is empty
("no samples currently have a value... may be temporarily unavailable")
rather than rendering it with no options and no explanation.

---

## 8. Auxiliary measurement layers

### 8.1 Aux layers are not separate datasets
**Decided and implemented (2026-09-13).** CRISPR gene-effect, drug
sensitivity, and somatic mutation status attach to an **existing**
dataset as `AuxLayer` entries (`backend/app/models/contract.py`), not as
new registry entries.

They are the same samples, keyed by the same `sample_id`s, measured a
second way — DepMap's CRISPR screen runs on the same cell lines as its
expression matrix; TARGET's exomes are the same patients as its RNA-seq.
Registering "DepMap CRISPR" beside "DepMap cell line RNA-seq" in the
sidebar would imply two cohorts where there is one, which is exactly the
kind of quiet misrepresentation §7.2 and §6.1 exist to prevent.

`Dataset.aux` defaults to `{}`, so every pre-existing loader is
unaffected. `require_aux()` fails with a message naming what the dataset
*does* carry rather than a bare `KeyError`, and `aux_sample_overlap()`
returns `(covered, total)` because **no aux layer covers a whole
dataset** (see 8.2–8.4) and reporting that split is mandatory.

**Aux loading is opt-in (`with_aux`, default True in app use, passed
False by tests).** Each layer reaches an additional network service; the
primary expression loader must not transitively depend on them. This was
found the hard way — the first implementation made every existing
`gdc_target` unit test start downloading 739 MAF files.

### 8.2 TARGET-ALL-P2 somatic mutations
**Implemented.** `backend/app/ingest/target_mutations.py`. 739
open-access Masked Somatic Mutation MAFs from GDC → a boolean
(gene × sample) matrix: 6,250 genes × 717 samples, 10,978 mutated
gene-sample pairs, ~6 min cold ingest.

**Coverage: 279 of 469 RNA-seq samples (59%).** Predicted from research
before building, then reproduced exactly by the finished loader — a
useful independent check that the barcode→submitter_id join is right.
Every analysis over this layer reports the exclusion (§6.1).

**Non-silent filter is an explicit allowlist**, not "not Silent":
missense, nonsense, nonstop, frameshift ins/del, in-frame ins/del,
splice-site, translation-start. Verified against a real file where 7 of
21 calls in one aliquot were `Silent` and 2 were `3'UTR` — counting
those as "this gene is mutated" would inflate every mutation rate with
variants that don't change the protein.

**Validated against known T-ALL biology**, which is the strongest
evidence the parsing is correct: the top recurrently-mutated genes come
out as NOTCH1 57%, FBXW7 21%, PHF6 14%, NRAS 10%, PTEN 10%, JAK3 6%,
RPL10 5% — the canonical T-ALL driver profile at literature-consistent
frequencies. Splitting DTX1 (a direct NOTCH1 target) by NOTCH1 mutation
status gives **p = 6.1e-3** (medians 18.86 mutated vs 8.71 wild-type,
159 vs 120), while MYCN (not NOTCH1-regulated) gives **p = 0.75** — the
values the running app serves. *Revised 2026-11-15:* first measured as
p = 4.9e-3 (medians 18.86 vs 9.69) and p = 0.87 on 2026-09-13; the §8.6a
deterministic file-selection fix (2026-09-22) changed 29 cases'
expression, and the old values reproduce exactly from the pre-fix
backup. Same conclusion either way.

Gene-level only for now; amino-acid-level resolution is deferred (the
loader keeps the columns that would allow it).

### 8.3 DepMap CRISPR gene effect
**Implemented.** `backend/app/ingest/depmap_crispr.py`. 17,916 genes ×
1,178 cell lines, from `CRISPRGeneEffect.csv` in the *same* Figshare
release as the expression matrix — deliberately the same release, since
pulling CRISPR from a different one would silently mismatch the two
layers' cell-line rosters. Inherits §7.5's 24Q4 staleness limitation for
the same reason.

**Coverage: 91 of 186 lymphoid lines (49%).**

Chronos scale: 0 = no effect, ≈-1 = median common-essential, so **more
negative means stronger dependency** — the opposite direction of
intuition for an expression value, which is why `AuxMatrix` carries
`value_label`/`value_description` rather than leaving a UI to guess.

**Validated:** RPL13A (essential ribosomal protein) has median gene
effect −1.995; A1BG (non-essential) −0.080. MYB expression correlates
with MYB dependency at r = −0.652, p = 2.5e-12 — lines expressing more
MYB depend on it more, a known oncogene-addiction relationship in
leukemia.

### 8.4 DepMap PRISM drug sensitivity
**Implemented.** `backend/app/ingest/depmap_prism.py`. 6,790 compounds ×
919 cell lines. Two things make this more involved than the CRISPR
layer, both verified against the real files:

1. **Separate Figshare release** — ships as `Repurposing Public YYQN`,
   not inside the `DepMap YYQN Public` bundle, so it needs its own
   title regex and release resolution. File names are matched by
   *suffix* so a future release's renamed files still resolve.
2. **Transposed relative to every other matrix here** — the raw matrix
   is (compound × cell-line) with compounds as rows.

The matrix index is a Broad compound id (`BRD:BRD-K05804044-001-18-5`),
meaningless to a reader, so the compound list is joined in as the
layer's `features` frame carrying drug name and mechanism of action
(that id resolves to AZ-628). **Coverage: 79 of 186 lymphoid lines
(42%).** The compound→name join was verified exhaustively in QA:
6,789 of 6,789 correct on both name and MOA, with one id absent from the
upstream list and correctly falling back to the raw BRD id.

**Sparsity caveat — the 6,790 figure oversells this layer.** Against the
79 covered lines, **only ~1,518 compounds (22%) were assayed on at least
3 of them**; the other 78% return "not enough samples" on any lookup.
Exactly one compound (AZ-628) covers all 79 — and that was the compound
development testing happened to sample, which is precisely how this went
unnoticed until QA picked random compounds and hit failures on 2 of 4.
`AuxMatrix.n_usable_features()` now measures this per layer and
`/aux-layers` reports `n_usable_features` alongside `n_features`, so the
UI can state what is actually analysable rather than a true-but-useless
total. CRISPR has no such problem (100% of 17,916 features usable),
which is why this is measured per layer rather than assumed.

**Follow-through (2026-09-13).** Reporting the number was necessary but
not sufficient — a user still had to *find* a usable compound, from a
free-text box pre-filled with a raw Broad id. Measured against the cached
matrix: 5,272 of 6,790 compounds (77.6%) have **zero** values on the 79
covered lines, so a typed guess failed most of the time, and the ids are
not something anyone recalls. Two changes close this:

- **`GET /datasets/{id}/aux-features`** — a searchable list of features
  that actually have data on *this* dataset, coverage-sorted, carrying
  drug name and MOA. Searches name, id, mechanism **and** target, because
  a researcher looks for "BRAF" or "MTOR inhibitor", not a trade name;
  matching name+id alone returned nothing for exactly those queries. The
  `usable_only` default hides features below the correlation's own
  minimum, which is the difference between a 6,790-entry list that mostly
  fails and a 1,518-entry list that works.
- **`MIN_AUX_SAMPLES`** is now defined once in `contract.py` and shared by
  the correlation, the usable-feature count, and the picker. These three
  previously each hardcoded `3`; a picker offering features the
  correlation would then reject is worse than no picker.

The zero-coverage error was also reworded. "Need at least 3 to correlate"
frames a never-assayed compound as a threshold the user narrowly missed;
it now says the compound was never assayed on this cohort and that no
gene choice will change it — true for 77.6% of the layer.

### 8.5 Analyses over aux layers
**Implemented.** `backend/app/analysis/aux_analysis.py`, surfaced as a
new "Multi-omics" pane:

- `expression_by_mutation_status` — expression split by mutant vs
  wild-type (GEPIA3's "Hotspot Mutation" equivalent). Unlike gene-vs-gene
  correlation, the same gene on both sides is a *legitimate* question
  ("do NOTCH1-mutant samples express NOTCH1 differently?") and is
  deliberately not blocked.
- `most_mutated_genes` — what to even split on, since a user can't guess
  which genes are recurrently mutated in a cohort.
- `expression_vs_aux_correlation` — expression against a CRISPR
  dependency or drug sensitivity on the same samples.

Every one reports `n_dataset_total` / `n_excluded` / `exclusion_reason`.
The UI disables (rather than hides) tabs whose layer a dataset lacks, and
auto-selects a tab the dataset can actually answer.

### 8.6 Exclusion reasons must distinguish *why* a sample is absent
**Decided (2026-09-13), after QA.** `expression_vs_aux_correlation`
reports `n_missing_layer` and `n_missing_feature` separately, not one
merged `n_excluded`.

A sample can be absent for two unrelated reasons: it has no aux layer at
all, or it has the layer but this particular feature was never assayed on
it. Merging them produced a UI that contradicted itself on the CRISPR
tab's own default value — "covers 91 of 186 samples (49%)" directly above
"173 excluded (no crispr gene effect data)", when in fact 95 lacked the
layer and **78 had it** but weren't screened for RPL13A (which DepMap
scores on only 349 of 1,178 lines). Blaming all 173 on a missing layer
misattributes 78 samples, which is precisely what §6.1's transparency
requirement exists to prevent — so a merged count wasn't a cosmetic
shortcut, it was the failure mode in miniature.

Two related QA fixes, same principle:
- A gene present in the mutation matrix but mutated in **zero** samples
  that also have expression data now raises rather than returning a
  degenerate 200 with `n_mutated=0, p_value=null`. 2,701 of TARGET's
  6,250 mutated genes are in that state.
- `_resolve_gene()` now normalises case/whitespace centrally. One
  endpoint was doing `.strip().upper()` on its own parameter while the
  shared resolver wasn't, so `?gene=dtx1` 404'd on an API the UI happens
  to uppercase for you.

### 8.6a One RNA-seq file per TARGET-ALL-P2 case, chosen deterministically
**Decided (2026-09-22), after independent verification during QA.**
`select_one_file_per_case()` in `app/ingest/gdc_target.py` picks exactly
one GDC file per case (patient) before the expression matrix is
assembled, and the same selection is applied to the mutation aux layer
(`app/ingest/target_mutations.py`).

**The bug.** `list_open_rna_files()` queries GDC with no `sample_type`
or `workflow_type` filter and no dedup. 56 of the 469 TARGET-ALL-P2
cases have more than one open-access RNA-seq file: 48 carry both a
Primary and a Recurrent sample (a genuine second disease timepoint —
relapse tissue, not a technical duplicate), 8 carry two files of the
*identical* sample_type (a true duplicate upload), and 1 case has only
Recurrent files. The pre-fix code (`if sample_path.exists(): continue`,
keyed on `case.submitter_id`) silently kept whichever file GDC happened
to list first — non-deterministic across ingests, and undetectable from
the outside since the total sample count didn't change. Found by
spot-checking a production reference value against a fresh local
ingest: identical n but a different statistic, traced via direct GDC
`/files` queries to this exact mechanism.

**Why picking one, not exposing both, was the right fix.** Every
grouping/clinical column (`vital_status`, `etp_status`, `mrd_status`)
and the mutation aux layer are keyed at the same patient level as the
expression matrix. `target_all_p2` is also the only dataset in this app
with `supports_survival=True`. Exposing a patient's Primary and
Recurrent samples as two separate rows would give that patient the same
`vital_status`/follow-up time twice — double-counting their survival
event in Kaplan-Meier/Cox, which assume one independent observation per
row — and double their weight in any group comparison or
differential-expression run that happens to include them. It would also
break the "469 samples = 469 patients" equivalence assumed everywhere
in the UI and in this document, for a subset of the cohort, silently.
One row per patient is a correctness requirement here, not a
simplification.

**Selection order:** Primary tumor over Recurrent (baseline disease
state, matching how published TARGET-ALL-P2 analyses treat this
cohort) → Bone Marrow over Peripheral Blood (the dominant, best
represented specimen site in this project) → lowest `file_id` (breaks
the 8 true-duplicate ties the same way on every ingest). 55 of the 56
affected cases resolve via the Primary preference; the 1 Recurrent-only
case falls back correctly rather than dropping out of the cohort.

Verified against live GDC before writing the fix: re-queried `/files`
directly and confirmed the exact partition (55 cases with an available
Primary sample, 1 without, 8 same-type duplicates). `select_one_file_per_case()`
is unit-tested against all three cases plus the ordinary single-file
case (`tests/test_gdc_target.py`), each validated by reintroducing the
bug and confirming the test fails first.

**Also fixed in the same change:** `scripts/ingest_dataset.py` was
re-reading its own previously-ingested database row instead of the
source loader on every subsequent run, once a dataset had been ingested
once — the database-wins precedence that's correct for serving requests
(§10.5) made "re-ingest" a no-op against the actual GDC source after the
first run. This is exactly why the fix above didn't reach data on disk
until this was also fixed: registry now exposes
`get_code_descriptor()`, which the ingest script uses instead, and which
always reaches the real loader regardless of what's already in the
database. Covered by two new tests in `test_registry_db.py`.

### 3.1a Tied scores at a percentile cutoff
**Decided (2026-09-13), after QA.** `binarize_by_cutoff()` takes a
`tie_policy` (`trim` default / `exclude` / `inclusive`), and the survival
response carries a `tie_note` reporting how many samples were tied and
which policy ran.

A percentile threshold can land exactly on a block of equal scores, and
the original inclusive comparisons (`scores <= low_threshold`) then swept
the **whole block** into one arm. This is not an edge case here:
`auc_signature_score` returns exactly `0.0` for every sample where the
gene never reaches the top-25% rank cut (§2.3 — a score of 0 is
meaningful, not missing), so for a single gene expressed in under ~75% of
the cohort the 25th percentile *is* 0.0. Measured on TARGET:

| gene | scores = 0.0 | inclusive | trim | exclude |
|---|---|---|---|---|
| MEF2C | 217 / 466 | 118 / **217** | 116 / 116 | 116 / **0** |
| TAL1 | 242 / 466 | 117 / **242** | 117 / 116 | 117 / **0** |
| LMO2 | 34 / 466 | 117 / 117 | 117 / 117 | 117 / 117 |

The UI meanwhile asserted "top vs. bottom quartile (middle 50%
excluded)", so the log-rank p and Cox HR were computed on arms that were
not what the caption claimed. For MEF2C the log-rank p moves from
7.51e-13 (inclusive) to 1.08e-08 (trim) — same data, materially
different strength of conclusion.

**Why this is a user-facing option rather than an internal default.**
There is no policy that is simply correct: `trim` matches the stated
definition but picks arbitrarily among equals; `exclude` never splits
equal scores across arms but **empties an arm entirely** for MEF2C and
TAL1, making the analysis impossible for exactly the genes that trigger
it; `inclusive` reproduces prior numbers. Choosing silently would hide a
judgement that belongs to the domain expert running the analysis — so
the app surfaces the tie, names the applied policy, prints the resulting
arm sizes beside the curve, and lets the user switch. LMO2 is identical
under all three policies, confirming the handling engages only where
ties actually occur.

`trim` is the default because it is the only policy that keeps the arms
consistent with what the interface says it is doing, and it never
silently produces an unusable result.

**Addendum (2026-09-23), after QA.** An emptied arm under `exclude` was
disclosed in the KM plot's tie note, but the Cox table below it was not
gated on the same condition — `lifelines` still fits a coefficient when
one arm is empty, and the result is a hazard ratio and 95% CI spanning
many orders of magnitude (observed live: target_all_p2/MEF2C, quartile +
exclude → HR=18341.107, CI 0.000–688692843566556.625, p=0.4294) rendered
in the same table styling as a real result, with no signal that the
model has nothing to compare. Fixed in `SurvivalPage.tsx`: the Cox
table and its summary p-value are now suppressed whenever fewer than 2
arms are present, replaced with a message naming the empty-arm cause and
pointing at the control that fixes it. The KM curve, tie note, and
sample-accounting line are unaffected — this only changes the Cox
model's ability to fabricate a data-free statistic. 2 new frontend
tests, validated by reintroducing the bug and confirming both fail.

### 8.7 Aux feature ids normalise like gene symbols — but not by upper-casing
**Decided (2026-09-13), after QA.** `_resolve_aux_feature()` mirrors
`_resolve_gene()`'s normalisation contract for an aux layer's own row
labels, but resolves by *exact match first, case-folded fallback second*
rather than by `.upper()`.

§8.6 fixed case handling for `gene` and missed `aux_feature` sitting in
the same endpoint, so one pane disagreed with itself: `?gene=myb`
resolved while `?aux_feature=myb` returned 404 *"myb not found in the
crispr_gene_effect layer"*. That phrasing is the real damage — it reads
as a claim about **data coverage**, so a user could reasonably conclude
MYB was never screened and abandon a valid analysis, when the true answer
is r=-0.65, p=2.5e-12 across 91 lines.

Upper-casing (the rule that is correct for gene symbols) is **not** safe
here: the drug layer is keyed by Broad compound id
(`BRD:BRD-A00047421-001-01-9`), not symbol. Hence exact-match-first, and
an explicit 422 rather than an arbitrary pick if a layer ever carried two
labels differing only in case. Verified zero case-collisions across all
three layers today (17,916 CRISPR + 6,250 mutation + 6,790 drug
features), so the fallback is unambiguous wherever it currently applies.

### 8.8 Gene annotation tolerates mygene.info's polymorphic `ensembl` field
**Decided (2026-09-13), after QA.** `GET /genes/{symbol}` crashed with a
500 (`AttributeError: 'list' object has no attribute 'get'`) for any
symbol mapping to multiple Ensembl genes.

mygene.info returns `ensembl` as a dict for most genes but as a **list of
dicts** when a symbol has several Ensembl mappings. Only 1 of 30 surveyed
T-ALL genes (ZAP70) has that shape — which is exactly why sampling missed
it — but ZAP70 is a core T-cell-receptor signalling gene a T-ALL
researcher would plausibly look up. Where several exist we take the first
stable `gene` id rather than inventing a merge: callers want one id to
link out with, and the alternatives are transcript-level variants of the
same locus.

Found by tailing the server log during the UI sweep, **not** from the
browser: the annotation panel degraded silently, so nothing looked wrong
on screen. Worth recording as a testing lesson — a UI-only pass cannot
see a 500 that the frontend swallows. The shaper now also degrades to
`None` for absent/empty/scalar shapes instead of raising or silently
returning a wrong value.

---

## 9. Testing

### 9.1 The API surface is tested at the HTTP layer
**Decided (2026-09-13).** `backend/tests/test_api.py` drives the app
through FastAPI's `TestClient` against a synthetic dataset registered
into the registry — no network, no real cohort load.

The rest of the suite tests analysis functions directly, which is the
right level for the statistics but left `app/api/main.py` at **26%
coverage** across 20 endpoints. That gap was not theoretical: *every*
bug found in that file during this session's QA — a case-sensitive
`aux_feature` beside a case-insensitive `gene`, a 500 on mygene.info's
list-shaped `ensembl`, `/group-values` returning an empty list instead
of a 404, `/correlation` crashing on identical genes — was an HTTP-shape
bug where the underlying analysis function was correct and the endpoint
around it was not. Testing one layer down could not have caught any of
them, and in fact did not: all were found by a human or an agent
clicking the UI.

Coverage of `app/api/main.py` went 26% → **74%**.

The tests were validated by reintroducing two of the real bugs and
confirming the suite fails — a test that passes without being able to
fail is worse than no test, because it reads as coverage. Reinstating
the case-sensitivity bug failed exactly the two case-variant cases;
reordering `/datasets/compare-multi` below the dynamic
`/datasets/{dataset_id}` route failed the shadowing test. Both restored
cleanly.

### 9.1a Tests cannot reach the network
**Decided (2026-11-15).** `backend/tests/conftest.py` refuses every
non-loopback connection and DNS lookup for every test.

CLAUDE.md already said tests must not hit the network; nothing enforced
it. `GET /datasets` loads every registered dataset that is not warming
and swallows load failures by design (§7.9), so
`test_dataset_listing_includes_provenance` downloaded every real cohort
whenever the machine had network access: ~1,150 s of a ~1,230 s CI run
(measured with Python 3.12 and an empty cache), with no failure to show
for it. Main's CI had been taking 10–20 minutes for this reason. With the
guard, those loaders fail in milliseconds, the endpoint's behaviour under
test is unchanged, and the full suite runs in ~7 s on an empty cache.
`tests/test_network_guard.py` proves the guard refuses outbound HTTP and
raw sockets and keeps loopback working (it fails with the guard
disabled). A future test that needs real data must use a vendored
fixture, as the contrast tests do.

### 9.2 What is deliberately not covered
- **Ingest loaders** (33–67%) hit the network and download large files.
  Their parsing logic is covered by `test_aux_ingest.py` against
  fixtures; the download paths are exercised by actually running them,
  which is how the mutation matrix was validated (§8.2).
- **Most of the frontend.** See §9.3 — there is now a component suite,
  but it is deliberately narrow (~30% of statements). Page scaffolding,
  dockview wiring, and chart rendering are still verified by
  browser-driving QA rather than unit tests.

### 9.3 The frontend suite targets shipped defects, not a coverage number
**Decided (2026-09-13).** Vitest + Testing Library, run with
`npm test` in `frontend/`. Every test corresponds to a bug that actually
reached the app, rather than to a coverage target:

- **`GeneAnnotation`** — the LYL1 white-screen (§8.8). Covers the bare
  string alias, normal lists, empty lists, four malformed shapes, and the
  fail-quiet path.
- **`histogramBinCount`** — the bin rule (§5.6), pinned against the real
  cohort sizes (52 / 186 / 469) and the zero-inflated shape of §3.1a.
  Extracted from an inline IIFE to make it testable, which also made it
  self-documenting.
- **`DatasetPanel`** — export identity (§5.7).

**Each test was validated by reintroducing its bug and confirming the
suite fails.** That check earned its keep immediately: the export tests
were *first* written against `ExportButton` alone, passed, and then
**failed to notice** when the real fix was reverted — because they only
proved the button forwards its props, while the actual defect was the
call site not passing the gene. Rewriting them against `DatasetPanel`
made them fail correctly. A test that cannot fail is worse than no test,
since it reads as coverage.

Coverage is ~30% of frontend statements and that is the intended shape:
the components that carried real defects are covered; the rest is page
scaffolding and chart rendering, where a browser QA pass finds more than
a jsdom assertion would. jsdom needs `ResizeObserver` and a non-zero
`clientWidth` stubbed (`src/test/setup.ts`) or any chart render throws.

---

## 10. Dataset metadata database

### 10.1 What goes in the database, and what deliberately does not
**Decided (2026-09-20).** PostgreSQL holds dataset *metadata*; the
expression values stay in parquet files. `app/db/models.py`.

The split is measured, not stylistic:

| matrix | shape | cells |
|---|---|---|
| TARGET expression | 60,660 x 469 | 28.4M |
| DepMap expression | 19,193 x 1,673 | 32.1M |
| DepMap CRISPR | 17,916 x 1,178 | 21.1M |
| PRISM drug | 6,790 x 919 | 6.2M |
| TARGET mutations | 6,250 x 717 | 4.5M |

92M dense numeric cells today. Every analysis reads whole rows or whole
columns of them -- a genome-wide differential scan touches all 60,660
genes and costs **0.03s** against parquet (measured). Stored relationally
that is 92M rows now and 600M+ as datasets are added, turning each scan
into a 60,000-row aggregate. Row storage buys nothing for dense numeric
data whose access pattern is "read it all", and costs an order of
magnitude. The metadata, by contrast, is ~102,000 small rows that are
genuinely filtered, grouped and joined -- exactly what SQL is for.

### 10.2 Why a database at all
Dataset identity was hardcoded: `registry.ensure_loaded()` held a literal
import list and each `ingest/*.py` carried its own accession and group
columns. Adding a dataset therefore required editing code and
redeploying. Moving the description into a table makes dataset addition a
data operation (`scripts/ingest_dataset.py`) instead of a code change --
which is the entire point, and the only reason this was worth doing.

The loader stays in code. Parsing a new upstream format is genuinely
programming, and pretending otherwise would move that work somewhere
worse.

### 10.3 Sample attributes are JSON, deliberately
The cohorts really do differ: TARGET carries `vital_status` /
`days_to_death` / `etp_status` / `mrd_status`, DepMap carries `lineage` /
`subtype` / `cell_line_name`, GDS4299 carries only `etp_status`. Typed
columns would mean a schema migration for every new cohort, defeating the
purpose of the change.

The tradeoff is real and worth stating plainly: column-level type
checking and NOT NULL constraints are given up. Postgres still indexes
and filters inside JSONB, so query capability is retained -- verified:
"ETP samples across all cohorts" returns gds4299=12, target_all_p2=19 in
2ms, a question that previously required loading ~500MB of matrices.
Validation moves to ingest time rather than write time.

### 10.5 The registry reads from the database, with code as the floor
**Decided (2026-09-21).** `registry.list_descriptors()` merges rows from
the `datasets` table (status `ready`) over the code-registered ones, so a
dataset that exists only as a row plus a parquet file is served without
any code change. That is the cutover that actually delivers "add a
dataset without a redeploy" -- the schema alone did not.

The descriptor's `loader` re-fetches its row by id rather than closing
over the ORM object: the session is closed by the time the loader runs,
and a fresh fetch also means a re-ingest is picked up on the next load
instead of serving a stale description.

**Code-registered datasets remain the fallback, deliberately.** An
unreachable, unmigrated, or empty database degrades to exactly the
previous behaviour and logs why -- verified by pointing `DATABASE_URL` at
a dead host: `/health`, `/datasets` and a real comparison all still
returned 200 off the code path. The previous behaviour is the floor,
never the casualty. Where both exist the database wins, since that is the
row an operator most recently wrote.

Two round-trip details that would break grouping silently if missed, both
covered by tests:
- Ingest flattens clinical fields into one JSON blob, but the analysis
  layer reads them **both** as top-level columns and from a nested
  `group_columns` dict. Both shapes are restored.
- Sample rows come back in arbitrary order, so they are re-ordered to
  follow the matrix columns -- the analysis layer aligns the two
  positionally.

Verified after cutover: a dataset loaded from the database is
bit-identical to the same dataset loaded through its ingest module (same
matrix values, same sample order, same 12 ETP / 40 non-ETP counts), and
every reference statistic is unchanged to 12 decimal places -- MYB CRISPR
r=-0.652360635454 n=91, DTX1/NOTCH1 p=4.940894660436e-03 (159/120/190),
quartile survival arms 116/116 with 217 tied. 165 tests pass on both
SQLite and Postgres.

### 10.4 Portability and ordering
`models.py` uses SQLAlchemy's generic `JSON` type rather than
Postgres-native `JSONB`, so the same models run on SQLite locally and in
tests (zero-setup, matching how `CACHE_DIR` and `CORS_ORIGINS` already
behave) and on Postgres in production. Alembic owns schema changes and
reads `DATABASE_URL` through the same resolver the app uses, so
migrations cannot drift to a different database than the one being
served.

Ingest writes matrices to their final location **before** marking a row
`ready`, so a crash leaves a `draft`/`failed` row rather than a `ready`
row pointing at a missing file. A failed ingest is recorded as `failed`
with its error rather than leaving the dataset silently absent.

**Verified against a live Postgres 16 instance (2026-09-21).** The
Alembic migration applied cleanly, all three datasets ingested in 17s
(707 samples / 101,449 features / 3 aux layers -- identical counts to
SQLite), the full suite passed with `DATABASE_URL` pointed at Postgres,
and the app served the unchanged reference value through it (MYB CRISPR
r=-0.652361, n=91). Cross-cohort JSON filtering behaves as designed:
`attributes->>'etp_status' = 'ETP'` returns gds4299=12, target_all_p2=19.
Re-ingest replaces rather than duplicates (cascade delete confirmed on
Postgres, not just SQLite).

Local compose maps Postgres to host port **5433**, not 5432 -- the
default is commonly taken by another local Postgres. The container still
listens on 5432 internally, so `DATABASE_URL` inside compose is
unaffected.

---

## 11. Designed perturbation experiments

Added 2026-11-15 for GSE324984 (MYCL gain and loss in prostate cancer
cell lines), the first dataset in this app that is a designed experiment
rather than a cohort.

### 11.1 Contrasts are part of the dataset, and stay within a cell line
**Decided.** A dataset may declare pre-specified two-group contrasts
(`ContrastSpec` in `app/models/contract.py`), computed at ingest and
stored (`contrasts` table + one parquet per contrast; alembic revision
`b3c1f0a2d4e5`). GSE324984 declares two:

| contrast | reference (n) | treatment (n) |
|---|---|---|
| `c42b_mycl_oe` | C4-2B empty vector (2) | C4-2B MYCL overexpression (2) |
| `h660_mycl_kd` | NCI-H660 shCtrl (2) | NCI-H660 shMYCL (2) |

Why fixed contrasts rather than "pick any two groups": the cell lines
differ far more than any perturbation changes them (log2-TPM r ≈ 0.99
between replicates, ≈ 0.987 within a line across arms, ≈ 0.89 across
lines). A C4-2B control against an H660 knockdown measures mostly
"adenocarcinoma vs NEPC line". Each contrast is fitted on its own four
samples only — filtering, normalisation and the variance prior included —
so one experiment cannot influence the other's numbers (pinned by
`test_other_experiments_samples_do_not_influence_a_contrast`).

The rest of the app still sees all 8 samples (Compare by `arm` shows the
four arms separately; Rank across all 8 is legitimate). The dataset notes
warn that correlations or PCA across all 8 samples mostly reflect cell
line. On Compare, every pairwise test whose two groups lie in different
strata of a contrast's `stratum_column` is marked `crosses_stratum`, and
the response carries a one-line `stratum_caveat` (GSE324984 by `arm`: 4
of the 6 pairs cross cell lines) — those p-values were previously shown
in one table with the within-line pairs and nothing to tell them apart.

### 11.2 Method: edgeR + limma-trend on the counts
**Decided.** `app/analysis/contrast.py` → `run_contrast()`.

    filterByExpr(group)  →  TMM  →  cpm(log=TRUE, prior.count=3)
    →  lmFit(~group)  →  eBayes(trend=TRUE)  →  BH over tested genes

Why not the app's Mann–Whitney: at 2 vs 2 it cannot go below p = 0.333
(§1.4a), so every gene would read "not significant" regardless of the
data. Replicated perturbation RNA-seq needs a test that shares variance
information across genes; limma's empirical-Bayes moderated t is a
standard one (prior df here: 6.97 for C4-2B, 13.75 for H660).

**It is a re-implementation, verified against R.** Each step mirrors the
named R function with its default arguments. On the full 78,691-gene
counts file, against R 4.2 / edgeR 3.40.2 / limma 3.54.2: identical
tested gene sets (16,855 and 16,987), identical TMM factors to 6 dp,
identical prior df, and every log2FC / AveExpr / t / p / q within
2e-12 absolute; 5,905 and 4,224 genes at q < 0.05 in both
implementations. A 2,720-gene subset with R's output is vendored in
`tests/fixtures/contrast/` so the comparison runs in the test suite
offline (`tests/test_contrast.py`), and the suite was checked to fail on
seven planted deviations (prior count, trend off, unweighted TMM, BH
off-by-one, using other experiments' samples, keeping pre-filter library
sizes, no filter).

**Edge cases matched to limma, found in review (2026-11-15):** when the
variance prior is infinite (very regular replicates) limma uses the plain
pooled variance as the prior scale without a trend, `exp(emean)` with
one — the first was wrong here and is now checked against R on a Poisson
fixture for both `trend=TRUE` and `FALSE`; an infinite prior df is stored
as `df_prior: null, df_prior_infinite: true` because JSON (and Postgres
JSON) cannot hold Infinity, which would have made every contrast endpoint
500; TMM returns factor 1 for a sample with no finite log-ratios, as R
does, instead of raising.

**Choices stated:** counts used as submitted (11% fractional, which edgeR
accepts — no rounding); genes failing `filterByExpr` are kept in the
table as `tested = false` and shown as "not tested", never as 0 or as
absent; q-values are genome-wide within the experiment, so the gene
panel and the Top-genes list always show the same q for the same gene.
A panel-only correction (25 genes) would give smaller q-values and is
defensible for a pre-specified list, but would mean two different q's for
one gene in one app — not done.

**Why not DESeq2:** it is the other standard choice and likely what the
authors used, so log2FC/p will differ somewhat from the paper's tables.
Adding it means a new Python dependency (`pydeseq2`); installing it into
the project environment was declined by this session's permission
settings, so it was not added. Implemented with numpy/scipy only, the
dataset ingests from the existing requirements. Revisit if the expert
wants the paper's own numbers reproduced.

**Checked against DESeq2 (R 4.2, DESeq2 1.38.3, counts rounded, default
Wald test, fitted per cell line), validation only — 2026-11-15.**
Genome-wide, the two agree on effect sizes: log2FC Pearson r = 0.966
(C4-2B) and 0.963 (H660) over all tested genes; signs agree for 98.7% /
98.4% of genes and for 100% of genes this app calls q < 0.05. limma-trend
is the more liberal of the two at n = 2: 5,905 vs 4,355 genes at q < 0.05
in C4-2B (overlap 4,316) and 4,224 vs 2,717 in H660 (overlap 2,635). For
the PP2A panel the headline is method-independent — PPP2R5B (+1.26 /
−0.89) and MYC (−0.45 / +0.79) are reciprocal and significant in both
experiments under both methods; CIP2A, SET, ANP32A, PPP2R1B (C4-2B) and
PPP2R2C, PPP2R5D (H660) are significant under both. Borderline genes
differ: PPP2R2A (q 0.034 vs padj 0.067) and PPP2R2C (0.016 vs 0.062) in
C4-2B are significant only here; ARPP19 in C4-2B (0.061 vs 0.048) and SET
in H660 (0.126 vs 0.014) only under DESeq2. A reader should treat those
four as not robust to the choice of method.

### 11.3 GSE324984 source-file problems, and what was done
**Decided.** `app/ingest/gse324984.py`.

- **Truncated TPM file.** `GSE324984_RNAseq_TPM_matrix.txt.gz` has
  exactly 65,535 gene rows — the legacy Excel row limit — and is a strict
  prefix of the 78,691-row counts file. The missing 13,156 rows are the
  tail of the Ensembl ordering (ENSG00000297233 onward). TPM values that
  are present are unchanged (each column sums to ≈ 997–998k, i.e. TPM
  computed before the cut). So: expression shown = submitted TPM for
  65,535 genes; statistics = complete counts. Top-genes rows with no TPM
  are flagged `in_expression_matrix = false`. All 23 panel genes are in
  both files.
- **One gene mis-named in the counts file.** ENSG00000073282 is "TP63" in
  the TPM file (and Ensembl) but "KDM5A" in the counts file; real KDM5A
  (ENSG00000073614) is present separately in both, and the row's counts
  match TP63's near-zero TPM. The Ensembl id is the key; the TPM file's
  name wins. More than 5 such conflicts fails ingest as misaligned files.
- **TPM and counts are different quantifications.** If the TPM were the
  counts divided by a fixed gene length, log(counts / TPM) would be
  exactly gene effect + sample effect. It is not: residual SD 0.066 over
  10,227 well-expressed genes, with 7.4% deviating by > 20% in some
  sample — consistent with an isoform-aware TPM (sample-specific
  effective lengths) beside gene-level counts. Consequence: for 2,029 of
  16,476 judged C4-2B genes (12%) and 910 of 16,699 H660 genes (5%) the
  mean TPM moves the opposite way from the count-based log2FC (largest
  case at q < 0.05: HSPA12A, log2FC +1.04, TPM 4.8 → 3.6). Neither is
  "wrong"; they measure differently. The test is the count-based one, so
  each contrast stores the TMM-normalised CPM it ran on
  (`norm_cpm:<sample>` columns, checked against R's library sizes ×
  TMM factors), the panel returns those alongside TPM, and
  `display_direction_agrees = false` flags the rows where the two
  disagree. All tested PP2A-panel genes agree.
- **Duplicate names — one locus per name, chosen over the complete
  counts (revised 2026-11-15 after review).** 1,595 TPM rows share a name
  with another Ensembl id: Y_RNA/snoRNA families, and PAR X/Y copies
  (GTPBP6, CD99, PLCXD1…) whose reads the aligner splits arbitrarily
  between copies — the X and Y copy of one gene showed opposite-signed
  fold-changes at q < 0.05. For each name one locus keeps the name: the
  highest mean TPM when every locus with that name has TPM (TPM is
  length-normalised: C4orf36's loci have equal mean counts, 24.7 vs 25.0,
  but TPM 1.77 vs 0.71), otherwise the highest mean count (TPM cannot rank
  a locus it lacks). The others show their Ensembl id, carry no alias for
  the shared name, and are flagged (`n_loci_with_name`, an identity note
  on panel rows). The first version chose by TPM among TPM rows only,
  which (a) sent "CAST" to a near-silent TPM-file locus while the
  expressed CAST (ENSG00000310517, ~1,000–5,000 counts) exists only past
  the truncation, (b) let the top-genes list (undeduplicated) and the
  panel show different loci under one symbol, and (c) stored the shared
  name as an *alias*, so the legacy-name fallback resolved "C4orf36" to
  the secondary locus. All three fixed and pinned by tests built on those
  cases. Label-independent and deterministic either way.
- **Genes only in the counts are reachable.** 13,156 genes have
  statistics but no TPM row (356 tested in C4-2B, 32 at q < 0.05; 265 in
  H660, 44). The gene panel resolves them by id or symbol and says the TPM
  column is empty and why; Compare/Rank (TPM-based) return a 404 that
  names the locus and points to the contrasts, rather than "not found".
- **Label check on both files.** The MYCL check now runs on the counts
  (as CPM) too — a column swap in the counts file alone would otherwise
  reach every statistic unchecked. Real files: MYCL counts 371/320 in
  controls vs 65,524/70,252 with overexpression.
- **Column mapping.** Matrix columns are library names
  (`C4-2B_Ctrl_1`, `shMYCL_H660_2`, …) in a different order from the
  GSMs; they are mapped through each GSM's "Library name" description,
  never by position.
- **Label check.** Ingest fails unless MYCL itself moves as the arms
  say. Measured mean TPM: C4-2B empty vector 5.9 → MYCL OE 1,399.1;
  H660 shCtrl 56.2 → shMYCL 6.6 (≈ 88% knockdown).

### 11.4 The gene panel endpoint
**Decided.** `POST /datasets/{id}/contrast-panel`, `GET
/datasets/{id}/contrasts`, `GET /datasets/{id}/contrasts/{cid}/top`.
Reads stored results only; computes no statistics at request time.

Every contrast summary carries the §6.1 sample accounting —
`n_dataset_total` 8, `n_excluded` 4, and the by-design reason (only that
cell line's two arms enter) — surfaced in the UI.

Unlike a signature score (§2.5), a panel row does not depend on the other
genes, so an unknown gene is listed in `missing` with its reason rather
than failing the request or being dropped. A legacy symbol the dataset
records (KIAA1524, FAM122A) resolves to the current gene with an
`alias_note` on the row (§7.4 update). Each row carries both replicate
values per arm so the reader sees what the fold-change rests on.

### 11.5 What the PP2A panel can and cannot say
**Provisional — for the domain expert.** Two replicates per arm, two
cell lines, one perturbation each. A gene that moves in opposite
directions with overexpression and knockdown, at q < 0.05 in both, is
the pattern MYCL-dependent regulation would produce in these lines; it
is an association, not evidence of direct regulation, and says nothing
about protein level or PP2A activity. Direction-only differences of
small size (|log2FC| < 0.3) at n = 2 deserve caution even when q is
small, since the variance prior is doing much of the work.

The panel's "Across both experiments" column states the direction and
which experiment reaches q < 0.05, by name: "Reciprocal (both q<0.05)",
"Opposite direction; q<0.05 only in C4-2B", "Same direction (both
q<0.05)", "… q≥0.05 in both", or "Tested only in NCI-H660" when the
other experiment's filter dropped the gene. Same-direction movement in
an overexpression and a knockdown is not the MYCL-dependent pattern,
whatever its q. Code: `describeAcross` in
`frontend/src/pages/ContrastsPage.tsx`; the panel CSV's `across_both_experiments` column
uses the same words.

---

## Change log

| Date | Change |
|---|---|
| 2026-08-25 | Document created. Existing decisions (§1–6) recorded retrospectively from the implementation and verified against source. |
| 2026-08-25 | §7.1 probe collapse → highest mean expression; §7.2 cross-dataset → side-by-side, no pooling. |
| 2026-08-25 | §6.1 sample accounting added, prompted by expert review of the n=466 survival count. |
| 2026-08-25 | §3.4 Cox hazard-ratio 95% CI now reported (was computed then dropped). |
| 2026-08-25 | GDS4299/GSE28703 ingested (`app/ingest/gds4299.py`). 21,596 genes × 52 samples (12 ETP / 40 non-ETP). Probe collapse and ETP labels verified — see §7.1, §7.1a. |
| 2026-08-27 | §7.2 built: `GET /datasets/compare-multi` + Compare page multi-select. Default = every dataset exposing the chosen grouping column; user can narrow. Compatibility gates on grouping column only, not assay type, since results are never pooled. Verified against known-good GDS4299 MYCN numbers through the new endpoint. |
| 2026-08-29 | §7.5 added: DepMap loader confirmed stuck on 24Q4 Public (Figshare has nothing newer; DepMap's own portal API is Cloudflare-gated, verified directly). Fix is disclosure, not access — `_release_age_warning()` surfaces staleness in `notes` and logs, doesn't attempt to bypass the bot gate. |
| 2026-08-29 | QA pass (3 parallel agents) on the §7.2 multi-dataset Compare feature + full-app sweep. Found and fixed: `exclusion_reason` wrongly asserted at 0-excluded (both `/compare` and `/compare-multi`); DepMap's §7.5 staleness note computed but never rendered anywhere in the UI (now threaded through `/datasets` → sidebar); Cox 95% CI computed and typed but never rendered on the Survival page (now a table column) — see §3.4 correction. |
| 2026-08-29 | Cheap-wins batch 1 shipped: §3.1 survival cutoff generalised to median/quartile/custom (`binarize_by_cutoff`), quartile/custom deliberately excluding a middle band with its own exclusion accounting; §7.4 partially resolved via a MyGene.info-backed gene annotation panel (`GET /genes/{symbol}`) — resolves aliases for display, explicitly NOT wired into gene-lookup-for-analysis yet (that's a separate, still-open half of §7.4). |
| 2026-08-29 | Cheap-wins batch 2 shipped: §7.6 gene-vs-gene correlation + co-expression network, new Correlation pane. Verified against real biology (MYCN/RAG1 anti-correlation on GDS4299, MYCN/MYCNOS co-expression on TARGET) as well as unit tests; confirmed symbol-mapping works correctly on both an Ensembl-keyed and a symbol-keyed dataset. |
| 2026-08-30 | QA pass (2 parallel agents) on §7.6. Backend logic/math fully confirmed via independent from-scratch recomputation. Found and fixed two UI bugs: results wiped on Correlation↔Co-expression tab switch (conditional-render unmount — now both tabs stay mounted, toggled via CSS); "Kendall" label clipped in the method toggle at ~700-900px pane widths (now wraps instead of clipping). |
| 2026-08-30 | Cheap-wins batch 3 shipped (final batch): §7.7 PCA (SVD-based, no scikit-learn, cross-checked against sklearn at test time) as a third Correlation-pane tab; §7.8 genome-wide top-differential-genes (vectorised Mann-Whitney, ~15x faster than a per-gene loop, verified bit-identical first) as a second Expression Compare tab. New `GET /datasets/{id}/group-values` endpoint added to support the differential-genes UI's group picker. Both tabs applied the mount-persistence pattern from the §7.6 QA fix proactively rather than re-discovering the bug. |
| 2026-08-30 | Final pre-ship QA (4 parallel agents, everything together). Found and fixed: `/correlation` crashed with an unhandled 500 when gene_a == gene_b (duplicate column names after transpose → scipy got a DataFrame, not a Series); `/group-values` silently returned an empty list for an unknown column instead of a 404. Everything else verified to full precision. Shipped as 8220f1d. |
| 2026-09-13 | §8 added: three new data pipelines as auxiliary measurement layers on existing datasets (not new datasets) — TARGET somatic mutations (279/469 coverage, validated against the canonical T-ALL driver profile), DepMap CRISPR gene effect (91/186, validated against known essential/non-essential genes), DepMap PRISM drug sensitivity (79/186, separate Figshare release, transposed matrix, compound→drug-name join). New `Dataset.aux` contract, `aux_analysis.py`, four endpoints, and a Multi-omics pane. Aux loading made opt-in after it was found to make existing unit tests hit the network. |
| 2026-09-13 | QA on §8 (2 parallel agents). Blocker check PASSED: all 739 TARGET MAFs re-downloaded from GDC independently and the mutation matrix reproduced by exact set identity (same barcodes, not just counts) across 7 genes / 717 samples; CRISPR sign convention confirmed two ways (essentials negative, 10/10 lineage oncogenes negative); all 6,789 PRISM compound joins verified against the upstream file. Found and fixed two reporting defects — §8.6 exclusion-reason conflation (visible on the CRISPR tab's default view) and the drug layer's 6,790 headline where only 22% of compounds are usable (now reported as `n_usable_features`). No computed statistic was wrong; every reference number reproduced unchanged after the fixes. |
| 2026-09-13 | E2E UI QA pass (4 parallel agents driving a real browser, Playwright MCP). Found and fixed: §8.7 `aux_feature` was case-sensitive while `gene` in the same endpoint was not, with an error that misread as a data-coverage claim; §8.8 `/genes/{symbol}` 500'd on mygene.info's list-shaped `ensembl` field (ZAP70). Six regression tests added (109 pass). Independently re-verified and found CORRECT, no change needed: DTX1/NOTCH1 biology direction (median 18.86 mutant vs 9.69 WT, p=4.9e-3), CRISPR sign convention (MYB r=-0.652, p=2.5e-12), §8.6 split exclusion reasons (95+78=173 exact), §3.1 quartile wording (25%/25%, 469=335+134), §3.4 Cox CI rendered with HR inside its own CI, §7.2 no cross-metric pooling (tpm and log2_intensity returned side by side, no pooled field). |
| 2026-09-13 | Testing-process note: the canonical backend port is **8420** (`vite.config.ts`, `README.md`, `docker-compose.yml`, `backend/Dockerfile` all agree). Starting uvicorn on any other port makes every frontend call 502 and makes every pane render empty — which looks like a catastrophic app failure and wasted most of one QA agent's run. Start the backend on 8420. |
| 2026-09-13 | §3.1a added after QA found quartile/custom survival arms silently violating their own definition when scores tie at a cutoff (MEF2C LOW n=217 of 466, captioned as a bottom quartile). Exposed as a user-chosen `tie_policy` with a reported `tie_note` rather than an internal default, at the user's direction: the domain expert should be able to see what happened and change it, since no policy is unconditionally correct (`exclude` empties an arm for the very genes that trigger it). §8.8 extended — `alias` has the same list-or-string polymorphism as `ensembl`; LYL1 (single alias "bHLHa18", and a member of this app's own ETP-TF5 preset) returned a bare string against a `list[str]` contract and white-screened the entire workspace via `.join` on a string. Fixed in the API, made structurally impossible in `GeneAnnotation`, and the poisoned LYL1 cache entry cleared. Multi-omics tab a11y: disabled tabs shared one accessible name ("Not available for this dataset") so a screen reader couldn't tell which layer was missing — now names the layer; disabled tabs no longer keep the active tint. 115 tests pass. |
| 2026-09-13 | Non-blocking QA findings cleared: §5.7 Expression Compare PNG export now carries the gene in filename and in-image title (two genes previously overwrote each other as `<dataset>-by-group.png`, with the gene nowhere in the image); PCA export title no longer a bare "PCA". §5.6 Signature Score histogram bins adapt (Freedman–Diaconis with a sqrt(n) floor, clamped 6–36) instead of a fixed 36 — GDS4299 goes from 39% empty bins to 0% while TARGET keeps its resolution at 22 bins. §7.2 `ExpressionUnit.LOG2_TPM` added so DepMap (RNA-seq log2 TPM) stops sharing the microarray "log2 intensity" label with GDS4299. Display/metadata only — no computed statistic changes. 115 tests pass, tsc clean, build clean. |
| 2026-09-13 | §8.4 follow-through on PRISM sparsity: new `GET /datasets/{id}/aux-features` searchable picker over features that actually have data on the selected dataset (searches name/id/MOA/target), replacing a free-text box pre-filled with a raw Broad id — 77.6% of compounds have zero measurements on the covered lines, so typed guesses mostly failed. `MIN_AUX_SAMPLES` centralised in `contract.py` so the picker, the usable-feature count, and the correlation cannot disagree. Zero-coverage error now says the compound was never assayed on this cohort rather than "need at least 3 to correlate". Drug tab default moved to AZ-628 (the one compound covering all 79 lines). 118 tests pass. |
| 2026-09-13 | §9 added. HTTP-level test suite for the API surface (`tests/test_api.py`, 24 tests against a synthetic in-registry dataset — no network). `app/api/main.py` coverage 26% → 74%; suite 118 → 142 tests. Motivated by the observation that every API bug found in this session's QA was an HTTP-shape bug invisible to the existing function-level tests. Suite validated by reintroducing two real bugs (case-sensitive aux_feature; compare-multi shadowed by the dynamic route) and confirming it fails on both. Documented the remaining deliberate gaps: ingest download paths, and the absence of any frontend test suite. |
| 2026-09-13 | §9.3 added: frontend test suite (Vitest + Testing Library, 21 tests, `npm test` in `frontend/`) — the last layer with zero automated coverage. Scoped to defects that actually shipped: the LYL1 white-screen (§8.8), the histogram bin rule (§5.6), and export identity (§5.7). Each validated by reintroducing its bug; the export tests initially passed while failing to catch a reverted fix, because they exercised `ExportButton` rather than the call site that carried the defect — rewritten against `DatasetPanel`. `histogramBinCount` extracted from an inline IIFE and exported to make it testable. Coverage ~30% by intent, not by omission. |
| 2026-09-20 | §7.9 added after a failed Railway deploy. Dataset warm-up moved off the startup path onto a background thread (`app/services/warmup.py`): uvicorn binds immediately, `/health` answers on liveness alone, `/readiness` reports per-dataset state, and `/datasets` returns `warming: true` rather than blocking on a loader. Root cause measured, not guessed — §8's aux layers pushed cold prewarm to ~848s against the 900s healthcheck budget (94% consumed, ~52s headroom from a home connection; GDC is slower to datacenter IPs and the MAF loader retries 4x per file). `healthcheckTimeout` reduced 900→120s because it no longer waits on data. Verified against an empty cache: /health 200 in 2s, /datasets in 16ms mid-download, DepMap ready at 90s serving its reference value while the other two still downloaded. 9 new tests (backend 151, frontend 24), each validated by reverting the fix. Warm cache is ~495MB — the Railway volume must exceed that. |
| 2026-09-20 | §10 added: dataset metadata moved into a database (`app/db/`, 4 tables, Alembic migration) so adding a dataset stops requiring a code change and redeploy. Postgres in production, SQLite locally/in tests via the same generic-typed models. Expression matrices deliberately stay as parquet — measured 92M dense cells across 5 matrices, and a genome-wide scan costs 0.03s against parquet vs a 60,000-row aggregate as SQL. `scripts/ingest_dataset.py` is the new add-a-dataset path; all 3 existing datasets migrated (102,162 metadata rows). Cross-dataset metadata queries now possible without loading matrices — "ETP samples across all cohorts" in 2ms. 7 new tests (158 total). NOT yet run against a live Postgres (Docker unavailable); schema compiles for the Postgres dialect. |
| 2026-09-21 | §10 verified end to end against live Postgres 16 (Docker). Migration applied, 3 datasets ingested in 17s with counts identical to SQLite, 158 tests passed with DATABASE_URL on Postgres, app served the unchanged MYB CRISPR reference value (r=-0.652361, n=91) through it, and JSONB cross-cohort filtering confirmed (ETP: gds4299=12, target_all_p2=19). Re-ingest replace-not-duplicate confirmed. Compose maps host port 5433 since 5432 is commonly occupied. |
| 2026-09-21 | §10.5: registry cut over to read datasets from the database, with code-registered datasets kept as the fallback floor. A dataset that exists only as a row + parquet now loads with no ingest module and no redeploy (covered by a test that seeds exactly that). Verified bit-identical to the code path and every reference statistic unchanged to 12 dp; graceful degradation confirmed by pointing DATABASE_URL at a dead host. 7 new tests, 165 total, passing on both SQLite and Postgres. |
| 2026-09-21 | §7.10 added after a production incident: `target_all_p2` failed to warm on the live Railway deploy (502s, "Excel file format cannot be determined") while `depmap`/`gds4299` warmed fine — one supplementary clinical-supplement fetch with no retry and no response-body validation took the whole 469-sample dataset down. New `app/ingest/_http.py` (`fetch_with_retry` + `validate_excel_bytes`) generalises the existing MAF-download retry pattern and adds response-body sniffing before parsing; wired into both `target_mrd.py` and `liu2017_etp.py` (the latter independently observed hanging 277s on the same class of failure while investigating). `gdc_target.load()` no longer lets either fetch be fatal — `etp_status`/`mrd_status` follow the same additive-never-fatal contract as the §8.1 aux layers, verified by mocking both to fail unconditionally and confirming the dataset still loads with its full matrix. Also closed a related gap: `/group-values` now reports `unavailable: true` for a declared column with zero current values, distinct from an unknown column, surfaced in the UI instead of a silently-empty picker. 9 new backend tests (174 total), frontend unaffected build/typecheck/test all clean. |
| 2026-09-21 | §5.8 added: the Differential Genes tab's "Rank genes" button and group selects went silently disabled with no explanation when the chosen grouping column has fewer than 2 distinct values (e.g. DepMap's `lineage`, filtered to Lymphoid-only upstream) -- the auto-select effect that populates groupA/groupB is gated on >= 2 options, so the pre-existing "must be different" hint (which needs groupA truthy) never fired either. Found during regression testing after the §7.10 deploy fix. Added an explicit message, distinguished from §7.10's `unavailable` case (declared-but-currently-empty) since the two conditions are mutually exclusive by construction. 3 new frontend tests (27 total), each validated by reverting the fix and confirming failure. Verified live: message shows/hides correctly across a single-valued and a two-valued column, no false positives. |
| 2026-09-22 | §8.6a added after independent verification during QA found a real data-accuracy bug: 56 of 469 `target_all_p2` cases had more than one open-access GDC RNA-seq file (48 Primary/Recurrent timepoint pairs, 8 true duplicate uploads, 1 Recurrent-only case), and the loader kept whichever file GDC's API happened to list first -- non-deterministic across ingests, invisible from sample counts alone. Fixed via `select_one_file_per_case()` in `app/ingest/gdc_target.py` (Primary > Recurrent > lowest file_id), applied identically to the mutation aux layer's 22 affected cases in `app/ingest/target_mutations.py`. Also fixed a second bug this uncovered: `scripts/ingest_dataset.py` resolved its loader via the database-wins `get_descriptor()`, so any re-ingest after the first silently re-read its own previous output and never reached GDC again -- added `get_code_descriptor()` to `app/registry.py`, which always bypasses the database, and pointed the ingest script at it. 8 new backend tests (192 total; each new test validated by reintroducing its bug and confirming failure first). Re-ingested from live GDC: 469 samples unchanged (no case lost), zero duplicate sample_ids, and 29 of the 56 previously-ambiguous cases' expression values changed from the prior (buggy) ingest, confirming the fix reached real data, not just the loader logic. Independent code-audit QA pass (separate agent) confirmed no circular-import risk, confirmed `get_code_descriptor()` is fully DB-independent, and found no other ingest module shares this bug pattern; also flagged that the §8.6a changelog entry was missing (added retroactively, this line) and a narrow gap where `ingest_dataset.py --all` would surface a raw `KeyError` for a future database-only dataset with no code loader -- both fixed same-day: `code_descriptor_ids()` added to `app/registry.py` so `--all` now skips and names any database-only dataset instead, covered by a new test. 186 backend tests passing. |
| 2026-09-23 | §3.1a addendum: the Cox table on the Survival page rendered a hazard ratio and 95% CI even when a tie policy emptied one arm entirely -- `lifelines` still fits a coefficient with nothing to compare it against, producing a CI spanning many orders of magnitude (target_all_p2/MEF2C, quartile+exclude: HR=18341.107, CI 0.000-688692843566556.625) styled identically to a real result. Found via live reproduction during the post-§8.6a QA pass. Fixed in `SurvivalPage.tsx`: Cox table and its summary p-value now suppress whenever fewer than 2 curve arms are present, replaced by a message naming the cause and the control that restores both arms. 2 new frontend tests (29 total), validated by reintroducing the bug and confirming both fail first. |
| 2026-11-15 | Two GEO datasets added. **GSE324984** (MYCL overexpression in C4-2B, MYCL knockdown in NCI-H660, n=2 per arm): new §11 designed-experiment contrasts, computed at ingest from the counts with an edgeR+limma-trend re-implementation (`app/analysis/contrast.py`) that matches R limma 3.54 / edgeR 3.40 to ≤2e-12 on all 78,691 genes; new `contrasts` table (alembic `b3c1f0a2d4e5`), `/contrasts`, `/contrast-panel`, `/contrasts/{id}/top` endpoints and a Perturbation Contrasts pane with the PP2A gene preset. Source-file defects found and handled (§11.3): TPM file truncated at 65,535 rows (Excel limit), one gene mis-named in the counts file (TP63/KDM5A), TPM and counts from different quantifications (12%/5% of tested genes disagree in direction, flagged per row). **GSE66586** (neuroblastoma lines, GPL6244): §7.1b, shared GEO helpers in `_geo.py` (GDS4299 verified bit-identical after the refactor); BE(2)c/SK-N-BE(2)c found to be subclone/parent by ATCC number, disclosed via `cell_line_identity`. §1.4a: every Mann–Whitney result now reports `min_p_achievable` (0.333 at 2 vs 2), shown in the UI. §7.4: legacy symbols get a 404 naming the current symbol; the panel resolves them with a per-row note. BH now skips NaN like `p.adjust`. Verified: 229 backend tests, 39 frontend tests, build clean; 414 served panel values match the raw GEO files; every new test checked to fail with its bug reintroduced; Perturbation Contrasts and Compare checked in the browser; every other pane's endpoint (signature, rank, correlation, co-expression, PCA, aux layers, survival refusal) checked on both datasets against the running backend, not clicked through in the browser. DESeq2 comparison not done (R DESeq2 install failed; `pydeseq2` not installed by permission). |
| 2026-11-15 | Two GEO datasets added for the L-MYC / PP2A question. **GSE324984** (MYCL overexpression in C4-2B, knockdown in NCI-H660; RNA-seq, n=2 per arm) with §11 — pre-specified within-cell-line contrasts computed at ingest by an edgeR + limma-trend re-implementation that matches R (limma 3.54.2 / edgeR 3.40.2) to ≤ 2e-12 on all 78,691 genes; new `contrasts` table (alembic `b3c1f0a2d4e5`), `/contrasts`, `/contrast-panel`, `/contrasts/{id}/top`, and a Perturbation Contrasts pane with the researcher's 23-gene PP2A preset. Source-file defects found and handled: TPM file truncated at Excel's 65,535 rows (13,156 genes missing), one gene mis-named in the counts file (TP63 as "KDM5A"), TPM and counts from different quantifications (12% / 5% of genes change direction between them — flagged per gene). **GSE66586** (neuroblastoma lines, 5 MYCN-amplified / 3 non-amplified / 2 non-malignant; microarray) with §7.1b — GEO helpers shared via `_geo.py` (GDS4299 bit-identical), legacy symbols renamed, and two non-independent line pairs disclosed (SK-N-BE(2)-C vs parental SK-N-BE(2) per ATCC numbers despite the sample title; SH-SY5Y ⊂ SK-N-SH). §1.4a added. Every value served for the panel (414) checked against the raw GEO files by an independent parser. |
| 2026-11-15 | Independent review of the above found and this fixed: §1.4a **scipy's vectorised Mann–Whitney fell back to the normal approximation for every gene when any gene had ties** (pre-existing; served p below the exact floor, could cross 0.05 at 3 vs 3) — replaced by an exact permutation test for ≤ 5,000 splits, which changes p-values for small comparisons in every dataset; §7.8 documented tie-break by |median difference| was never applied; §11.3 shared-name loci chosen over the truncated TPM file (CAST resolved to a silent locus), panel and top list disagreeing on loci, shared names stored as aliases (C4orf36 mis-resolved), count-only genes unreachable; mixed-case symbols (C1orf112) never resolvable on Ensembl-indexed datasets including TARGET; §11.2 limma infinite-prior edge case and Infinity in JSON; counts-file label check; contrast sample accounting; cross-cell-line pairwise tests on Compare now flagged (`crosses_stratum`). Each fix has a test validated by reintroducing the bug. Final: backend 259 passed (1 skipped), frontend 43 passed, tsc and production build clean; Perturbation Contrasts pane checked in the browser against the API. |
| 2026-11-15 | Whole-app verification pass after the above, against the running system. API sweep: every endpoint × all 5 datasets (107 checks: strict JSON, sample accounting adds up, every p ≥ its floor) — found and fixed a 500 on DepMap Compare by subtype (mixed exact/approximate pairs stored the None floor as NaN). Existing results re-run through the committed HEAD code side by side on the same data: GDS4299, TARGET, DepMap CRISPR (MYB r = −0.652361, n = 91), differential, signature and correlation are bit-identical; only DepMap's small-subtype pairwise p-values change, as intended by §1.4a. (Two older documented numbers did not match the served ones; both causes confirmed, neither is a data error: §7.1's GDS4299 p-values are Welch t-tests, not the app's Mann–Whitney — recomputed both ways, §7.1 now says so; the 2026-09-13 DTX1-by-NOTCH1 values (medians 18.86 vs 9.69, p = 0.00494) reproduce exactly from `data/cache_backup_pre_fix/`, and the current 18.86 vs 8.71, p = 0.00608 is the result of the §8.6a deterministic file-selection fix, which changed 29 cases' expression.) Top-gene lists: all 2,000 served rows (both contrasts, both directions, n = 500) match R limma to ≤ 1.5e-12 and contain exactly R's top-500 sets; the C4-2B and H660 top-genes CSVs match the API row for row. Browser: every pane exercised on the new datasets. Exports checked by capturing the downloads in-page: the panel CSV (46 rows) matches the API value for value; the top-genes CSV matches row order; the PNG export lost all colour and MYCL's "+8.24 ›" label because inline `style="fill: var(--…)"` overrode the inlined literal paint — fixed in the shared exporter (also affects Compare/Survival PNG strokes), with a test. Min-p notice now merges mirrored sizes ("4 vs 3" = "3 vs 4"). §11.2 DESeq2 concordance added. Backend 260 passed (1 skipped), frontend 45 passed, tsc and build clean. |
| 2026-11-15 | §5.7 export audit — every PNG export in the app exported and inspected in the browser (Compare box plot, Perturbation Contrasts panel, Correlation scatter, PCA, Survival). Found and fixed, each with a test validated by reintroducing the defect: PCA coloured its first and third groups in two near-identical blues (`--accent` #2455A4 vs `--cool` #2C6E9B), so GSE66586's non-amplified and non-malignant lines were indistinguishable on screen and in the export — palette reordered (blue, red, amber, …); PCA and Survival PNGs had no legend (HTML legends are not part of the image), so the image did not say which colour was which group or which curve was HIGH — a text colour key is now written into the image; Correlation, PCA and Survival PNG titles did not name the dataset — now do, with the pooled-experiments note for designed datasets; Compare PNGs printed p = 0.0952 with no small-sample floor — the floor is now in the image; long titles and notes ran off the image edge — wrapped. Frontend 51 passed, tsc and build clean. |
| 2026-11-15 | Remaining usability/identity gaps closed, each verified in the browser and pinned by a test validated by reintroducing the defect. Signature Score listed bare GSM ids; it now reads the same scores from `/rank-signature` (identical to `/signature-score` on all 5 datasets, max difference 0.0, both methods) and shows each sample's cell line and groups. A renamed gene (C1orf112, now FIRRM) returned a bare "not found" on GSE324984: the gene-annotation panel, which already knew the current symbol from MyGene.info, now says the searched name is a previous one and offers "Search FIRRM instead" — an explicit choice, never a silent re-query (§7.4). Compare's skipped-dataset line printed a double full stop. Differential and Sample Ranking CSVs did not name their dataset (same group pair from two datasets overwrote one file): dataset id now in filename and every row; differential CSV carries the exact-test floor. Differential, Sample Ranking, panel and top-genes CSVs and every PNG export checked against the API (corrected: the co-expression CSV was not yet checked at this point — see next line). Backend 260 passed (1 skipped), frontend 56 passed, tsc and build clean, 107-check API sweep clean. |
| 2026-11-15 | Co-expression network CSV — the one CSV export not yet checked — exported in the browser: matched the API on all 20 rows but carried no dataset in its filename or rows (one gene's network from two datasets overwrote one file); fixed like the others, with a test validated by reverting the fix. The app has exactly five CSV exports (Compare differential, Perturbation Contrasts panel and top genes, Sample Ranking, co-expression), all now verified value-for-value against the API; Correlation scatter, PCA and Survival export PNG only, all verified. §7.1 and §8.2 now lead with the values the running app serves (Mann–Whitney on GDS4299; post-§8.6a DTX1/MYCN), the earlier figures kept as dated history. Frontend 57 passed, backend 260 passed. |
| 2026-11-15 | §9.1a: backend tests are now network-free by construction (`tests/conftest.py` refuses non-loopback connections). Found while waiting on PR #1's CI: one listing test was silently downloading every real cohort, ~19 of CI's ~20 minutes. Suite on an empty cache: 264 passed, 1 skipped in ~7 s. |
| 2026-11-15 | §11.5: the panel's "Across both experiments" column said "Opposite direction, not significant in both" for CIP2A, SET and PPP2R1B, which are q < 0.002 in C4-2B — readable as "significant in neither" — and gave "Same direction in both" to PPP2R2C (q < 0.05 in both) and ARPP19 (in neither) alike. Now names which experiment reaches q < 0.05, and says "Tested only in NCI-H660" for MYCN/PPP2R2B instead of a dash. Found in user-standpoint testing; same session re-checked the 23-gene panel against an independent R edgeR/limma run on the raw GEO counts (44 gene × contrast pairs, max |ΔlogFC| 4.5e-15, max relative Δp and Δq 1.5e-12) and all 184 replicate TPMs against the raw TPM file (0 mismatches). Test validated by restoring the old wording. |
| 2026-11-15 | §5.1/5.2/5.7, from user-standpoint testing of GSE66586 in Expression Compare: box plots of log-scale units (log2 intensity, log2 TPM) no longer start the axis at 0 or offer log-of-log; ticks are round and of one precision; the axis names the unit; groups 1 and 3 no longer share near-identical blues; small-n sample dots are visible; box-plot PNGs print FDR q beside p; export date is local. `/compare-by-mutation` now returns `expression_unit` so the Multi-omics box plot gets the same axis. Statistics unchanged: PPP2R5B by MYCN status re-derived by hand (exact Mann–Whitney 2/56, 2/21, 2/10; BH 0.107/0.143/0.2) and matches the served table. Each fix has a test validated by reintroducing the defect; frontend 67 passed, backend 264 passed. |
| 2026-11-15 | §5.7 correlation/PCA scatter, from user-standpoint testing: correlation scatter coloured by group with a pooled-groups note (GSE66586 MYCN vs PPP2R5B r = −0.847 re-derived from the raw series matrix); PCA on TARGET coloured nothing by `etp_status` (top-level column ignored) — both endpoints now use `_sample_groups`, verified against `/group-values` (146 non-ETP / 25 near-ETP / 19 ETP); ungrouped points grey instead of the first group's blue; units and round ticks on scatter axes. Tests validated by reintroducing each defect; backend 265 passed, frontend 72 passed. |
| 2026-11-15 | §2.4a: gene-set tag input split nothing — a typed or pasted list became one unresolvable tag, and unconfirmed text was dropped on clicking Score. Now parsed like the gene panel; tests validated by reintroducing each defect. Found scoring the PP2A inhibitors on GSE66586 as a user. |
| 2026-11-15 | §2.4: "log2 mean" signature scores re-logged data already on a log2 scale (DepMap, GDS4299, GSE66586), giving mean log2(log2(x)+1). Fixed for both endpoints; GSE66586 scores now match the raw GEO probe means exactly. AUC scoring and survival unaffected. Tests validated by reintroducing the defect; backend 270 passed. |
| 2026-11-15 | §3 Survival, from user-standpoint testing: TARGET MEF2C median split re-run in R `survival` on the same scores — n 466 (233/233), log-rank p 1.797e-13 vs app 1.798e-13, Cox HR 7.003 (3.904–12.562), p 6.628e-11, identical. Day-axis ticks now round (0, 1000 … 4000, was 1096, 2192 …); a caption under the Cox table says the signature_score hazard ratio is per 1.0 of a 0–1 AUC score — its whole range — since "HR 27.7" for ETP-TF5 otherwise reads as a per-unit effect. Censoring marks are not drawn (the API does not return censoring times); noted, not done. |
