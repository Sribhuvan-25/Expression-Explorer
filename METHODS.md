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
  only reporting one of the two causes.
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
**Decided and implemented (2026-08-27).** Datasets are shown
**side by side, each with its own axis and its own statistics.** No
pooling of values across datasets, ever, by default.
`GET /datasets/compare-multi` (`backend/app/api/main.py`) runs the
existing single-dataset `expression_by_group` independently per dataset
and returns one block per dataset — points, pairwise tests, Kruskal-Wallis,
assay type, expression unit, sample-exclusion accounting — with no shared
axis or merged statistic anywhere in the response.

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
