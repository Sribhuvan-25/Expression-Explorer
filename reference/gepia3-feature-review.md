# GEPIA3 feature review — gap analysis against Expression Explorer

Reviewed live at https://gepia3.bioinfoliu.com on 2026-08-27, by walking every
menu item in the running app (not from documentation). GEPIA3 is a public
pan-cancer expression portal built on TCGA (33 cancer types, tumor + matched
normal) and GTEx (54 normal tissue types), plus several bolt-on data sources
for drugs, networks, and RNA-level alterations.

This is not a build plan — it's the raw comparison so we can decide together
what's actually worth adding versus what's out of scope for this project's
purpose (ETP-ALL / T-ALL research, not a general pan-cancer portal).

## What GEPIA3 covers, end to end

### Data backbone
- **TCGA**: 33 cancer types, tumor samples + (where available) paired
  peritumor/normal — e.g. BRCA 1108 tumor / 133 normal, LIHC 379/59.
- **GTEx**: 54 normal tissue types (not cancer-matched), used as an
  alternative/additional normal baseline where TCGA's own normal count is
  thin or absent (e.g. ACC has 92 tumor, 0 TCGA-normal, but 161 GTEx Adrenal
  Gland samples fill that gap).
- All expression on a single processing pipeline, log2(TPM+1), so every
  cancer type is comparable on one scale within TCGA/GTEx.

### Single-gene features (Gene Card / Profile Summary)
- **Gene Card**: one gene's annotation (Ensembl ID, aliases, RefSeq summary,
  links out to NCBI/OMIM/COSMIC/HPA/DrugBank/Xena/cBioPortal) plus an
  **interactive anatomical bodymap** — a clickable body diagram colored by
  median tumor vs. normal expression per organ/tissue.
- **Differential Genes**: genome-wide (not per-query-gene) ranked
  differential-expression list for one chosen cancer type, filterable by
  biotype (protein coding / lncRNA / TCR / IG genes) and gene-vs-isoform
  level, limma or DESeq2.
- **Prognostic Genes**: genome-wide scan for genes significantly associated
  with survival in one chosen cancer type (the inverse of a single-gene
  survival lookup).
- **Find Drugs**: real-world TCGA clinical drug-treatment records — did
  receiving a named drug (from a curated list of ~75 real chemo/targeted
  agents) correlate with survival, per cancer type — with a patient-count
  heatmap (cancer type × drug).

### Expression Analysis
- **Expression DIY** (4 modes): Profile DIY (differential expression across
  multiple selected cancer types at once, log2FC + q-value cutoffs,
  limma/DESeq2, choice of paired-peritumor / +GTEx / tumor-only baseline),
  Boxplot (with high/low group split by a second gene), Stage Plot
  (expression vs. clinical tumor stage — major stage / sub-stage),
  Multiple Gene Comparison (a gene set plotted together, not one gene).
- **Correlation Analysis**: gene-vs-gene (or curated immune-cell geneset
  vs. gene) correlation — Pearson/Spearman/Kendall — across any mix of TCGA
  tumor, TCGA peritumor, and GTEx normal as the sample pool.
- **Dimensionality Reduction**: PCA/embedding view of a gene set across
  selected sample pools.
- **Hotspot Mutation**: gene expression compared between mutated vs.
  wildtype samples (gene-level or amino-acid-level mutation calls,
  Wilcoxon/ANOVA) — requires MAF-level mutation calling data.

### Survival Analysis
- **Univariable**: single gene/geneset, Overall Survival *or*
  Progression-Free Interval as the endpoint, Median/Quartile/**Custom**
  cutoff (not just median), HR + 95% CI, across multiple cancer types.
- **Multivariable**: several genes as separate Cox covariates at once
  (not one combined signature score), across multiple cancer types.

### Drug Analysis
- **TCGA Drug Response**: does gene expression differ between RECIST
  response categories (Complete/Partial Response vs. Stable/Progressive
  Disease) for patients on a given drug?
- **TCGA Drug vs Survival**: same real-world drug-treatment survival
  question as "Find Drugs" above.
- **Cell Line Screen**: correlate a gene's expression in cancer cell lines
  against drug sensitivity (IC50/AUC/IC90/EC50/Einf) from GDSC1/GDSC2/
  CREAMMIST/CTRP screens, filterable by tissue/cancer type of the cell line.
- **Cell CRISPR Screen**: compare CRISPR gene-knockout fitness effects
  under drug exposure vs. untreated, from 4 published high-throughput
  screens.

### Network Analysis
- **Expression Network**: top-N co-expressed genes for a query gene
  (positive or anti-correlated), across any TCGA/GTEx sample pool.
- **Alteration Network, Protein Interaction, SL/SV (synthetic
  lethality/viability), eQTL, Comprehensive Analysis** — all require data
  we don't have any pipeline for (genotype/SNP data, curated PPI databases,
  CRISPR co-dependency screens).

### RNA Alterations
- **Allele-Specific Expression, Alternative Promoter, Gene Fusion** — all
  require raw-read-level or structural-variant-level processing (not just
  a quantified expression matrix), sourced partly from ICGC rather than
  TCGA/GEO.

### Linked companion tool: GEPIA2021 (separate domain)
- **Deconvolution Analysis**: infers immune/stromal cell-type proportions
  from bulk RNA-seq (CIBERSORT/EPIC/quanTIseq), then supports four
  downstream views: Proportion (boxplot + ANOVA across cell types or
  cohorts), Correlation (two cell types' proportions against each other),
  Sub-expression (a gene's expression *within* one inferred cell type,
  with cell-type-level differential expression), Survival (split by
  inferred cell-type proportion, Kaplan-Meier + log-rank).

---

## Gap analysis: what we have vs. what's missing

| Capability | We have it? | Notes |
|---|---|---|
| Multi-cancer-type / multi-dataset selection for one query | **Yes, just added** | Compare page now defaults to every dataset exposing the chosen grouping column, with checkboxes to narrow — same pattern GEPIA3 uses for cancer-type multi-select. |
| Pan-cancer breadth (33 TCGA cancer types) | No, by design | We're scoped to T-ALL/ETP-ALL (TARGET-ALL-P2, GDS4299) + DepMap lymphoid lines. Adding pan-cancer breadth is a scope decision, not a bug — flagging so it's explicit. |
| Gene annotation + external DB links (Gene Card) | No | Cheap, high-value: gene symbol → Ensembl ID, aliases, RefSeq summary, links to NCBI/GeneCards/COSMIC/DrugBank. All derivable from a public gene-info source (NCBI `gene_info`/MyGene.info), no new sample data needed. Also solves METHODS.md §7.4 (gene-symbol aliases), which is already an open item. |
| Genome-wide differential-expression browse (not just per-query-gene) | No | "Show me the top DE genes for this dataset/grouping" — inverse of current gene-first UX. Computable today from data we already have (no new ingestion), just a different query shape. |
| Genome-wide prognostic-gene scan | No | Same idea for survival — "which genes matter" instead of "does this one gene matter." Only applicable to `target_all_p2` today (only survival-capable dataset). Computationally expensive (one Cox fit per gene) — needs to be precomputed/cached, not live-queried. |
| Custom survival cutoff (vs. median-only) | No | We binarize by median only (`backend/app/analysis/survival.py`). Adding quartile/custom-percentile is a small, contained change. |
| Progression-Free Interval as a second survival endpoint | No | We only have Overall Survival (`vital_status`/`days_to_death`/`days_to_last_follow_up`). Would need PFI-relevant fields in the underlying clinical data — GDC's TARGET clinical files should be checked for whether PFI fields exist before assuming this is free. |
| Multi-gene Cox (each gene its own covariate) | Partial | We fit Cox on one signature score plus optional extra covariates. GEPIA3 also supports several genes as independent covariates without collapsing to a score — a different modeling choice, not strictly better, worth a METHODS.md discussion before adding. |
| Gene-vs-gene correlation analysis | No | Real gap. Pearson/Spearman correlation between two genes (or a gene and a curated immune-cell signature) across a dataset. Fully computable from our existing matrices — no new data needed. Reasonably high value for T-ALL biology (e.g. does MYCN correlate with a specific TF). |
| Curated immune-cell gene signatures (Naive/Effector/Treg T-cell etc.) | No | GEPIA3 ships these as example genesets for correlation/DE/dimensionality reduction. We'd need to either source a published signature list (e.g. from a cited immunology paper) or skip this — not something to invent ourselves. |
| Dimensionality reduction (PCA) across samples | No | Computable today from any of our matrices. Useful for "do ETP and non-ETP actually separate on more than one gene" sanity checks. |
| Tumor-stage-stratified expression (Stage Plot) | No | Requires clinical stage fields. Would need to check whether GDC's TARGET clinical data / GDS4299 supplementary tables carry a stage-equivalent field (pediatric ALL typically doesn't stage the way solid tumors do — may not apply to our cohorts at all). |
| Mutation-status-stratified expression (Hotspot Mutation) | No | Requires MAF-level somatic mutation calls per sample — a new data type entirely, not derivable from expression matrices. Would need a new ingest pipeline (GDC does have MAF files for TARGET-ALL-P2's exome data). |
| Drug-treatment survival / drug response (real-world clinical) | No | Requires per-patient treatment/drug records — GDC clinical supplements may have this for TARGET-ALL-P2 (worth checking), but it's a new ingest surface, not a derived analysis. |
| Cell-line drug sensitivity (GDSC/CTRP-style) | No | Directly relevant to our existing `depmap` dataset — DepMap's own portal also publishes companion drug-sensitivity screens (PRISM repurposing, GDSC) keyed by the same `ModelID`s we already ingest. This is the single most natural extension: same cell lines, same identifiers, just a second file to pull from the same Figshare release or a linked one. |
| Co-expression network (top-N correlated genes) | No | Computable from existing matrices (correlate the query gene's row against every other row, rank). No new data needed. |
| Immune-cell deconvolution (CIBERSORT/EPIC/quanTIseq) | No | Meaningful for T-ALL (leukemia = abnormal hematopoietic/immune cell state), but this is a modeling addition (a deconvolution algorithm + a reference signature matrix), not just a query — heavier lift than most of the above. |
| Gene fusion / allele-specific expression / eQTL / SL-SV / PPI network | No | All require data types we have no pipeline for at all (structural variants, genotypes, curated interaction databases, CRISPR co-dependency screens). Flagging for completeness, not recommending — none of these fit a project scoped to expression + survival in T-ALL. |

---

## Recommendation: what's actually worth adding here

Grouped by effort vs. value, given this project is scoped to T-ALL/ETP-ALL
research rather than a general pan-cancer portal:

**Low effort, clearly in scope (computable from data we already have):**
1. Gene-vs-gene correlation analysis (Pearson/Spearman) within a dataset —
   also incidentally useful for validating grouping choices.
2. Co-expression network (top-N correlated genes for a query gene).
3. Dimensionality reduction (PCA) across a dataset's samples.
4. Custom/quartile survival cutoff, not just median.
5. Gene annotation panel (Ensembl ID, aliases, summary, external DB links)
   — also directly closes METHODS.md §7.4 (gene-symbol alias resolution).
6. Genome-wide "top differentially expressed genes" browse per
   dataset/grouping (inverse of current single-gene query).

**Medium effort, needs a new (but accessible) data source:**
7. DepMap cell-line drug sensitivity (PRISM/GDSC) — same `ModelID`s we
   already have, likely the same or an adjacent Figshare release.
8. Genome-wide prognostic-gene scan on `target_all_p2` — needs
   precomputation/caching since it's ~20k Cox fits, not a live query.

**Needs a decision or new data before it's buildable — flagging, not
starting:**
9. Immune-cell deconvolution — real value for leukemia biology, but a
   genuinely new analytical component (needs a reference signature matrix
   and an algorithm choice), not a small addition.
10. Mutation-status-stratified expression — GDC does have MAF files for
    TARGET-ALL-P2; would need a new ingest module.
11. Real-world drug-treatment survival — would need to confirm GDC's
    TARGET clinical supplements actually carry per-patient drug records
    before treating this as available.

**Out of scope, not recommended:**
- Pan-cancer breadth (33 TCGA types) — a different product, not a gap.
- Gene fusion, eQTL, SL/SV, protein-interaction network, allele-specific
  expression, tumor-stage plots, CRISPR-drug screens — each needs a data
  type this project has no pipeline for and no clear T-ALL research
  motivation cited yet.
