# Implementation plan: TARGET mutation data + DepMap CRISPR/drug-sensitivity

Research done against live sources (GDC API, Figshare API, real file
downloads), not documentation alone. This is the plan only — no
implementation yet, per instruction to research and plan before touching
code. Cheap wins (correlation, co-expression network, PCA, custom
survival cutoff, gene annotation, top-DE-genes browse) are being built and
deployed first; these two pipelines come after.

---

## Pipeline 1: TARGET-ALL-P2 mutation status (Hotspot Mutation equivalent)

### What was verified
- GDC API confirms **739 open-access MAF files** for TARGET-ALL-P2
  (`data_type=Masked Somatic Mutation`, `data_format=MAF`), one file per
  aliquot — same per-file-per-sample shape as the RNA-seq loader already
  handles, not one combined project MAF.
- Real file downloaded and inspected directly (not assumed from docs):
  standard MAF columns — `Hugo_Symbol`, `Entrez_Gene_Id`,
  `Variant_Classification`, `Variant_Type`, `Chromosome`,
  `Start_Position`, `Reference_Allele`, `Tumor_Seq_Allele1/2`,
  `Tumor_Sample_Barcode` (format `TARGET-10-PAPZPJ-09A-01D` — a GDC
  barcode, needs prefix-matching down to `TARGET-10-PAPZPJ` to join our
  `sample_id`, not an exact match).
- **717 distinct cases have MAF files; our existing `target_all_p2`
  RNA-seq matrix has 469 samples. Overlap is 279 (60% of our RNA-seq
  cohort).** 190 of our RNA-seq samples have no matching mutation data at
  all, and 438 MAF cases aren't in our RNA-seq matrix. This was computed
  directly against the cached matrix, not estimated.
- One sampled file had only 21 called mutations — low per-sample mutation
  burden is expected and well-documented for pediatric ALL, not a parsing
  bug to investigate later.

### What this means for the design
A "gene expression, split by mutation status" feature run on this data
would silently operate on 60% of the cohort unless built with the same
sample-exclusion transparency already used elsewhere in this app
(`n_dataset_total`/`n_excluded`/`exclusion_reason` — see METHODS.md 6.1).
This is not optional here; skipping it would repeat exactly the mistake
the domain-expert review already caught once (the n=466 TARGET survival
confusion).

### Proposed shape
- New ingest module `backend/app/ingest/target_mutations.py`:
  - Lists the 739 open MAF files via the GDC `/files` endpoint (same
    `_client()` pattern as `gdc_target.py`), downloads each
    (`/data/{file_id}`, gzip), extracts `Hugo_Symbol` +
    `Variant_Classification` + `Tumor_Sample_Barcode` (trimmed to the
    `TARGET-##-XXXXXX` submitter_id prefix).
  - Builds a **gene-level boolean mutation matrix**: rows=sample_id,
    columns=gene symbol, True if that sample has >=1 non-silent variant
    in that gene. (Amino-acid-level, GEPIA3's other mode, needs
    `Start_Position`/`HGVSp_Short` retained per-call — worth deferring to
    a v2, since gene-level is the immediately useful case and the same
    download step captures both; no reason to throw away the extra
    columns even if v1 only uses gene-level.)
  - Cache as parquet, same pattern as `gdc_target.py`'s per-sample
    checkpoint files, for the same reason (739 files, avoid holding all
    in memory, survive a partial crash).
- New analysis function `backend/app/analysis/mutation.py`:
  `expression_by_mutation_status(matrix, mutation_matrix, samples,
  gene_id, mutation_gene)` — same shape as `expression_by_group`, reuses
  `pairwise_tests`/`kruskal_wallis` from `compare.py` rather than
  duplicating the stats.
- New endpoint, e.g. `GET /datasets/{id}/compare-by-mutation?gene=X&mutation_gene=Y`
  — deliberately a *separate* endpoint from `/compare`, not a mode flag
  on it, because the grouping source (a metadata column vs. a derived
  mutation call) is different enough in provenance and failure mode
  (silent 40% cohort loss) to want its own explicit exclusion accounting
  rather than being folded into the existing one.
- Applies **only to `target_all_p2`** — `gds4299` and `depmap` have no
  mutation data source available at all; the dataset registry's existing
  `group_columns`/`supports_survival`-style capability flag pattern
  extends naturally to a new `supports_mutation_status: bool` on
  `DatasetDescriptor`.

### Open questions to settle before coding (not blocking research, just decisions)
1. Gene-level only for v1, or build amino-acid-level from day one since
   the same download gets both? (Leaning gene-level v1, given effort
   should go to more datasets/features rather than a level of resolution
   nothing has asked for yet.)
2. Which non-silent `Variant_Classification` values count as "mutated"
   for the boolean matrix (standard MAF practice excludes Silent, Intron,
   5'/3'-UTR, etc. — needs an explicit, documented allowlist, not an
   implicit one).

---

## Pipeline 2: DepMap CRISPR gene-effect + drug sensitivity

### What was verified
- **CRISPRGeneEffect.csv** ships inside the *same* `DepMap 24Q4 Public`
  Figshare bundle (`figshare.com/articles/dataset/DepMap_24Q4_Public/27993966`)
  our `depmap.py` already downloads from — rows=`ModelID` (confirmed same
  join key as `Model.csv`/expression matrix), columns=`SYMBOL (ENTREZID)`
  (confirmed same header format as
  `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, which the existing
  loader already parses with the same transpose-then-column-strip
  pattern). This is the smallest-effort addition of the two — no new
  release, no new join logic, same file-list dict this loader already
  builds.
- **Repurposing_Public_24Q2** (PRISM drug screen) is a **separate**
  Figshare article (`figshare.com/articles/dataset/Repurposing_Public_24Q2/25917643`,
  confirmed public, confirmed real via direct download), not bundled with
  "DepMap Public" — the existing `_RELEASE_TITLE` regex
  (`^DepMap \d\dQ\d Public$`) will never match it; needs its own
  `_resolve_release_urls()`-style lookup with a different title pattern
  (`^Repurposing Public \d\dQ\d$`, confirmed against the real title).
  - `Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv` confirmed to carry
    a `depmap_id` column in `ACH-XXXXXX` format — same join key.
  - The actual sensitivity matrix
    (`Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`) is
    **oriented opposite** to every matrix we currently handle: rows are
    Broad compound IDs (`BRD:BRD-...`), columns are `ModelID`s. Compound
    ID → human drug name requires a join against
    `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv`. This is
    a materially different shape from "gene × sample" — it's
    "drug × cell-line," and a gene only enters by correlating its
    expression (from the existing expression matrix, same `ModelID`
    columns) against a drug's sensitivity row. That correlation is
    exactly the "Correlation Analysis" cheap-win being built now, just
    with the second variable being a drug-sensitivity vector instead of
    a second gene's expression — meaning this pipeline's real payoff
    depends on the correlation feature landing first, which lines up
    with the agreed cheap-wins-first ordering.

### Proposed shape
- Extend `backend/app/ingest/depmap.py`'s existing file dict to also
  pull `CRISPRGeneEffect.csv` from the same already-resolved release
  (near-zero extra ingest complexity — same bundle, same download
  function, same cache directory).
- New module `backend/app/ingest/depmap_repurposing.py` for the PRISM
  data: separate release-resolution (own title regex), separate cache
  dir, joins to the same `Model.csv`/`ModelID`s the main `depmap.py`
  loader already exposes via `samples["sample_id"]`.
- Neither of these becomes a new *registered dataset* in the
  `DatasetDescriptor` sense — they're **auxiliary matrices attached to
  the existing `depmap` dataset**, not new sample cohorts (same 186
  lymphoid cell lines, same `sample_id`s, different measured quantity per
  cell line). The registry's `Dataset` object may need a second optional
  matrix slot (e.g. `Dataset.aux_matrices: dict[str, pd.DataFrame]`) or a
  dedicated `DepMapAuxData` side-object, rather than forcing CRISPR
  effect scores or drug sensitivity through the existing single
  `matrix: pd.DataFrame` gene-expression contract — worth a short design
  pass on `app/models/contract.py` before writing this, not a decision to
  make silently inside the ingest module.
- New endpoint(s), keyed by the same "gene X vs Y" shape as the
  correlation-analysis cheap win: e.g.
  `GET /datasets/depmap/dependency?gene=X` (CRISPR effect score
  distribution/rank across the 186 lines) and
  `GET /datasets/depmap/drug-correlation?gene=X&compound=Y` (expression
  vs. sensitivity correlation, reusing whatever correlation function the
  cheap-win phase builds).

### Open questions to settle before coding
1. Contract shape for auxiliary per-cell-line matrices (extend
   `Dataset`, or a new sibling type) — needs a short design decision, not
   a guess baked into the first PR.
2. Whether to bundle PRISM ingestion in the same pass as CRISPR, or land
   CRISPR first (it's materially cheaper — same bundle, same shape) and
   treat PRISM as its own follow-on given the orientation mismatch and
   extra compound-name join.

---

## Sequencing recommendation

Given the agreed order (cheap wins → deploy → new pipelines), inside
"new pipelines" specifically:
1. **DepMap CRISPR gene-effect** first — smallest marginal effort (same
   release, same bundle, same join key already in use).
2. **TARGET mutation status** second — new ingest surface but a single,
   well-understood file format (MAF) and a clear reuse of existing stats
   functions; the 60%-overlap caveat needs the exclusion-transparency
   pattern applied, which is a known recipe in this codebase already.
3. **DepMap PRISM drug sensitivity** last — genuinely new matrix
   orientation and an extra compound-name join, and its main payoff (gene
   vs. drug correlation) is gated on the correlation-analysis cheap win
   existing first anyway.
