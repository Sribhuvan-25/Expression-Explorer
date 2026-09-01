// In dev, Vite's proxy (vite.config.ts) forwards /api to localhost:8420.
// In a real deployment, set VITE_API_BASE_URL at build time (see
// frontend/Dockerfile) to wherever the backend actually runs — this is
// never hardcoded to a specific host.
const BASE = import.meta.env.VITE_API_BASE_URL ?? "/api";

export interface DatasetSummary {
  dataset_id: string;
  display_name: string;
  group_columns: string[];
  supports_survival: boolean;
  // Provenance, included so the sidebar can show what a dataset actually
  // is without a follow-up request per dataset. Optional: a dataset that
  // failed to load is still listed, just without these.
  n_samples?: number;
  assay_type?: string;
  accession?: string;
  // Provenance caveats -- e.g. DepMap's release-staleness disclosure
  // (METHODS.md 7.5). Was computed correctly by the backend but had no
  // reader anywhere in the frontend until this field was added: the
  // per-dataset /datasets/{id} endpoint carried it, but nothing in the UI
  // ever called that endpoint, so a real, correctly-computed warning was
  // silently invisible to every user.
  notes?: string;
}

export interface GeneInfo {
  symbol: string;
  name: string | null;
  summary: string | null;
  aliases: string[];
  ensembl_gene_id: string | null;
  entrez_gene_id: string | number | null;
}

export interface DatasetSource {
  dataset_id: string;
  display_name: string;
  accession: string;
  repository: string;
  assay_type: string;
  expression_unit: string;
  gene_identifier: string;
  n_samples: number;
  disease_area: string;
  notes?: string;
  group_columns: string[];
  supports_survival: boolean;
}

export interface ComparePoint {
  sample_id: string;
  group: string;
  value: number;
}

export interface PairwiseTest {
  group_a: string;
  group_b: string;
  n_a: number;
  n_b: number;
  u_stat: number;
  p_value: number;
  q_value?: number;
}

export interface CompareResult {
  gene: string;
  group_column: string;
  points: ComparePoint[];
  pairwise_tests: PairwiseTest[];
  kruskal_wallis: { h_stat: number | null; p_value: number | null };
  // How much of the cohort this comparison actually ran on. A grouping
  // column populated for only part of a dataset (TARGET classifies
  // etp_status for 190 of 469) otherwise silently narrows the analysis.
  n_dataset_total: number;
  n_excluded: number;
  // null when n_excluded is 0 -- no reason to give when nothing was
  // excluded (previously always a string, even claiming an exclusion
  // reason on a 100%-covered grouping column).
  exclusion_reason: string | null;
}

export interface CompareMultiDatasetResult {
  dataset_id: string;
  display_name: string;
  skipped: boolean;
  skip_reason?: string;
  assay_type?: string;
  expression_unit?: string;
  n_dataset_total?: number;
  n_excluded?: number;
  exclusion_reason?: string | null;
  points?: ComparePoint[];
  pairwise_tests?: PairwiseTest[];
  kruskal_wallis?: { h_stat: number | null; p_value: number | null };
}

export interface CompareMultiResult {
  gene: string;
  group_column: string;
  datasets: CompareMultiDatasetResult[];
}

export interface SurvivalPoint {
  days: number;
  survival_probability: number;
}

export interface SurvivalCurve {
  n: number;
  points: SurvivalPoint[];
}

export interface SurvivalResult {
  genes: string[];
  n: number;
  curves: Record<string, SurvivalCurve>;
  logrank_p_value?: number;
  cox: {
    n: number;
    coefficients: Record<
      string,
      {
        coef: number;
        "exp(coef)": number;
        "exp(coef) lower 95%": number;
        "exp(coef) upper 95%": number;
        p: number;
      }
    >;
    log_likelihood_p_value: number;
  } | null;
  // Patients with no follow-up time at all can't be placed on a survival
  // curve and are dropped -- surfaced so "n = 466" on a 469-sample dataset
  // doesn't read as a different cohort.
  n_dataset_total: number;
  n_excluded: number;
  exclusion_reason: string | null;
}

export interface CorrelationPoint {
  sample_id: string;
  x: number;
  y: number;
}

export interface CorrelationResult {
  gene_a: string;
  gene_b: string;
  method: "pearson" | "spearman" | "kendall";
  n: number;
  coefficient: number;
  p_value: number;
  points: CorrelationPoint[];
}

export interface CoExpressionEntry {
  gene: string;
  symbol: string;
  coefficient: number;
}

export interface CoExpressionResult {
  gene: string;
  method: "pearson" | "spearman";
  direction: "positive" | "negative";
  network: CoExpressionEntry[];
}

export interface PCAPoint {
  sample_id: string;
  pc1: number;
  pc2?: number;
  [key: `pc${number}`]: number | undefined;
  group_columns?: Record<string, string | null>;
}

export interface PCAResult {
  n_components: number;
  n_samples: number;
  genes_used: string[];
  genes_missing: string[];
  genes_zero_variance: string[];
  genes_unrecognized: string[];
  variance_explained: number[];
  points: PCAPoint[];
}

export interface DifferentialGeneRow {
  feature_id: string;
  symbol: string;
  p_value: number;
  q_value: number;
  median_diff: number;
  median_a: number;
  median_b: number;
}

export interface DifferentialResult {
  group_column: string;
  group_a: string;
  group_b: string;
  n_a: number;
  n_b: number;
  n_genes_tested: number;
  top_n: number;
  genes: DifferentialGeneRow[];
}

export interface RankRow {
  sample_id: string;
  value: number;
  [metadataKey: string]: string | number | null;
}

export interface RankResult {
  n: number;
  rows: RankRow[];
}

// Carries the HTTP status as a real property instead of encoding it into
// the message text -- callers that need to distinguish "this request was
// invalid, don't retry it" (4xx) from "this was transient, worth another
// try" (5xx/network) read error.status directly rather than parsing the
// leading characters of a message meant for display.
export class ApiError extends Error {
  status: number;
  constructor(status: number, message: string) {
    super(message);
    this.name = "ApiError";
    this.status = status;
  }
}

// FastAPI error responses come in two shapes: our own HTTPException(status,
// "message") calls always produce {detail: string}; a request that fails
// Pydantic validation before a handler ever runs (a missing/malformed
// field) produces {detail: [{loc, msg, type, ...}, ...]}. Both used to be
// shown to the user as the raw response body -- '404 Not Found:
// {"detail":"Gene \'X\' not found in dataset."}' -- which is legible with
// effort but is JSON plumbing, not a sentence a domain expert should have
// to parse. This extracts the human-readable message from either shape,
// falling back to the raw body only if the response isn't the JSON error
// shape the backend actually sends (e.g. a proxy timeout page, or the
// backend being down entirely).
async function parseErrorMessage(res: Response): Promise<string> {
  const body = await res.text();
  try {
    const parsed = JSON.parse(body);
    const detail = parsed?.detail;
    if (typeof detail === "string") return detail;
    if (Array.isArray(detail)) {
      return detail
        .map((e) => {
          const field = Array.isArray(e?.loc) ? e.loc.filter((p: unknown) => p !== "body").join(".") : null;
          const msg = typeof e?.msg === "string" ? e.msg : "Invalid value";
          return field ? `${field}: ${msg}` : msg;
        })
        .join("; ");
    }
  } catch {
    // Not JSON -- fall through to the generic message below.
  }
  return body || `${res.status} ${res.statusText}`;
}

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${BASE}${path}`, {
    ...init,
    headers: { "Content-Type": "application/json", ...init?.headers },
  });
  if (!res.ok) {
    throw new ApiError(res.status, await parseErrorMessage(res));
  }
  return res.json();
}

export const api = {
  listDatasets: () => request<{ datasets: DatasetSummary[] }>("/datasets"),
  datasetInfo: (datasetId: string) => request<DatasetSource>(`/datasets/${datasetId}`),
  geneInfo: (symbol: string) => request<GeneInfo>(`/genes/${encodeURIComponent(symbol)}`),
  compare: (datasetId: string, gene: string, groupColumn: string) =>
    request<CompareResult>(
      `/datasets/${datasetId}/compare?gene=${encodeURIComponent(gene)}&group_column=${encodeURIComponent(groupColumn)}`,
    ),
  compareMulti: (gene: string, groupColumn: string, datasetIds: string[]) =>
    request<CompareMultiResult>(
      `/datasets/compare-multi?gene=${encodeURIComponent(gene)}&group_column=${encodeURIComponent(groupColumn)}&dataset_ids=${encodeURIComponent(datasetIds.join(","))}`,
    ),
  signatureScore: (datasetId: string, genes: string[], method: "auc" | "log2_mean" = "auc") =>
    request<{ genes: string[]; method: string; scores: Record<string, number> }>(
      `/datasets/${datasetId}/signature-score`,
      { method: "POST", body: JSON.stringify({ genes, method }) },
    ),
  survival: (
    datasetId: string,
    genes: string[],
    covariates: string[] = [],
    cutoff?: { method: "median" | "quartile" | "custom"; highPct?: number; lowPct?: number },
  ) =>
    request<SurvivalResult>(`/datasets/${datasetId}/survival`, {
      method: "POST",
      body: JSON.stringify({
        genes,
        covariates,
        cutoff_method: cutoff?.method ?? "median",
        cutoff_high_pct: cutoff?.highPct ?? 50,
        cutoff_low_pct: cutoff?.lowPct ?? 50,
      }),
    }),
  rankByGene: (datasetId: string, gene: string) =>
    request<RankResult & { gene: string }>(`/datasets/${datasetId}/rank?gene=${encodeURIComponent(gene)}`),
  rankBySignature: (datasetId: string, genes: string[], method: "auc" | "log2_mean" = "auc") =>
    request<RankResult & { genes: string[]; method: string }>(`/datasets/${datasetId}/rank-signature`, {
      method: "POST",
      body: JSON.stringify({ genes, method }),
    }),
  correlation: (datasetId: string, geneA: string, geneB: string, method: "pearson" | "spearman" | "kendall" = "pearson") =>
    request<CorrelationResult>(
      `/datasets/${datasetId}/correlation?gene_a=${encodeURIComponent(geneA)}&gene_b=${encodeURIComponent(geneB)}&method=${method}`,
    ),
  coExpression: (
    datasetId: string,
    gene: string,
    options?: { topN?: number; method?: "pearson" | "spearman"; direction?: "positive" | "negative" },
  ) =>
    request<CoExpressionResult>(
      `/datasets/${datasetId}/co-expression?gene=${encodeURIComponent(gene)}&top_n=${options?.topN ?? 20}&method=${options?.method ?? "pearson"}&direction=${options?.direction ?? "positive"}`,
    ),
  pca: (datasetId: string, genes: string[], nComponents: number = 2) =>
    request<PCAResult>(`/datasets/${datasetId}/pca`, {
      method: "POST",
      body: JSON.stringify({ genes, n_components: nComponents }),
    }),
  differential: (datasetId: string, groupColumn: string, groupA: string, groupB: string, topN: number = 50) =>
    request<DifferentialResult>(
      `/datasets/${datasetId}/differential?group_column=${encodeURIComponent(groupColumn)}&group_a=${encodeURIComponent(groupA)}&group_b=${encodeURIComponent(groupB)}&top_n=${topN}`,
    ),
  groupValues: (datasetId: string, groupColumn: string) =>
    request<{ group_column: string; values: { value: string; n: number }[] }>(
      `/datasets/${datasetId}/group-values?group_column=${encodeURIComponent(groupColumn)}`,
    ),
};
