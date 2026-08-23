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
    coefficients: Record<string, { coef: number; "exp(coef)": number; p: number }>;
    log_likelihood_p_value: number;
  } | null;
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

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await fetch(`${BASE}${path}`, {
    ...init,
    headers: { "Content-Type": "application/json", ...init?.headers },
  });
  if (!res.ok) {
    const body = await res.text();
    throw new Error(`${res.status} ${res.statusText}: ${body}`);
  }
  return res.json();
}

export const api = {
  listDatasets: () => request<{ datasets: DatasetSummary[] }>("/datasets"),
  datasetInfo: (datasetId: string) => request<DatasetSource>(`/datasets/${datasetId}`),
  compare: (datasetId: string, gene: string, groupColumn: string) =>
    request<CompareResult>(
      `/datasets/${datasetId}/compare?gene=${encodeURIComponent(gene)}&group_column=${encodeURIComponent(groupColumn)}`,
    ),
  signatureScore: (datasetId: string, genes: string[], method: "auc" | "log2_mean" = "auc") =>
    request<{ genes: string[]; method: string; scores: Record<string, number> }>(
      `/datasets/${datasetId}/signature-score`,
      { method: "POST", body: JSON.stringify({ genes, method }) },
    ),
  survival: (datasetId: string, genes: string[], covariates: string[] = []) =>
    request<SurvivalResult>(`/datasets/${datasetId}/survival`, {
      method: "POST",
      body: JSON.stringify({ genes, covariates }),
    }),
  rankByGene: (datasetId: string, gene: string) =>
    request<RankResult & { gene: string }>(`/datasets/${datasetId}/rank?gene=${encodeURIComponent(gene)}`),
  rankBySignature: (datasetId: string, genes: string[], method: "auc" | "log2_mean" = "auc") =>
    request<RankResult & { genes: string[]; method: string }>(`/datasets/${datasetId}/rank-signature`, {
      method: "POST",
      body: JSON.stringify({ genes, method }),
    }),
};
