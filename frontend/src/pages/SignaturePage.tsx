import { useEffect, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, EmptyState, ErrorState, DataTable } from "../components/ui";

const PRESETS: Record<string, string[]> = {
  "ETP-TF5 (paper)": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
  "ZMIZ1-5 (paper)": ["MEF2C", "BCL2", "MYB", "MYCN", "ZMIZ1"],
};

export function SignaturePage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);

  const [genes, setGenes] = useState<string[]>(PRESETS["ETP-TF5 (paper)"]);
  const [method, setMethod] = useState<"auc" | "log2_mean">("auc");
  const [activePreset, setActivePreset] = useState<string | null>("ETP-TF5 (paper)");
  const [query, setQuery] = useState<{ datasetId: string; genes: string[]; method: "auc" | "log2_mean" } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["signature", query],
    queryFn: () => api.signatureScore(query!.datasetId, query!.genes, query!.method),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) setQuery({ datasetId, genes, method });
  };

  const rows = data
    ? Object.entries(data.scores)
        .filter(([, v]) => v != null && !Number.isNaN(v))
        .sort((a, b) => b[1] - a[1])
    : [];

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Signature Score"
        title="Score any gene set across samples"
        description="AUCell-equivalent rank-based scoring (or a simpler log2-mean), matching the paper's method. Enter a custom gene set or start from a preset."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <Panel title="Gene set">
          {datasets.length > 1 && (
            <div className="mb-4 max-w-[28ch]">
              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Dataset
              </label>
              <select
                value={datasetId}
                onChange={(e) => setDatasetId(e.target.value)}
                className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent"
              >
                {datasets.map((d) => (
                  <option key={d.dataset_id} value={d.dataset_id}>
                    {d.display_name}
                  </option>
                ))}
              </select>
            </div>
          )}
          <div className="mb-3 flex flex-wrap gap-1.5">
            {Object.keys(PRESETS).map((name) => (
              <PresetButton
                key={name}
                label={name}
                active={activePreset === name}
                onClick={() => {
                  setGenes(PRESETS[name]);
                  setActivePreset(name);
                }}
              />
            ))}
            <PresetButton
              label="Custom"
              active={activePreset === null}
              onClick={() => setActivePreset(null)}
            />
          </div>

          <GeneTagInput
            genes={genes}
            onChange={(g) => {
              setGenes(g);
              setActivePreset(null);
            }}
          />

          <div className="mt-4 flex items-center justify-between">
            <div className="flex items-center gap-3">
              <label className="font-mono text-[10px] uppercase tracking-wider text-ink-mute">Method</label>
              <div className="flex overflow-hidden rounded-[3px] border border-rule">
                {(["auc", "log2_mean"] as const).map((m) => (
                  <button
                    key={m}
                    type="button"
                    onClick={() => setMethod(m)}
                    className={`px-3 py-1 font-mono text-[11.5px] transition-colors ${
                      method === m ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
                    }`}
                  >
                    {m === "auc" ? "AUCell (paper)" : "log2 mean"}
                  </button>
                ))}
              </div>
            </div>
            <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length === 0}>
              Score {genes.length} gene{genes.length === 1 ? "" : "s"}
            </PrimaryButton>
          </div>
        </Panel>

        {error && <ErrorState message={(error as Error).message} />}

        {!query && !error && (
          <EmptyState
            icon={
              <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                <circle cx="12" cy="12" r="8.5" stroke="currentColor" strokeWidth="1.7" />
                <path d="M12 7.5V12L15 14" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
              </svg>
            }
            title="No signature scored yet"
            description="Choose a preset or build a custom gene set, then run to rank every sample by signature score."
          />
        )}

        {data && (
          <Panel
            title={`Scores — ${data.genes.join(", ")}`}
            action={<span className="font-mono text-[10.5px] text-ink-mute">n = {rows.length}</span>}
          >
            <DataTable
              columns={[
                { key: "rank", label: "#", align: "right" },
                { key: "sample", label: "Sample" },
                { key: "score", label: "Score", align: "right" },
              ]}
              rows={rows.slice(0, 25).map(([sample, score], i) => ({
                rank: i + 1,
                sample: <span className="font-mono">{sample}</span>,
                score: score.toFixed(4),
              }))}
            />
            {rows.length > 25 && (
              <p className="mt-2 text-[11.5px] text-ink-mute">Showing top 25 of {rows.length} samples.</p>
            )}
          </Panel>
        )}
      </div>
    </div>
  );
}
