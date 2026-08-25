import { useEffect, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, EmptyState, ErrorState, DataTable } from "../components/ui";

const PRESETS: Record<string, string[]> = {
  "ETP-TF5": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
  "ZMIZ1-5": ["MEF2C", "BCL2", "MYB", "MYCN", "ZMIZ1"],
};

/**
 * A compact histogram of every score in the cohort, with the range covered
 * by the visible top-N shaded. Without it the page shows 25 numbers out of
 * a few hundred and gives no sense of whether the top scores are a
 * distinct population or the tail of one continuous distribution -- which
 * is the actual question when picking candidates.
 */
function ScoreDistribution({ scores, shown }: { scores: number[]; shown: number }) {
  if (scores.length < 4) return null;
  const min = Math.min(...scores);
  const max = Math.max(...scores);
  if (!(max > min)) return null;

  const BINS = 36;
  const counts = new Array(BINS).fill(0);
  for (const s of scores) {
    const idx = Math.min(BINS - 1, Math.floor(((s - min) / (max - min)) * BINS));
    counts[idx] += 1;
  }
  const peak = Math.max(...counts);
  // Threshold score of the last visible row -- everything at or above it is
  // in the table.
  const cutoff = [...scores].sort((a, b) => b - a)[shown - 1] ?? min;
  const cutoffFrac = (cutoff - min) / (max - min);

  return (
    <div className="mb-4">
      <div className="mb-1.5 flex items-baseline justify-between">
        <span className="font-mono text-[10px] uppercase tracking-wider text-ink-mute">
          Score distribution · {scores.length} samples
        </span>
        <span className="font-mono text-[10px] text-ink-mute">
          shaded = top {shown} shown below
        </span>
      </div>
      <div className="flex h-14 items-end gap-[1.5px] rounded-[3px] border border-rule bg-ground px-2 pt-2">
        {counts.map((c, i) => {
          // A bin is "in the table" if its right edge is past the cutoff --
          // using the left edge (i / BINS) excluded the bin the cutoff
          // itself falls inside, so the shaded region under-covered the
          // labeled top-N by 1-2 samples.
          const inTable = (i + 1) / BINS > cutoffFrac;
          return (
            <span
              key={i}
              title={`${c} sample${c === 1 ? "" : "s"}`}
              className={`flex-1 rounded-t-[1px] ${inTable ? "bg-accent" : "bg-rule-firm"}`}
              style={{ height: `${peak ? Math.max((c / peak) * 100, c > 0 ? 4 : 0) : 0}%` }}
            />
          );
        })}
      </div>
      <div className="mt-1 flex justify-between font-mono text-[9.5px] text-ink-mute">
        <span>{min.toFixed(3)}</span>
        <span>{max.toFixed(3)}</span>
      </div>
    </div>
  );
}

export function SignaturePage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);

  const [genes, setGenes] = useState<string[]>(PRESETS["ETP-TF5"]);
  const [method, setMethod] = useState<"auc" | "log2_mean">("auc");
  const [activePreset, setActivePreset] = useState<string | null>("ETP-TF5");
  const [query, setQuery] = useState<{ datasetId: string; genes: string[]; method: "auc" | "log2_mean" } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["signature", query],
    queryFn: () => api.signatureScore(query!.datasetId, query!.genes, query!.method),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) setQuery({ datasetId, genes, method });
  };

  // A results table from the PREVIOUS query is only safe to keep on screen
  // while it still describes the current form state. The moment the user
  // changes dataset or method, `data` reflects a query that no longer
  // matches what's selected -- e.g. 186 DepMap rows staying on screen,
  // still labeled "n = 186", after switching to TARGET, with nothing
  // marking them stale. Comparing against the query that produced `data`
  // (not the live form state) is what lets a fresh run clear this safely.
  const resultsAreStale = !!data && !!query && (query.datasetId !== datasetId || query.method !== method);

  const rows = data
    ? Object.entries(data.scores)
        .filter(([, v]) => v != null && !Number.isNaN(v))
        .sort((a, b) => b[1] - a[1])
    : [];
  // Bars are scaled across the whole cohort, not just the visible 25, so
  // the top-25 view doesn't imply the full range is on screen.
  const scoreMax = rows.length ? rows[0][1] : 1;
  const scoreMin = rows.length ? rows[rows.length - 1][1] : 0;

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Signature Score"
        title="Score any gene set across samples"
        description="AUCell-equivalent rank-based scoring, or a simpler log2-mean. Enter a custom gene set or start from a preset."
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
                    {m === "auc" ? "AUCell" : "log2 mean"}
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
            {resultsAreStale && (
              <p className="mb-3 rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink">
                Dataset or method changed since this ran — these results are from{" "}
                {datasets.find((d) => d.dataset_id === query!.datasetId)?.display_name ?? query!.datasetId}
                {" · "}
                {query!.method === "auc" ? "AUCell" : "log2 mean"}. Run again to refresh.
              </p>
            )}
            <div className={resultsAreStale ? "opacity-50" : undefined}>
            <ScoreDistribution scores={rows.map(([, v]) => v)} shown={Math.min(25, rows.length)} />
            <DataTable
              columns={[
                { key: "rank", label: "#", align: "right" },
                { key: "sample", label: "Sample" },
                { key: "bar", label: "Score" },
                { key: "score", label: "", align: "right" },
              ]}
              rows={rows.slice(0, 25).map(([sample, score], i) => ({
                rank: i + 1,
                sample: <span className="font-mono">{sample}</span>,
                // An inline bar makes the spread between ranks readable at a
                // glance; a column of 4-decimal numbers alone does not show
                // whether rank 1 and rank 25 differ by a little or a lot.
                bar: (
                  <span className="flex items-center">
                    <span
                      className="inline-block h-[7px] rounded-[2px] bg-accent"
                      style={{
                        width: `${scoreMax > scoreMin ? ((score - scoreMin) / (scoreMax - scoreMin)) * 100 : 100}%`,
                        minWidth: "2px",
                      }}
                    />
                  </span>
                ),
                score: score.toFixed(4),
              }))}
            />
            {rows.length > 25 && (
              <p className="mt-2 text-[11.5px] text-ink-mute">
                Showing top 25 of {rows.length} samples — the distribution above covers all {rows.length}.
              </p>
            )}
            </div>
          </Panel>
        )}
      </div>
    </div>
  );
}
