import { useEffect, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { PageHeader, Panel, PrimaryButton, Stat, EmptyState, ErrorState, DataTable } from "../components/ui";
import { BoxPlot } from "../components/BoxPlot";
import { ExportButton } from "../components/ExportButton";

export function ComparePage() {
  const chartRef = useRef<SVGSVGElement>(null);
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId] = useState<string>("");
  const [geneInput, setGeneInput] = useState("MYCN");
  const [groupColumn, setGroupColumn] = useState<string>("");
  const [query, setQuery] = useState<{ datasetId: string; gene: string; groupColumn: string } | null>(null);

  // Default to the first dataset, and its first group column, once the
  // registry responds — nothing here names a specific dataset_id.
  useEffect(() => {
    if (!datasetId && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);

  const activeDataset = datasets.find((d) => d.dataset_id === datasetId);

  useEffect(() => {
    if (activeDataset && activeDataset.group_columns.length > 0) {
      setGroupColumn((prev) =>
        activeDataset.group_columns.includes(prev) ? prev : activeDataset.group_columns[0],
      );
    }
  }, [activeDataset]);

  const { data, isFetching, error } = useQuery({
    queryKey: ["compare", query],
    queryFn: () => api.compare(query!.datasetId, query!.gene, query!.groupColumn),
    enabled: !!query,
  });

  const run = () => {
    const gene = geneInput.trim().toUpperCase();
    if (gene && datasetId && groupColumn) setQuery({ datasetId, gene, groupColumn });
  };

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Expression Compare"
        title="Gene expression across sample groups"
        description="Box plot with pairwise Mann–Whitney tests and FDR correction, matching the paper's group-comparison methodology."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <Panel title="Query">
          <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:flex-wrap">
            <div className="min-w-[14ch] flex-1">
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
            <div className="min-w-[12ch] flex-1">
              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Gene symbol
              </label>
              <input
                value={geneInput}
                onChange={(e) => setGeneInput(e.target.value)}
                onKeyDown={(e) => e.key === "Enter" && run()}
                placeholder="e.g. MYCN"
                className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 font-mono text-[13px] italic text-ink outline-none focus:border-accent"
              />
            </div>
            <div className="min-w-[14ch] flex-1">
              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Group by
              </label>
              <select
                value={groupColumn}
                onChange={(e) => setGroupColumn(e.target.value)}
                disabled={!activeDataset || activeDataset.group_columns.length === 0}
                className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent disabled:opacity-50"
              >
                {activeDataset?.group_columns.map((g) => (
                  <option key={g} value={g}>
                    {g}
                  </option>
                ))}
              </select>
            </div>
            <PrimaryButton onClick={run} loading={isFetching} disabled={!datasetId || !groupColumn}>
              Run comparison
            </PrimaryButton>
          </div>
          {activeDataset && activeDataset.group_columns.length === 0 && (
            <p className="mt-2.5 text-[12px] text-warn">
              This dataset exposes no grouping columns yet — nothing to compare against.
            </p>
          )}
        </Panel>

        {error && <ErrorState message={(error as Error).message} />}

        {!query && !error && (
          <EmptyState
            icon={
              <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                <path d="M4 20V10M12 20V4M20 20V14" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
              </svg>
            }
            title="No query run yet"
            description="Choose a dataset, a gene symbol, and a grouping, then run a comparison to see the distribution and statistics."
          />
        )}

        {data && (
          <>
            <Panel
              title={`${data.gene} by ${data.group_column}`}
              action={
                <div className="flex items-center gap-2.5">
                  <span className="font-mono text-[10.5px] text-ink-mute">n = {data.points.length}</span>
                  <ExportButton
                    svgRef={chartRef}
                    filename={`${data.gene}-by-${data.group_column}`}
                    title={`${data.gene} by ${data.group_column}`}
                    subtitle={`n = ${data.points.length} · Mann–Whitney with FDR correction`}
                    statLines={[
                      `Kruskal-Wallis p = ${
                        data.kruskal_wallis.p_value != null ? data.kruskal_wallis.p_value.toExponential(2) : "—"
                      }`,
                      ...data.pairwise_tests.map(
                        (t) =>
                          `${t.group_a} vs ${t.group_b}: p = ${t.p_value < 0.001 ? t.p_value.toExponential(2) : t.p_value.toFixed(4)}`,
                      ),
                    ]}
                  />
                </div>
              }
            >
              <BoxPlot points={data.points} valueLabel={`${data.gene} (expression)`} svgRef={chartRef} />
            </Panel>

            <Panel title="Statistics">
              <div className="mb-4 flex gap-8">
                <Stat
                  label="Kruskal-Wallis p"
                  value={data.kruskal_wallis.p_value != null ? data.kruskal_wallis.p_value.toExponential(2) : "—"}
                  tone={data.kruskal_wallis.p_value != null && data.kruskal_wallis.p_value < 0.05 ? "hot" : "neutral"}
                />
              </div>
              {data.pairwise_tests.length > 0 ? (
                <DataTable
                  columns={[
                    { key: "pair", label: "Comparison" },
                    { key: "n", label: "n (a / b)", align: "right" },
                    { key: "p", label: "p-value", align: "right" },
                    { key: "q", label: "FDR q", align: "right" },
                  ]}
                  rows={data.pairwise_tests.map((t) => ({
                    pair: `${t.group_a} vs ${t.group_b}`,
                    n: `${t.n_a} / ${t.n_b}`,
                    p: t.p_value < 0.001 ? t.p_value.toExponential(2) : t.p_value.toFixed(4),
                    q:
                      t.q_value != null
                        ? t.q_value < 0.001
                          ? t.q_value.toExponential(2)
                          : t.q_value.toFixed(4)
                        : "—",
                  }))}
                />
              ) : (
                <p className="text-[12.5px] text-ink-mute">Only one group present — no pairwise test to show.</p>
              )}
            </Panel>
          </>
        )}
      </div>
    </div>
  );
}
