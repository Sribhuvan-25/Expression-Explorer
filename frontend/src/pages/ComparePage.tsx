import { useEffect, useMemo, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api, type CompareMultiDatasetResult } from "../lib/api";
import { PageHeader, Panel, PrimaryButton, Stat, EmptyState, ErrorState, DataTable } from "../components/ui";
import { BoxPlot } from "../components/BoxPlot";
import { ExportButton } from "../components/ExportButton";
import { GroupFilter, useGroupFilter } from "../components/GroupFilter";
import { GeneAnnotation } from "../components/GeneAnnotation";
import { exportRowsAsCsv } from "../lib/exportTable";

// Every group-by column any registered dataset has ever declared, so the
// "Group by" dropdown offers a column even before a dataset exposing it is
// checked -- picking the column first, then seeing which datasets qualify,
// is the natural order (the alternative, deriving options only from
// checked datasets, means unchecking everything empties the dropdown).
function allGroupColumns(datasets: { group_columns: string[] }[]): string[] {
  const seen = new Set<string>();
  for (const d of datasets) for (const g of d.group_columns) seen.add(g);
  return Array.from(seen);
}

function DatasetPanel({ result }: { result: CompareMultiDatasetResult }) {
  // One ref per panel, not a ref shared across every dataset's chart --
  // each DatasetPanel mounts its own BoxPlot, and a shared ref would just
  // get overwritten by whichever panel mounts last, silently breaking
  // export for every dataset above it.
  const chartRef = useRef<SVGSVGElement>(null);
  const { groupCounts, selected: selectedGroups, setSelected: setSelectedGroups } = useGroupFilter(
    result.points ?? [],
  );
  const filteredPoints = useMemo(
    () => (result.points ?? []).filter((p) => selectedGroups.has(p.group)),
    [result.points, selectedGroups],
  );

  if (result.skipped) {
    return (
      <Panel title={result.display_name}>
        <p className="text-[12.5px] text-ink-mute">Skipped — {result.skip_reason}.</p>
      </Panel>
    );
  }

  const points = result.points ?? [];
  const pairwise = result.pairwise_tests ?? [];

  return (
    <Panel
      title={result.display_name}
      action={
        <div className="flex items-center gap-2.5">
          <span className="font-mono text-[10px] uppercase tracking-wider text-ink-mute">
            {result.assay_type?.replace(/_/g, "-")} · {result.expression_unit?.replace(/_/g, " ")}
          </span>
          <span className="font-mono text-[10.5px] text-ink-mute">
            n = {filteredPoints.length}
            {filteredPoints.length !== points.length && ` of ${points.length}`}
          </span>
          <ExportButton
            svgRef={chartRef}
            filename={`${result.dataset_id}-by-group`}
            title={result.display_name}
            subtitle={`n = ${filteredPoints.length} · Mann–Whitney with FDR correction`}
            statLines={[
              `Kruskal-Wallis p = ${
                result.kruskal_wallis?.p_value != null ? result.kruskal_wallis.p_value.toExponential(2) : "—"
              }`,
              ...pairwise.map(
                (t) =>
                  `${t.group_a} vs ${t.group_b}: p = ${t.p_value < 0.001 ? t.p_value.toExponential(2) : t.p_value.toFixed(4)}`,
              ),
            ]}
          />
        </div>
      }
    >
      <div className="mb-3">
        <GroupFilter groupCounts={groupCounts} selected={selectedGroups} onChange={setSelectedGroups} />
      </div>
      {filteredPoints.length > 0 ? (
        <BoxPlot points={filteredPoints} valueLabel="expression" svgRef={chartRef} />
      ) : (
        <p className="py-8 text-center text-[12.5px] text-ink-mute">
          No groups selected — check at least one group above to see the chart.
        </p>
      )}
      {(result.n_excluded ?? 0) > 0 && (
        <p className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
          This grouping covers {points.length} of {result.n_dataset_total} samples in the dataset —{" "}
          {result.n_excluded} excluded ({result.exclusion_reason}).
        </p>
      )}

      <div className="mt-4 border-t border-rule pt-3">
        <div className="mb-3 flex gap-8">
          <Stat
            label="Kruskal-Wallis p"
            value={result.kruskal_wallis?.p_value != null ? result.kruskal_wallis.p_value.toExponential(2) : "—"}
            tone={result.kruskal_wallis?.p_value != null && result.kruskal_wallis.p_value < 0.05 ? "hot" : "neutral"}
          />
        </div>
        {pairwise.length > 0 ? (
          <DataTable
            columns={[
              { key: "pair", label: "Comparison" },
              { key: "n", label: "n (a / b)", align: "right" },
              { key: "p", label: "p-value", align: "right" },
              { key: "q", label: "FDR q", align: "right" },
            ]}
            rows={pairwise.map((t) => ({
              pair: `${t.group_a} vs ${t.group_b}`,
              n: `${t.n_a} / ${t.n_b}`,
              p: t.p_value < 0.001 ? t.p_value.toExponential(2) : t.p_value.toFixed(4),
              q: t.q_value != null ? (t.q_value < 0.001 ? t.q_value.toExponential(2) : t.q_value.toFixed(4)) : "—",
            }))}
          />
        ) : (
          <p className="text-[12.5px] text-ink-mute">Only one group present — no pairwise test to show.</p>
        )}
      </div>
    </Panel>
  );
}

type Tab = "single-gene" | "differential";

function SingleGeneCompareTab() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [geneInput, setGeneInput] = useState("MYCN");
  const [groupColumn, setGroupColumn] = useState<string>("");
  // null = "follow the default" (every dataset that exposes groupColumn).
  // Becomes an explicit Set the moment the user checks/unchecks anything,
  // so a manual narrowing survives until they change groupColumn again.
  const [selectedDatasetIds, setSelectedDatasetIds] = useState<Set<string> | null>(null);
  const [query, setQuery] = useState<{ gene: string; groupColumn: string; datasetIds: string[] } | null>(null);

  const groupColumnOptions = useMemo(() => allGroupColumns(datasets), [datasets]);

  useEffect(() => {
    if (!groupColumn && groupColumnOptions.length > 0) setGroupColumn(groupColumnOptions[0]);
  }, [groupColumnOptions, groupColumn]);

  const compatibleDatasets = useMemo(
    () => datasets.filter((d) => d.group_columns.includes(groupColumn)),
    [datasets, groupColumn],
  );
  const incompatibleDatasets = useMemo(
    () => datasets.filter((d) => !d.group_columns.includes(groupColumn)),
    [datasets, groupColumn],
  );

  // Resetting to "follow the default" whenever groupColumn changes is
  // deliberate: a dataset checked under the old grouping may not even
  // expose the new one, and silently carrying forward a selection that no
  // longer makes sense is worse than re-defaulting to "everything
  // compatible" and letting the user narrow again.
  useEffect(() => {
    setSelectedDatasetIds(null);
  }, [groupColumn]);

  const effectiveSelected = selectedDatasetIds ?? new Set(compatibleDatasets.map((d) => d.dataset_id));

  const toggleDataset = (id: string) => {
    const next = new Set(effectiveSelected);
    if (next.has(id)) next.delete(id);
    else next.add(id);
    setSelectedDatasetIds(next);
  };

  const { data, isFetching, error } = useQuery({
    queryKey: ["compare-multi", query],
    queryFn: () => api.compareMulti(query!.gene, query!.groupColumn, query!.datasetIds),
    enabled: !!query,
  });

  const run = () => {
    const gene = geneInput.trim().toUpperCase();
    const datasetIds = Array.from(effectiveSelected);
    if (gene && groupColumn && datasetIds.length > 0) setQuery({ gene, groupColumn, datasetIds });
  };

  const resultsAreStale =
    !!data &&
    !!query &&
    (query.gene !== geneInput.trim().toUpperCase() ||
      query.groupColumn !== groupColumn ||
      query.datasetIds.length !== effectiveSelected.size ||
      query.datasetIds.some((id) => !effectiveSelected.has(id)));

  const shownResults = data?.datasets ?? [];
  const runnable = shownResults.filter((r) => !r.skipped);
  const skipped = shownResults.filter((r) => r.skipped);

  return (
    <>
        <Panel title="Query">
          <div className="flex flex-col gap-4 @sm:flex-row @sm:flex-wrap @sm:items-end">
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
                disabled={groupColumnOptions.length === 0}
                className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent disabled:opacity-50"
              >
                {groupColumnOptions.map((g) => (
                  <option key={g} value={g}>
                    {g}
                  </option>
                ))}
              </select>
            </div>
            <PrimaryButton onClick={run} loading={isFetching} disabled={!groupColumn || effectiveSelected.size === 0}>
              Run comparison
            </PrimaryButton>
          </div>

          <div className="mt-4 border-t border-rule pt-3">
            <p className="mb-2 font-mono text-[10px] uppercase tracking-wider text-ink-mute">
              Datasets — defaults to every dataset that exposes "{groupColumn || "…"}"
            </p>
            <div className="flex flex-col gap-1">
              {compatibleDatasets.map((d) => (
                <label
                  key={d.dataset_id}
                  className="flex cursor-pointer items-center gap-2 rounded-[3px] px-1.5 py-1 text-[12.5px] hover:bg-surface-hover"
                >
                  <input
                    type="checkbox"
                    checked={effectiveSelected.has(d.dataset_id)}
                    onChange={() => toggleDataset(d.dataset_id)}
                    className="h-3.5 w-3.5 shrink-0 accent-accent"
                  />
                  <span className="min-w-0 flex-1 text-ink-soft">{d.display_name}</span>
                  <span className="shrink-0 font-mono text-[10.5px] text-ink-mute">
                    {d.n_samples != null && `${d.n_samples} samples`}
                    {d.assay_type && ` · ${d.assay_type.replace(/_/g, "-")}`}
                  </span>
                </label>
              ))}
            </div>
            {incompatibleDatasets.length > 0 && (
              <p className="mt-2 text-[11.5px] text-ink-mute">
                Not shown (no "{groupColumn}" column):{" "}
                {incompatibleDatasets.map((d) => d.display_name).join(", ")}.
              </p>
            )}
          </div>

          {groupColumnOptions.length === 0 && (
            <p className="mt-2.5 text-[12px] text-warn">No dataset exposes any grouping columns yet.</p>
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
            description="Enter a gene symbol and a grouping, then run a comparison to see the distribution and statistics for every compatible dataset."
          />
        )}

        {data && (
          <>
            {resultsAreStale && (
              <p className="rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink">
                Query settings changed since this ran — these results are from {query!.gene} by {query!.groupColumn}{" "}
                across {query!.datasetIds.length} dataset{query!.datasetIds.length === 1 ? "" : "s"}. Run again to
                refresh.
              </p>
            )}
            <div className={resultsAreStale ? "flex flex-col gap-5 opacity-50" : "flex flex-col gap-5"}>
              {/* Keyed off the RUN query's gene, not the live input -- an
                  annotation lookup on every keystroke would be wasteful and
                  would flicker between genes as the user types. */}
              <GeneAnnotation gene={query!.gene} />
              {runnable.length === 0 && skipped.length > 0 && (
                <EmptyState
                  icon={
                    <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                      <path d="M4 20V10M12 20V4M20 20V14" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
                    </svg>
                  }
                  title="No dataset returned results"
                  description="Every selected dataset was skipped — see the reasons above."
                />
              )}
              {shownResults.map((r) => (
                <DatasetPanel key={r.dataset_id} result={r} />
              ))}
            </div>
          </>
        )}
    </>
  );
}

function DifferentialGenesTab() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);

  const activeDataset = datasets.find((d) => d.dataset_id === datasetId);
  const [groupColumn, setGroupColumn] = useState<string>("");
  const [groupA, setGroupA] = useState<string>("");
  const [groupB, setGroupB] = useState<string>("");

  useEffect(() => {
    if (activeDataset && activeDataset.group_columns.length > 0) {
      setGroupColumn((prev) => (activeDataset.group_columns.includes(prev) ? prev : activeDataset.group_columns[0]));
    }
  }, [activeDataset]);

  // The actual group VALUES for the chosen column (e.g. "ETP" / "non-ETP")
  // aren't in the dataset registry response (that only lists column NAMES
  // like "etp_status") -- /group-values reads them directly off the
  // dataset's own samples, independent of any particular gene existing.
  const { data: groupValuesResult } = useQuery({
    queryKey: ["group-values", datasetId, groupColumn],
    queryFn: () => api.groupValues(datasetId, groupColumn),
    enabled: !!datasetId && !!groupColumn,
    retry: false,
  });
  const groupOptions = useMemo(
    () => groupValuesResult?.values.map((v) => v.value) ?? [],
    [groupValuesResult],
  );
  const countByGroup = useMemo(
    () => new Map((groupValuesResult?.values ?? []).map((v) => [v.value, v.n])),
    [groupValuesResult],
  );

  useEffect(() => {
    if (groupOptions.length >= 2) {
      setGroupA((prev) => (groupOptions.includes(prev) ? prev : groupOptions[0]));
      setGroupB((prev) => (groupOptions.includes(prev) && prev !== groupOptions[0] ? prev : groupOptions[1]));
    }
  }, [groupOptions]);

  const [query, setQuery] = useState<{ datasetId: string; groupColumn: string; groupA: string; groupB: string } | null>(
    null,
  );

  const { data, isFetching, error } = useQuery({
    queryKey: ["differential", query],
    queryFn: () => api.differential(query!.datasetId, query!.groupColumn, query!.groupA, query!.groupB),
    enabled: !!query,
  });

  const run = () => {
    if (datasetId && groupColumn && groupA && groupB && groupA !== groupB) {
      setQuery({ datasetId, groupColumn, groupA, groupB });
    }
  };

  return (
    <>
      <Panel title="Query">
        <div className="flex flex-col gap-4 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[20ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Dataset</label>
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
          <div className="min-w-[14ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Group by</label>
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
          <div className="min-w-[12ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Group A</label>
            <select
              value={groupA}
              onChange={(e) => setGroupA(e.target.value)}
              disabled={groupOptions.length < 2}
              className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent disabled:opacity-50"
            >
              {groupOptions.map((g) => (
                <option key={g} value={g}>
                  {g} (n={countByGroup.get(g) ?? "?"})
                </option>
              ))}
            </select>
          </div>
          <div className="min-w-[12ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Group B</label>
            <select
              value={groupB}
              onChange={(e) => setGroupB(e.target.value)}
              disabled={groupOptions.length < 2}
              className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent disabled:opacity-50"
            >
              {groupOptions.map((g) => (
                <option key={g} value={g}>
                  {g} (n={countByGroup.get(g) ?? "?"})
                </option>
              ))}
            </select>
          </div>
          <PrimaryButton onClick={run} loading={isFetching} disabled={!groupA || !groupB || groupA === groupB}>
            Rank genes
          </PrimaryButton>
        </div>
        {groupA && groupA === groupB && (
          <p className="mt-2.5 text-[12px] text-warn">Group A and Group B must be different.</p>
        )}
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <path d="M6 6h5M6 12h9M6 18h13" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
            </svg>
          }
          title="No ranking run yet"
          description="Choose a dataset, a grouping, and two groups to compare — every gene in the dataset is ranked by how strongly it separates them."
        />
      )}

      {data && (
        <Panel
          title={`Top ${data.genes.length} genes: ${data.group_a} vs ${data.group_b}`}
          action={
            <button
              type="button"
              onClick={() =>
                exportRowsAsCsv(
                  `${data.group_a}-vs-${data.group_b}-differential`,
                  [
                    { key: "rank", label: "#" },
                    { key: "symbol", label: "gene" },
                    { key: "p_value", label: "p-value" },
                    { key: "q_value", label: "FDR q" },
                    { key: "median_diff", label: `median(${data.group_a}) - median(${data.group_b})` },
                  ],
                  data.genes.map((g, i) => ({ ...g, rank: i + 1 })),
                )
              }
              title="Download as CSV"
              className="flex items-center gap-1 rounded-[3px] border border-rule-firm bg-surface px-2 py-1 font-mono text-[10.5px] text-ink-mute transition-colors hover:border-accent hover:text-accent"
            >
              <svg width="10" height="10" viewBox="0 0 10 10" fill="none">
                <path d="M5 1v6M2.5 5L5 7.5L7.5 5M1.5 8.5h7" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" />
              </svg>
              CSV
            </button>
          }
        >
          <p className="mb-3 text-[12px] text-ink-mute">
            n = {data.n_a} ({data.group_a}) vs {data.n_b} ({data.group_b}) · {data.n_genes_tested.toLocaleString()}{" "}
            genes tested, Mann–Whitney with genome-wide FDR correction.
          </p>
          <DataTable
            columns={[
              { key: "rank", label: "#", align: "right" },
              { key: "gene", label: "Gene" },
              { key: "p", label: "p-value", align: "right" },
              { key: "q", label: "FDR q", align: "right" },
              { key: "diff", label: `Δ median (${data.group_a} − ${data.group_b})`, align: "right" },
            ]}
            rows={data.genes.map((g, i) => ({
              rank: i + 1,
              gene: <span className="font-mono">{g.symbol}</span>,
              p: g.p_value < 0.001 ? g.p_value.toExponential(2) : g.p_value.toFixed(4),
              q: g.q_value < 0.001 ? g.q_value.toExponential(2) : g.q_value.toFixed(4),
              diff: g.median_diff.toFixed(3),
            }))}
          />
        </Panel>
      )}
    </>
  );
}

export function ComparePage() {
  const [tab, setTab] = useState<Tab>("single-gene");

  return (
    <div className="mx-auto max-w-[1100px]">
      <PageHeader
        eyebrow="Expression Compare"
        title="Gene expression across sample groups"
        description="Box plot with pairwise Mann–Whitney tests and FDR correction, faceted per dataset — values are never pooled across assay types or units."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <div className="flex overflow-hidden rounded-[3px] border border-rule self-start">
          {([
            { key: "single-gene" as const, label: "Single gene" },
            { key: "differential" as const, label: "Top differential genes" },
          ]).map((t) => (
            <button
              key={t.key}
              type="button"
              onClick={() => setTab(t.key)}
              className={`px-3 py-2 text-[12.5px] font-medium transition-colors ${
                tab === t.key ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
              }`}
            >
              {t.label}
            </button>
          ))}
        </div>

        {/* Both tabs stay mounted, toggled with CSS rather than
            conditionally rendered -- a conditional render unmounts the
            inactive tab and wipes its local state, exactly the bug found
            and fixed on the Correlation pane in the last QA round (see
            METHODS.md 7.6). Applying the same fix here from the start
            rather than re-discovering it. */}
        <div className={tab === "single-gene" ? "contents" : "hidden"}>
          <SingleGeneCompareTab />
        </div>
        <div className={tab === "differential" ? "contents" : "hidden"}>
          <DifferentialGenesTab />
        </div>
      </div>
    </div>
  );
}
