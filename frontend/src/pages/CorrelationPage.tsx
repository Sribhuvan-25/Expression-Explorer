import { useEffect, useMemo, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { DataTable, EmptyState, ErrorState, GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, Stat } from "../components/ui";
import { ScatterPlot } from "../components/ScatterPlot";
import { ExportButton } from "../components/ExportButton";
import { exportRowsAsCsv } from "../lib/exportTable";

type Tab = "correlation" | "network" | "pca";

const PCA_PRESETS: Record<string, string[]> = {
  "ETP-TF5": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
};

function CorrelationTab({ datasetId }: { datasetId: string }) {
  const chartRef = useRef<SVGSVGElement>(null);
  const [geneA, setGeneA] = useState("MYCN");
  const [geneB, setGeneB] = useState("MYC");
  const [method, setMethod] = useState<"pearson" | "spearman" | "kendall">("pearson");
  const [query, setQuery] = useState<{ datasetId: string; geneA: string; geneB: string; method: typeof method } | null>(
    null,
  );

  const { data, isFetching, error } = useQuery({
    queryKey: ["correlation", query],
    queryFn: () => api.correlation(query!.datasetId, query!.geneA, query!.geneB, query!.method),
    enabled: !!query,
  });

  const run = () => {
    const a = geneA.trim().toUpperCase();
    const b = geneB.trim().toUpperCase();
    if (a && b && datasetId) setQuery({ datasetId, geneA: a, geneB: b, method });
  };

  return (
    <>
      <Panel title="Query">
        <div className="flex flex-col gap-4 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Gene A</label>
            <input
              value={geneA}
              onChange={(e) => setGeneA(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 font-mono text-[13px] italic text-ink outline-none focus:border-accent"
            />
          </div>
          <div className="min-w-[12ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Gene B</label>
            <input
              value={geneB}
              onChange={(e) => setGeneB(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 font-mono text-[13px] italic text-ink outline-none focus:border-accent"
            />
          </div>
          <div className="w-full @sm:w-auto @sm:min-w-[27ch] @sm:flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Method</label>
            {/* "Kendall" is the widest label of the three -- a container that
                can be squeezed by its @sm:flex-row siblings (as the prior
                overflow-hidden version was) clips it instead of reflowing,
                since none of the three buttons can shrink past their text.
                min-w-[27ch] keeps the group from being squeezed that far in
                the first place, and dropping overflow-hidden means if it
                ever is squeezed anyway, the buttons wrap onto a second row
                instead of silently cutting off "Kendall" mid-word. */}
            <div className="flex flex-wrap overflow-hidden rounded-[3px] border border-rule">
              {(["pearson", "spearman", "kendall"] as const).map((m) => (
                <button
                  key={m}
                  type="button"
                  onClick={() => setMethod(m)}
                  className={`flex-1 whitespace-nowrap px-2 py-2 font-mono text-[11.5px] capitalize transition-colors ${
                    method === m ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
                  }`}
                >
                  {m}
                </button>
              ))}
            </div>
          </div>
          <PrimaryButton
            onClick={run}
            loading={isFetching}
            disabled={!geneA.trim() || !geneB.trim() || geneA.trim().toUpperCase() === geneB.trim().toUpperCase()}
          >
            Run correlation
          </PrimaryButton>
        </div>
        {/* A gene correlated against itself isn't a meaningful query (the
            backend now rejects it cleanly rather than crashing -- see
            METHODS.md 7.6 -- but disabling Run here avoids the round trip
            entirely once the mistake is visible on screen). */}
        {geneA.trim() && geneA.trim().toUpperCase() === geneB.trim().toUpperCase() && (
          <p className="mt-2.5 text-[12px] text-warn">Gene A and Gene B must be different.</p>
        )}
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <circle cx="6" cy="17" r="1.6" fill="currentColor" />
              <circle cx="10" cy="11" r="1.6" fill="currentColor" />
              <circle cx="14" cy="14" r="1.6" fill="currentColor" />
              <circle cx="18" cy="7" r="1.6" fill="currentColor" />
            </svg>
          }
          title="No correlation run yet"
          description="Enter two gene symbols to see how their expression tracks across this dataset's samples."
        />
      )}

      {data && (
        <Panel
          title={`${data.gene_a} vs ${data.gene_b}`}
          action={
            <div className="flex items-center gap-2.5">
              <span className="font-mono text-[10.5px] text-ink-mute">n = {data.n}</span>
              <ExportButton
                svgRef={chartRef}
                filename={`${data.gene_a}-vs-${data.gene_b}-correlation`}
                title={`${data.gene_a} vs ${data.gene_b}`}
                subtitle={`n = ${data.n} · ${data.method}`}
                statLines={[
                  `${data.method} r = ${data.coefficient.toFixed(3)}`,
                  `p = ${data.p_value < 0.001 ? data.p_value.toExponential(2) : data.p_value.toFixed(4)}`,
                ]}
              />
            </div>
          }
        >
          <ScatterPlot points={data.points} xLabel={data.gene_a} yLabel={data.gene_b} svgRef={chartRef} />
          <div className="mt-4 flex gap-8 border-t border-rule pt-3">
            <Stat
              label={`${data.method[0].toUpperCase()}${data.method.slice(1)} r`}
              value={data.coefficient.toFixed(3)}
              tone={Math.abs(data.coefficient) > 0.5 ? "accent" : "neutral"}
            />
            <Stat
              label="p-value"
              value={data.p_value < 0.001 ? data.p_value.toExponential(2) : data.p_value.toFixed(4)}
              tone={data.p_value < 0.05 ? "hot" : "neutral"}
            />
          </div>
        </Panel>
      )}
    </>
  );
}

function NetworkTab({ datasetId }: { datasetId: string }) {
  const [gene, setGene] = useState("MYCN");
  const [topN, setTopN] = useState(20);
  const [method, setMethod] = useState<"pearson" | "spearman">("pearson");
  const [direction, setDirection] = useState<"positive" | "negative">("positive");
  const [query, setQuery] = useState<{
    datasetId: string;
    gene: string;
    topN: number;
    method: typeof method;
    direction: typeof direction;
  } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["co-expression", query],
    queryFn: () => api.coExpression(query!.datasetId, query!.gene, { topN: query!.topN, method: query!.method, direction: query!.direction }),
    enabled: !!query,
  });

  const run = () => {
    const g = gene.trim().toUpperCase();
    if (g && datasetId) setQuery({ datasetId, gene: g, topN, method, direction });
  };

  return (
    <>
      <Panel title="Query">
        <div className="flex flex-col gap-4 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Gene</label>
            <input
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 font-mono text-[13px] italic text-ink outline-none focus:border-accent"
            />
          </div>
          <div className="min-w-[10ch]">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Top N</label>
            <input
              type="number"
              min={1}
              max={200}
              value={topN}
              onChange={(e) => setTopN(Number(e.target.value))}
              className="w-20 rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent"
            />
          </div>
          <div className="min-w-[14ch]">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Method</label>
            <div className="flex overflow-hidden rounded-[3px] border border-rule">
              {(["pearson", "spearman"] as const).map((m) => (
                <button
                  key={m}
                  type="button"
                  onClick={() => setMethod(m)}
                  className={`px-3 py-2 font-mono text-[11.5px] capitalize transition-colors ${
                    method === m ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
                  }`}
                >
                  {m}
                </button>
              ))}
            </div>
          </div>
          <div className="min-w-[16ch]">
            <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">Direction</label>
            <div className="flex overflow-hidden rounded-[3px] border border-rule">
              {(["positive", "negative"] as const).map((d) => (
                <button
                  key={d}
                  type="button"
                  onClick={() => setDirection(d)}
                  className={`px-3 py-2 font-mono text-[11.5px] capitalize transition-colors ${
                    direction === d ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
                  }`}
                >
                  {d}
                </button>
              ))}
            </div>
          </div>
          <PrimaryButton onClick={run} loading={isFetching} disabled={!gene.trim()}>
            Build network
          </PrimaryButton>
        </div>
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <circle cx="12" cy="12" r="2.5" stroke="currentColor" strokeWidth="1.5" />
              <circle cx="4" cy="6" r="1.8" stroke="currentColor" strokeWidth="1.5" />
              <circle cx="20" cy="6" r="1.8" stroke="currentColor" strokeWidth="1.5" />
              <circle cx="4" cy="18" r="1.8" stroke="currentColor" strokeWidth="1.5" />
              <circle cx="20" cy="18" r="1.8" stroke="currentColor" strokeWidth="1.5" />
              <path d="M6 7.5L10.5 10.5M18 7.5L13.5 10.5M6 16.5L10.5 13.5M18 16.5L13.5 13.5" stroke="currentColor" strokeWidth="1.3" />
            </svg>
          }
          title="No network built yet"
          description="Enter a gene to find the samples-wide genes most (or least) correlated with it in this dataset."
        />
      )}

      {data && (
        <Panel
          title={`Top ${data.network.length} ${data.direction === "positive" ? "co-expressed with" : "anti-correlated with"} ${data.gene}`}
          action={
            <button
              type="button"
              onClick={() =>
                exportRowsAsCsv(
                  `${data.gene}-${data.direction}-network`,
                  [
                    { key: "rank", label: "#" },
                    { key: "symbol", label: "gene" },
                    { key: "coefficient", label: `${data.method} r` },
                  ],
                  data.network.map((n, i) => ({ ...n, rank: i + 1 })),
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
          <DataTable
            columns={[
              { key: "rank", label: "#", align: "right" },
              { key: "gene", label: "Gene" },
              { key: "r", label: `${data.method} r`, align: "right" },
            ]}
            rows={data.network.map((n, i) => ({
              rank: i + 1,
              gene: <span className="font-mono">{n.symbol}</span>,
              r: n.coefficient.toFixed(3),
            }))}
          />
        </Panel>
      )}
    </>
  );
}

// Colors PCA points by whichever group column the caller picked, if any
// -- ETP-TF5 on GDS4299 is the canonical "does this actually separate the
// groups" sanity check, so defaulting to etp_status when present (falling
// back to the dataset's first group column otherwise) means the first
// thing a user sees after hitting Run is already colored meaningfully,
// not an undifferentiated cloud of same-colored dots.
function PCATab({ datasetId }: { datasetId: string }) {
  const chartRef = useRef<SVGSVGElement>(null);
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const activeDataset = datasetList?.datasets.find((d) => d.dataset_id === datasetId);

  const [genes, setGenes] = useState<string[]>(PCA_PRESETS["ETP-TF5"]);
  const [activePreset, setActivePreset] = useState<string | null>("ETP-TF5");
  const [colorBy, setColorBy] = useState<string>("");
  useEffect(() => {
    if (!activeDataset) return;
    const preferred = activeDataset.group_columns.includes("etp_status") ? "etp_status" : activeDataset.group_columns[0];
    setColorBy((prev) => (activeDataset.group_columns.includes(prev) ? prev : preferred ?? ""));
  }, [activeDataset]);

  const [query, setQuery] = useState<{ datasetId: string; genes: string[] } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["pca", query],
    queryFn: () => api.pca(query!.datasetId, query!.genes),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length >= 2 && datasetId) setQuery({ datasetId, genes });
  };

  const plotPoints = useMemo(
    () =>
      (data?.points ?? []).map((p) => ({
        sample_id: p.sample_id,
        x: p.pc1,
        y: p.pc2 ?? 0,
        group: colorBy ? (p.group_columns?.[colorBy] ?? null) : null,
      })),
    [data, colorBy],
  );

  return (
    <>
      <Panel title="Query">
        <div className="mb-3 flex flex-wrap gap-1.5">
          {Object.keys(PCA_PRESETS).map((name) => (
            <PresetButton
              key={name}
              label={name}
              active={activePreset === name}
              onClick={() => {
                setGenes(PCA_PRESETS[name]);
                setActivePreset(name);
              }}
            />
          ))}
          <PresetButton label="Custom" active={activePreset === null} onClick={() => setActivePreset(null)} />
        </div>
        <GeneTagInput
          genes={genes}
          onChange={(g) => {
            setGenes(g);
            setActivePreset(null);
          }}
        />
        <div className="mt-4 flex flex-wrap items-end gap-4">
          {activeDataset && activeDataset.group_columns.length > 0 && (
            <div className="min-w-[16ch]">
              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Color by
              </label>
              <select
                value={colorBy}
                onChange={(e) => setColorBy(e.target.value)}
                className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent"
              >
                <option value="">None</option>
                {activeDataset.group_columns.map((g) => (
                  <option key={g} value={g}>
                    {g}
                  </option>
                ))}
              </select>
            </div>
          )}
          <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length < 2}>
            Run PCA
          </PrimaryButton>
        </div>
        <p className="mt-2.5 text-[12px] text-ink-mute">Needs at least 2 genes present in the dataset.</p>
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <circle cx="7" cy="8" r="1.6" fill="currentColor" />
              <circle cx="14" cy="6" r="1.6" fill="currentColor" />
              <circle cx="10" cy="13" r="1.6" fill="currentColor" />
              <circle cx="17" cy="15" r="1.6" fill="currentColor" />
              <circle cx="6" cy="17" r="1.6" fill="currentColor" />
            </svg>
          }
          title="No PCA run yet"
          description="Enter a gene set to see whether these samples separate structurally on more than one axis — e.g. does an ETP signature actually cluster ETP vs non-ETP samples."
        />
      )}

      {data && (
        <Panel
          title={`PCA — ${data.genes_used.length} gene${data.genes_used.length === 1 ? "" : "s"}`}
          action={
            <div className="flex items-center gap-2.5">
              <span className="font-mono text-[10.5px] text-ink-mute">n = {data.n_samples}</span>
              <ExportButton
                svgRef={chartRef}
                filename={`pca-${data.genes_used.join("-")}`}
                title="PCA"
                subtitle={`n = ${data.n_samples} · ${data.genes_used.length} genes`}
                statLines={[
                  `PC1 variance explained = ${(data.variance_explained[0] * 100).toFixed(1)}%`,
                  ...(data.variance_explained[1] != null
                    ? [`PC2 variance explained = ${(data.variance_explained[1] * 100).toFixed(1)}%`]
                    : []),
                ]}
              />
            </div>
          }
        >
          <ScatterPlot points={plotPoints} xLabel="PC1" yLabel="PC2" svgRef={chartRef} />
          <div className="mt-4 flex gap-8 border-t border-rule pt-3">
            <Stat label="PC1 variance" value={`${(data.variance_explained[0] * 100).toFixed(1)}%`} tone="accent" />
            {data.variance_explained[1] != null && (
              <Stat label="PC2 variance" value={`${(data.variance_explained[1] * 100).toFixed(1)}%`} tone="accent" />
            )}
          </div>
          {(data.genes_missing.length > 0 || data.genes_unrecognized.length > 0 || data.genes_zero_variance.length > 0) && (
            <p className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
              {data.genes_unrecognized.length > 0 && `Not found: ${data.genes_unrecognized.join(", ")}. `}
              {data.genes_zero_variance.length > 0 &&
                `No variance across samples, excluded: ${data.genes_zero_variance.join(", ")}. `}
            </p>
          )}
        </Panel>
      )}
    </>
  );
}

export function CorrelationPage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);

  const [tab, setTab] = useState<Tab>("correlation");

  return (
    <div className="mx-auto max-w-[1000px]">
      <PageHeader
        eyebrow="Correlation"
        title="Correlation, co-expression, and structure"
        description="How two genes track together, the top genes most correlated with a query gene, or whether a gene set separates samples structurally (PCA)."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <div className="flex items-center gap-4">
          {datasets.length > 1 && (
            <div className="min-w-[24ch]">
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
          {/* flex-wrap + no overflow-hidden here deliberately -- a fixed-clip
              row is exactly what caused the "Kendall" label to get cut off
              inside CorrelationTab's own method toggle (see METHODS.md
              7.6's QA fix). Three tab labels, one of them "Co-expression
              network", is even more likely to overflow a squeezed
              container than that was. */}
          <div className="flex flex-wrap gap-px overflow-hidden rounded-[3px] border border-rule">
            {([
              { key: "correlation" as const, label: "Correlation" },
              { key: "network" as const, label: "Co-expression network" },
              { key: "pca" as const, label: "PCA" },
            ]).map((t) => (
              <button
                key={t.key}
                type="button"
                onClick={() => setTab(t.key)}
                className={`whitespace-nowrap px-3 py-2 text-[12.5px] font-medium transition-colors ${
                  tab === t.key ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
                }`}
              >
                {t.label}
              </button>
            ))}
          </div>
        </div>

        {/* Both tabs stay mounted at all times, hidden with CSS rather than
            conditionally rendered -- a conditional render unmounts the
            inactive tab, which wipes its local useState (query inputs AND
            results) the moment you switch away, then remounts it fresh
            (defaults, no results) on switching back. Found in QA: gene
            inputs looked "preserved" only because they reset to the same
            defaults each time, not because state actually survived. */}
        {datasetId && (
          <div className={tab === "correlation" ? "contents" : "hidden"}>
            <CorrelationTab datasetId={datasetId} />
          </div>
        )}
        {datasetId && (
          <div className={tab === "network" ? "contents" : "hidden"}>
            <NetworkTab datasetId={datasetId} />
          </div>
        )}
        {datasetId && (
          <div className={tab === "pca" ? "contents" : "hidden"}>
            <PCATab datasetId={datasetId} />
          </div>
        )}
      </div>
    </div>
  );
}
