import { useEffect, useMemo, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { requestOpenPane } from "../lib/paneEvents";
import { api } from "../lib/api";
import { DataTable, EmptyState, ErrorState, FieldLabel, Gene, GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, QueryPanel, Stat, SummaryParts } from "../components/ui";
import { useContextHasGenes, useLinkedDataset, useLinkedGene, useLinkedGenes } from "../lib/workbench";
import { ScatterPlot, scatterColorKey } from "../components/ScatterPlot";
import { ExportButton } from "../components/ExportButton";
import { exportRowsAsCsv } from "../lib/exportTable";
import { formatUnit } from "../lib/units";

type Tab = "correlation" | "network" | "pca";

const PCA_PRESETS: Record<string, string[]> = {
  "ETP-TF5": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
};

// Exported figures travel without the on-page warning (METHODS.md 5.7, 11.1).
const POOLED_EXPORT_NOTE =
  "Note: pools separate designed experiments (e.g. cell lines); mostly reflects differences between them, not the perturbation.";

/** Exported for tests. Shown (and printed in the PNG) whenever the
 * points span more than one group of the colouring column. */
/**
 * The grouping worth colouring by: the first column with 2-8 groups among
 * the plotted samples. DepMap's first column, lineage, is one group (all
 * Lymphoid), so its default colouring showed nothing; cell_line_name is
 * one group per sample. With none, no colouring.
 */
export function pickColorColumn(
  columns: string[],
  points: { group_columns?: Record<string, string | null> }[],
  preferred?: string,
): string {
  const usable = (col: string) => {
    const n = new Set(points.map((p) => p.group_columns?.[col]).filter((v) => v != null)).size;
    return n >= 2 && n <= 8;
  };
  if (preferred && columns.includes(preferred) && usable(preferred)) return preferred;
  return columns.find(usable) ?? "";
}

/** "PC1 (73.2% of variance)" -- the axis alone said "PC1". */
export function pcLabel(k: number, varianceExplained: number[] | undefined): string {
  const v = varianceExplained?.[k - 1];
  return v == null ? `PC${k}` : `PC${k} (${(v * 100).toFixed(1)}% of variance)`;
}

export function poolsGroupsNote(n: number, nGroups: number, column: string): string {
  return `r is computed across all ${n} samples together, spanning ${nGroups} ${column} groups. A difference between groups alone can produce a strong r; check whether the trend also holds within each colour.`;
}

export function CorrelationTab({ datasetId, datasetName = "", pooled = false }: { datasetId: string; datasetName?: string; pooled?: boolean }) {
  const chartRef = useRef<SVGSVGElement>(null);
  const [geneA, setGeneA, geneALink] = useLinkedGene("MYCN");
  const [geneB, setGeneB] = useState("MYC");
  const [method, setMethod] = useState<"pearson" | "spearman" | "kendall">("pearson");
  const [query, setQuery] = useState<{ datasetId: string; geneA: string; geneB: string; method: typeof method } | null>(
    null,
  );
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const activeDataset = datasetList?.datasets.find((d) => d.dataset_id === datasetId);
  // r pools every sample, so a difference between groups alone can make it
  // strong (GSE66586 MYCN vs PPP2R5B: r = -0.85, mostly amplified vs not).
  // Colouring by the dataset's first grouping by default makes that
  // visible without the reader having to think to ask.
  const [colorBy, setColorBy] = useState<string>("");
  useEffect(() => {
    if (!activeDataset) return;
    setColorBy((prev) => (activeDataset.group_columns.includes(prev) ? prev : activeDataset.group_columns[0] ?? ""));
  }, [activeDataset]);

  const { data, isFetching, error } = useQuery({
    queryKey: ["correlation", query],
    queryFn: () => api.correlation(query!.datasetId, query!.geneA, query!.geneB, query!.method),
    enabled: !!query,
  });
  // Until the user picks a colouring, choose one that separates the
  // plotted samples (see pickColorColumn).
  const colorChosen = useRef(false);
  useEffect(() => {
    if (!data || !activeDataset || colorChosen.current) return;
    setColorBy(pickColorColumn(activeDataset.group_columns, data.points));
  }, [data, activeDataset]);

  const run = () => {
    const a = geneA.trim().toUpperCase();
    const b = geneB.trim().toUpperCase();
    if (a && b && datasetId) setQuery({ datasetId, geneA: a, geneB: b, method });
  };

  const plotPoints = useMemo(
    () =>
      (data?.points ?? []).map((p) => ({
        sample_id: p.sample_id,
        x: p.x,
        y: p.y,
        group: colorBy ? (p.group_columns?.[colorBy] ?? null) : null,
      })),
    [data, colorBy],
  );
  const nColorGroups = new Set(plotPoints.map((p) => p.group).filter((g) => g != null)).size;
  // The response says which scale the points are on (log2(TPM+1) for TPM
  // data, METHODS.md 7.6a); the dataset's own unit is the fallback.
  const shownUnit = data?.value_unit ?? activeDataset?.expression_unit;
  const unitSuffix = shownUnit ? ` (${formatUnit(shownUnit)})` : "";

  return (
    <>
      <QueryPanel
        done={!!data && !!query}
        runKey={query}
        dirty={!!query && !!data && (geneA.trim().toUpperCase() !== query.geneA || query.datasetId !== datasetId)}
        summary={
          query && (
            <SummaryParts
              parts={[
                <><Gene>{query.geneA}</Gene> vs <Gene>{query.geneB}</Gene></>,
                <span className="capitalize">{query.method}</span>,
                datasetName,
              ]}
            />
          )
        }
      >
        <div className="flex flex-col gap-3 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1">
            <FieldLabel link={geneALink}>Gene A</FieldLabel>
            <input
              aria-label="Gene A"
              value={geneA}
              onChange={(e) => setGeneA(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="wb-field wb-gene"
            />
          </div>
          <div className="min-w-[12ch] flex-1">
            <label className="wb-label">Gene B</label>
            <input
              aria-label="Gene B"
              value={geneB}
              onChange={(e) => setGeneB(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="wb-field wb-gene"
            />
          </div>
          <div className="w-full @sm:w-auto @sm:min-w-[27ch] @sm:flex-1">
            <label className="wb-label">Method</label>
            {/* "Kendall" is the widest label of the three -- a container that
                can be squeezed by its @sm:flex-row siblings (as the prior
                overflow-hidden version was) clips it instead of reflowing,
                since none of the three buttons can shrink past their text.
                min-w-[27ch] keeps the group from being squeezed that far in
                the first place, and dropping overflow-hidden means if it
                ever is squeezed anyway, the buttons wrap onto a second row
                instead of silently cutting off "Kendall" mid-word. */}
            <div className="flex flex-wrap overflow-hidden rounded-[3px] border border-field-edge bg-surface">
              {(["pearson", "spearman", "kendall"] as const).map((m) => (
                <button
                  aria-pressed={method === m}
                  key={m}
                  type="button"
                  onClick={() => setMethod(m)}
                  className={`flex-1 whitespace-nowrap h-[28px] px-2 font-mono text-[11.5px] capitalize transition-colors ${
                    method === m ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
                  }`}
                >
                  {m}
                </button>
              ))}
            </div>
          </div>
          {activeDataset && activeDataset.group_columns.length > 0 && (
            <div className="min-w-[16ch]">
              <label className="wb-label">
                Color by
              </label>
              <select
                aria-label="Color by"
                value={colorBy}
                onChange={(e) => {
                  colorChosen.current = true;
                  setColorBy(e.target.value);
                }}
                className="wb-field"
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
      </QueryPanel>

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
                filename={`${datasetId}-${data.gene_a}-vs-${data.gene_b}-${data.method}-correlation`}
                title={`${data.gene_a} vs ${data.gene_b}${datasetName ? ` — ${datasetName}` : ""}`}
                subtitle={`n = ${data.n} · ${data.method}`}
                statLines={[
                  `${data.method} r = ${data.coefficient.toFixed(3)}`,
                  `p = ${data.p_value < 0.001 ? data.p_value.toExponential(2) : data.p_value.toFixed(4)}`,
                  ...(colorBy && scatterColorKey(plotPoints) ? [`Colour by ${colorBy}: ${scatterColorKey(plotPoints)}`] : []),
                  ...(nColorGroups > 1 ? [poolsGroupsNote(data.n, nColorGroups, colorBy)] : []),
                  ...(pooled ? [POOLED_EXPORT_NOTE] : []),
                ]}
              />
            </div>
          }
        >
          <div className="grid gap-x-5 @4xl:grid-cols-[minmax(0,1fr)_minmax(240px,300px)]">
          <div className="min-w-0">
          <ScatterPlot
            points={plotPoints}
            xLabel={`${data.gene_a}${unitSuffix}`}
            yLabel={`${data.gene_b}${unitSuffix}`}
            svgRef={chartRef}
          />
          </div>
          <div className="min-w-0 @4xl:border-l @4xl:border-rule @4xl:pl-5">
          <div className="mt-4 flex flex-wrap gap-x-8 gap-y-3 border-t border-rule pt-3 @4xl:mt-0 @4xl:flex-col @4xl:border-t-0 @4xl:pt-0">
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
          {nColorGroups > 1 && (
            <p data-testid="pools-groups-note" className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
              {poolsGroupsNote(data.n, nColorGroups, colorBy)}
            </p>
          )}
          </div>
          </div>
        </Panel>
      )}
    </>
  );
}

export function NetworkTab({ datasetId, datasetName = "" }: { datasetId: string; datasetName?: string }) {
  const [gene, setGene, geneLink] = useLinkedGene("MYCN");
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
      <QueryPanel
        done={!!data && !!query}
        runKey={query}
        dirty={!!query && !!data && (gene.trim().toUpperCase() !== query.gene || query.datasetId !== datasetId)}
        summary={
          query && (
            <SummaryParts
              parts={[
                <Gene>{query.gene}</Gene>,
                `top ${query.topN} ${query.direction}`,
                <span className="capitalize">{query.method}</span>,
                datasetName,
              ]}
            />
          )
        }
      >
        <div className="flex flex-col gap-3 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1">
            <FieldLabel link={geneLink}>Gene</FieldLabel>
            <input
              aria-label="Gene"
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="wb-field wb-gene"
            />
          </div>
          <div className="min-w-[10ch]">
            <label className="wb-label">Top N</label>
            <input
              aria-label="Top N"
              type="number"
              min={1}
              max={200}
              value={topN}
              onChange={(e) => setTopN(Number(e.target.value))}
              className="wb-field !w-20"
            />
          </div>
          <div className="min-w-[14ch]">
            <label className="wb-label">Method</label>
            <div className="flex overflow-hidden rounded-[3px] border border-field-edge bg-surface">
              {(["pearson", "spearman"] as const).map((m) => (
                <button
                  aria-pressed={method === m}
                  key={m}
                  type="button"
                  onClick={() => setMethod(m)}
                  className={`h-[28px] px-3 font-mono text-[11.5px] capitalize transition-colors ${
                    method === m ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
                  }`}
                >
                  {m}
                </button>
              ))}
            </div>
          </div>
          <div className="min-w-[16ch]">
            <label className="wb-label">Direction</label>
            <div className="flex overflow-hidden rounded-[3px] border border-field-edge bg-surface">
              {(["positive", "negative"] as const).map((d) => (
                <button
                  aria-pressed={direction === d}
                  key={d}
                  type="button"
                  onClick={() => setDirection(d)}
                  className={`h-[28px] px-3 font-mono text-[11.5px] capitalize transition-colors ${
                    direction === d ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
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
      </QueryPanel>

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
                  // Dataset in the name and every row, as for the other CSVs
                  // (METHODS.md 5.7): one gene's network from two datasets
                  // used to overwrite one file with no provenance inside.
                  `${datasetId}-${data.gene}-${data.method}-${data.direction}-network`,
                  [
                    { key: "dataset_id", label: "dataset_id" },
                    { key: "rank", label: "#" },
                    { key: "symbol", label: "gene" },
                    { key: "coefficient", label: `${data.method} r` },
                  ],
                  data.network.map((n, i) => ({ ...n, rank: i + 1, dataset_id: datasetId })),
                )
              }
              title="Download as CSV"
              className="wb-btn-quiet"
            >
              <svg width="10" height="10" viewBox="0 0 10 10" fill="none">
                <path d="M5 1v6M2.5 5L5 7.5L7.5 5M1.5 8.5h7" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" />
              </svg>
              CSV
            </button>
          }
        >
          {data.n != null && (
            <p data-testid="network-caption" className="mb-3 text-[11.5px] text-ink-mute">
              {data.method === "pearson" ? "Pearson" : "Spearman"} r across {data.n} samples
              {data.value_unit ? `, on ${formatUnit(data.value_unit)}` : ""}. Ranked by r only, with no p-value or
              multiple-testing correction across the genes searched
              {data.n < 20 ? `; with ${data.n} samples a high r is easily reached by chance` : ""}.
            </p>
          )}
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
function PCATab({ datasetId, datasetName = "", pooled = false }: { datasetId: string; datasetName?: string; pooled?: boolean }) {
  const chartRef = useRef<SVGSVGElement>(null);
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const activeDataset = datasetList?.datasets.find((d) => d.dataset_id === datasetId);

  // PCA needs two genes: a single context gene leaves the preset alone.
  const [activePreset, setActivePreset] = useState<string | null>(useContextHasGenes(2) ? null : "ETP-TF5");
  const [genes, setGenes, genesLink] = useLinkedGenes(PCA_PRESETS["ETP-TF5"], {
    minGenes: 2,
    onAdopt: () => setActivePreset(null),
  });
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
  const colorChosen = useRef(false);
  useEffect(() => {
    if (!data || !activeDataset || colorChosen.current) return;
    setColorBy(pickColorColumn(activeDataset.group_columns, data.points, "etp_status"));
  }, [data, activeDataset]);

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
      <QueryPanel
        done={!!data && !!query}
        runKey={query}
        dirty={!!query && !!data && (genes.join(",") !== query.genes.join(",") || query.datasetId !== datasetId)}
        summary={
          query && (
            <SummaryParts
              parts={[
                <>
                  {query.genes.map((g, i) => (
                    <span key={g}>
                      {i > 0 && ", "}
                      <Gene>{g}</Gene>
                    </span>
                  ))}
                </>,
                datasetName,
              ]}
            />
          )
        }
      >
        <FieldLabel link={genesLink} noun="genes">Gene set</FieldLabel>
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
              <label className="wb-label">
                Color by
              </label>
              <select
                aria-label="Color by"
                value={colorBy}
                onChange={(e) => {
                  colorChosen.current = true;
                  setColorBy(e.target.value);
                }}
                className="wb-field"
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
      </QueryPanel>

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
                filename={`${datasetId}-pca-${data.genes_used.join("-")}`}
                title={`PCA — ${data.genes_used.slice(0, 4).join(", ")}${
                  data.genes_used.length > 4 ? ` +${data.genes_used.length - 4} more` : ""
                }${datasetName ? ` — ${datasetName}` : ""}`}
                subtitle={`n = ${data.n_samples} · ${data.genes_used.length} genes`}
                statLines={[
                  `PC1 variance explained = ${(data.variance_explained[0] * 100).toFixed(1)}%`,
                  ...(data.variance_explained[1] != null
                    ? [`PC2 variance explained = ${(data.variance_explained[1] * 100).toFixed(1)}%`]
                    : []),
                  ...(colorBy && scatterColorKey(plotPoints) ? [`Colour by ${colorBy}: ${scatterColorKey(plotPoints)}`] : []),
                  ...(pooled ? [POOLED_EXPORT_NOTE] : []),
                ]}
              />
            </div>
          }
        >
          <ScatterPlot
            points={plotPoints}
            xLabel={pcLabel(1, data.variance_explained)}
            yLabel={pcLabel(2, data.variance_explained)}
            svgRef={chartRef}
          />
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

  const [datasetId, setDatasetId, datasetLink] = useLinkedDataset(datasets.map((d) => d.dataset_id));

  const [tab, setTab] = useState<Tab>("correlation");
  const activeDataset = datasets.find((d) => d.dataset_id === datasetId);
  const isPooled = (activeDataset?.n_contrasts ?? 0) > 0;

  return (
    <div className="mx-auto max-w-[1000px]">
      <PageHeader
        eyebrow="Correlation"
        title="Correlation, co-expression, and structure"
        description="How two genes track together, the top genes most correlated with a query gene, or whether a gene set separates samples structurally (PCA)."
      />

      <div className="flex flex-col gap-3 px-5 py-4 max-sm:px-3">
        <div className="flex flex-wrap items-end gap-x-4 gap-y-3">
          {datasets.length > 1 && (
            <div className="min-w-[24ch] max-sm:w-full">
              <FieldLabel link={datasetLink} noun="dataset" display={datasets.find((d) => d.dataset_id === datasetLink.contextValue)?.display_name}>
                Dataset
              </FieldLabel>
              <select
                aria-label="Dataset"
                value={datasetId}
                onChange={(e) => setDatasetId(e.target.value)}
                className="wb-field"
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
          <div className="flex flex-wrap gap-px overflow-hidden rounded-[3px] border border-field-edge bg-surface">
            {([
              { key: "correlation" as const, label: "Correlation" },
              { key: "network" as const, label: "Co-expression network" },
              { key: "pca" as const, label: "PCA" },
            ]).map((t) => (
              <button
                aria-pressed={tab === t.key}
                key={t.key}
                type="button"
                onClick={() => setTab(t.key)}
                className={`whitespace-nowrap h-[30px] px-3 text-[12.5px] font-medium transition-colors ${
                  tab === t.key ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
                }`}
              >
                {t.label}
              </button>
            ))}
          </div>
        </div>

        {/* A dataset that declares design contrasts is several separate
            experiments (GSE324984: two cell lines, r = 0.89 across lines vs
            0.99 within). Every analysis on this page pools all its samples,
            so the result mostly tracks the experiment, not the
            perturbation -- and nothing in the correlation responses says so. */}
        {(() => {
          const active = datasets.find((d) => d.dataset_id === datasetId);
          if (!active || !((active.n_contrasts ?? 0) > 0)) return null;
          return (
            <div
              data-testid="pooled-experiments-note"
              className="rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink"
            >
              {active.display_name} combines separate designed experiments. Correlation, co-expression and PCA here run
              across all its samples together, so they mostly reflect differences between the experiments (e.g. cell
              lines), not the perturbation.{" "}
              <button
                type="button"
                onClick={() => requestOpenPane("contrasts")}
                className="font-medium text-accent underline decoration-accent/40 underline-offset-2 hover:decoration-accent"
              >
                Open Perturbation Contrasts
              </button>{" "}
              for within-experiment effects.
            </div>
          );
        })()}

        {/* Both tabs stay mounted at all times, hidden with CSS rather than
            conditionally rendered -- a conditional render unmounts the
            inactive tab, which wipes its local useState (query inputs AND
            results) the moment you switch away, then remounts it fresh
            (defaults, no results) on switching back. Found in QA: gene
            inputs looked "preserved" only because they reset to the same
            defaults each time, not because state actually survived. */}
        {datasetId && (
          <div className={tab === "correlation" ? "contents" : "hidden"}>
            <CorrelationTab datasetId={datasetId} datasetName={activeDataset?.display_name} pooled={isPooled} />
          </div>
        )}
        {datasetId && (
          <div className={tab === "network" ? "contents" : "hidden"}>
            <NetworkTab datasetId={datasetId} datasetName={activeDataset?.display_name} />
          </div>
        )}
        {datasetId && (
          <div className={tab === "pca" ? "contents" : "hidden"}>
            <PCATab datasetId={datasetId} datasetName={activeDataset?.display_name} pooled={isPooled} />
          </div>
        )}
      </div>
    </div>
  );
}
