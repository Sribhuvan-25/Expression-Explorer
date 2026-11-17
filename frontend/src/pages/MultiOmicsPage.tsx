import { useEffect, useMemo, useRef, useState } from "react";
import { formatP } from "../lib/format";
import { useQuery } from "@tanstack/react-query";
import { api, type AuxLayerSummary } from "../lib/api";
import { EmptyState, ErrorState, FieldLabel, FigureCaption, Gene, PageHeader, Panel, PrimaryButton, QueryPanel, Stat, SummaryParts } from "../components/ui";
import { useLinkedDataset, useLinkedGene } from "../lib/workbench";
import { formatUnit } from "../lib/units";
import { BoxPlot } from "../components/BoxPlot";
import { ScatterPlot } from "../components/ScatterPlot";
import { ExportButton } from "../components/ExportButton";
import { auxCaption, mutationCaption } from "../lib/figureCaptions";
import { auxAxisLabel, auxExportFilename, auxResultTitle } from "../lib/auxLabels";

type Tab = "mutation" | "crispr" | "drug";

const TAB_LABELS: Record<Tab, string> = {
  mutation: "Expression by mutation",
  crispr: "Expression vs CRISPR",
  drug: "Expression vs drug",
};

// Which aux layer each tab needs. A tab whose layer isn't present on the
// selected dataset is disabled rather than hidden -- a user who came
// looking for drug analysis should see that it exists and just isn't
// available for THIS cohort, not be left wondering where it went.
const TAB_LAYER: Record<Tab, AuxLayerSummary["layer"]> = {
  mutation: "mutation_status",
  crispr: "crispr_gene_effect",
  drug: "drug_sensitivity",
};

function CoverageNote({ layer }: { layer: AuxLayerSummary }) {
  const pct = layer.n_dataset_total > 0 ? (layer.n_samples_covered / layer.n_dataset_total) * 100 : 0;
  // Only mention the usable/advertised split when it actually matters.
  // CRISPR is 100% usable, so noting it there would be noise; PRISM is
  // ~22%, where staying silent would leave a user picking compounds that
  // mostly fail with no idea why.
  const sparse = layer.n_usable_features < layer.n_features * 0.9;
  return (
    <p className="text-[11.5px] text-ink-mute">
      {layer.value_label} — covers {layer.n_samples_covered} of {layer.n_dataset_total} samples (
      {pct.toFixed(0)}%). {layer.value_description}
      {sparse && (
        <>
          {" "}
          <span className="text-warn">
            {layer.layer === "mutation_status" ? (
              // Usable means >= 3 mutated samples for this layer, not
              // "measured" (a mutation call is never missing).
              <>
                Note: only {layer.n_usable_features.toLocaleString()} of {layer.n_features.toLocaleString()} genes
                are mutated in at least 3 of these samples; the list below offers the most-mutated ones.
              </>
            ) : (
              <>
                Note: only {layer.n_usable_features.toLocaleString()} of {layer.n_features.toLocaleString()} were
                measured on enough of these samples to analyse — most lookups will not have enough data.
              </>
            )}
          </span>
        </>
      )}
    </p>
  );
}

function MutationTab({ datasetId, layer }: { datasetId: string; layer: AuxLayerSummary }) {
  const chartRef = useRef<SVGSVGElement>(null);
  const [gene, setGene, geneLink] = useLinkedGene("MYCN");
  const [mutatedGene, setMutatedGene] = useState("");
  const [query, setQuery] = useState<{ datasetId: string; gene: string; mutatedGene: string } | null>(null);

  const { data: topMutated } = useQuery({
    queryKey: ["mutated-genes", datasetId],
    queryFn: () => api.mutatedGenes(datasetId, 30),
    retry: false,
  });

  useEffect(() => {
    if (!mutatedGene && topMutated?.genes.length) setMutatedGene(topMutated.genes[0].gene);
  }, [topMutated, mutatedGene]);

  const { data, isFetching, error } = useQuery({
    queryKey: ["compare-by-mutation", query],
    queryFn: () => api.compareByMutation(query!.datasetId, query!.gene, query!.mutatedGene),
    enabled: !!query,
  });

  const run = () => {
    const g = gene.trim().toUpperCase();
    if (g && mutatedGene) setQuery({ datasetId, gene: g, mutatedGene });
  };

  const caption = data
    ? mutationCaption({
        gene: data.gene,
        mutatedGene: data.mutated_gene,
        expressionUnit: data.expression_unit,
        nMutated: data.n_mutated,
        nWildtype: data.n_wildtype,
        pValue: data.test.p_value,
        nDatasetTotal: data.n_dataset_total,
        nExcluded: data.n_excluded,
        exclusionReason: data.exclusion_reason,
      })
    : "";

  return (
    <>
      <QueryPanel
        done={!!data && !!query}
        runKey={query}
        dirty={!!query && !!data && gene.trim().toUpperCase() !== query.gene}
        summary={
          query && (
            <SummaryParts
              parts={[
                <><Gene>{query.gene}</Gene> expression</>,
                <>split by <Gene>{query.mutatedGene}</Gene> mutation</>,
              ]}
            />
          )
        }
      >
        <div className="flex flex-col gap-3 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1 @sm:max-w-[30ch]">
            <FieldLabel link={geneLink}>Expression of</FieldLabel>
            <input
              aria-label="Expression of"
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="wb-field wb-gene"
            />
          </div>
          <div className="min-w-[18ch] flex-1 @sm:max-w-[36ch]">
            <label className="wb-label">
              Split by mutation in
            </label>
            <select
              aria-label="Split by mutation in"
              value={mutatedGene}
              onChange={(e) => setMutatedGene(e.target.value)}
              disabled={!topMutated?.genes.length}
              className="wb-field"
            >
              {topMutated?.genes.map((g) => (
                <option key={g.gene} value={g.gene}>
                  {g.gene} ({g.n_mutated} mutated)
                </option>
              ))}
            </select>
          </div>
          <PrimaryButton onClick={run} loading={isFetching} disabled={!gene.trim() || !mutatedGene}>
            Compare
          </PrimaryButton>
        </div>
        <div className="mt-3 border-t border-rule pt-2.5">
          <CoverageNote layer={layer} />
        </div>
      </QueryPanel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <path d="M4 20V10M12 20V4M20 20V14" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
            </svg>
          }
          title="No comparison run yet"
          description="Pick a gene to measure and a recurrently-mutated gene to split on — samples are grouped by whether they carry a non-silent mutation."
        />
      )}

      {data && (
        <Panel
          title={`${data.gene} by ${data.mutated_gene} mutation status`}
          action={
            <div className="flex items-center gap-2.5">
              <span className="font-mono text-[10.5px] text-ink-mute">
                n = {data.n_mutated + data.n_wildtype}
                {data.n_excluded > 0 && ` of ${data.n_dataset_total}`}
              </span>
              <ExportButton
                svgRef={chartRef}
                filename={`${datasetId}-${data.gene}-by-${data.mutated_gene}-mutation`}
                title={`${data.gene} by ${data.mutated_gene} mutation status`}
                subtitle={`${data.n_mutated} mutated vs ${data.n_wildtype} wild-type`}
                statLines={[
                  `Mann-Whitney p = ${
                    formatP(data.test.p_value)
                  }`,
                ]}
                caption={caption}
              />
            </div>
          }
        >
          {data.points.length > 0 ? (
            <>
              <BoxPlot points={data.points} valueLabel={`${data.gene} expression`} unit={data.expression_unit} svgRef={chartRef} />
              <FigureCaption text={caption} />
            </>
          ) : (
            <p className="py-8 text-center text-[12.5px] text-ink-mute">No samples to plot.</p>
          )}
          <div className="mt-4 flex flex-wrap gap-8 border-t border-rule pt-3">
            <Stat label="Mutated" value={String(data.n_mutated)} />
            <Stat label="Wild-type" value={String(data.n_wildtype)} />
            <Stat
              label="Mann-Whitney p"
              value={formatP(data.test.p_value)}
              tone={data.test.p_value != null && data.test.p_value < 0.05 ? "hot" : "neutral"}
            />
          </div>
          {data.n_excluded > 0 && (
            <p className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
              Showing {data.n_mutated + data.n_wildtype} of {data.n_dataset_total} samples —{" "}
              {data.n_excluded} excluded ({data.exclusion_reason}).
            </p>
          )}
        </Panel>
      )}
    </>
  );
}

/**
 * Searchable picker over the features an aux layer actually has data for
 * on THIS dataset.
 *
 * Replaces a free-text box that was pre-filled with a raw Broad compound
 * id. Two things made that unusable for the drug layer specifically:
 * 5,272 of PRISM's 6,790 compounds (77.6%) have zero measurements on the
 * 79 covered lines, so a typed guess failed most of the time; and the ids
 * ("BRD:BRD-K05804044-001-18-5") are not something anyone recalls, so
 * there was no way to make an informed guess in the first place. The list
 * is server-filtered to features that clear the correlation's own minimum,
 * so anything offered here will actually run.
 *
 * Free text is still accepted -- a user who has a specific id keeps the
 * fast path, and the CRISPR layer (where every feature is a gene symbol
 * and 100% are usable) loses nothing.
 */
function AuxFeaturePicker({
  datasetId,
  auxLayerKey,
  value,
  onChange,
  onSubmit,
  placeholder,
}: {
  datasetId: string;
  auxLayerKey: "crispr_gene_effect" | "drug_sensitivity";
  value: string;
  onChange: (v: string) => void;
  onSubmit: () => void;
  placeholder: string;
}) {
  const [open, setOpen] = useState(false);
  const [search, setSearch] = useState("");
  const boxRef = useRef<HTMLDivElement>(null);

  // Debounced so typing doesn't fire a request per keystroke.
  const [debounced, setDebounced] = useState("");
  useEffect(() => {
    const t = setTimeout(() => setDebounced(search), 200);
    return () => clearTimeout(t);
  }, [search]);

  useEffect(() => {
    if (!open) return;
    const onDocClick = (e: MouseEvent) => {
      if (boxRef.current && !boxRef.current.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("mousedown", onDocClick);
    return () => document.removeEventListener("mousedown", onDocClick);
  }, [open]);

  const { data, isFetching } = useQuery({
    queryKey: ["aux-features", datasetId, auxLayerKey, debounced],
    queryFn: () => api.auxFeatures(datasetId, auxLayerKey, debounced, 50),
    enabled: open,
  });

  return (
    <div ref={boxRef} className="relative">
      <div className="flex gap-1">
        <input
          aria-label="Feature"
          value={value}
          onChange={(e) => onChange(e.target.value)}
          onKeyDown={(e) => e.key === "Enter" && onSubmit()}
          placeholder={placeholder}
          className="wb-field min-w-0 font-mono"
        />
        <button
          type="button"
          onClick={() => setOpen((o) => !o)}
          title="Browse features with data in this dataset"
          className="shrink-0 rounded-[3px] border border-rule bg-surface px-2 font-mono text-[11px] text-ink-mute transition-colors hover:border-accent hover:text-accent"
        >
          Browse
        </button>
      </div>

      {open && (
        <div className="absolute z-20 mt-1 max-h-[280px] w-full overflow-y-auto rounded-[3px] border border-rule-firm bg-surface shadow-lg">
          <div className="sticky top-0 border-b border-rule bg-surface p-2">
            <input
              aria-label="Search features"
              autoFocus
              value={search}
              onChange={(e) => setSearch(e.target.value)}
              placeholder="Search name, id, or mechanism…"
              className="wb-field wb-field-sm"
            />
            {data && (
              <p className="mt-1.5 text-[10.5px] text-ink-mute">
                {data.n_matching} with data on {data.n_samples_covered} covered samples
                {data.n_matching > data.features.length && ` · showing top ${data.features.length}`}
              </p>
            )}
          </div>
          {isFetching && <p className="px-3 py-2 text-[12px] text-ink-mute">Searching…</p>}
          {data && data.features.length === 0 && !isFetching && (
            <p className="px-3 py-2 text-[12px] text-ink-mute">
              Nothing with usable data matches “{search}”.
            </p>
          )}
          {data?.features.map((f) => (
            <button
              key={f.feature_id}
              type="button"
              onClick={() => {
                onChange(f.feature_id);
                setOpen(false);
              }}
              className="flex w-full flex-col gap-0.5 border-b border-rule px-3 py-2 text-left transition-colors last:border-b-0 hover:bg-accent-soft"
            >
              <span className="flex items-baseline justify-between gap-2">
                <span className="truncate text-[12.5px] font-medium text-ink">{f.label}</span>
                <span className="shrink-0 font-mono text-[10.5px] text-ink-mute">
                  n={f.n_samples_with_value}
                </span>
              </span>
              {f.moa && <span className="truncate text-[10.5px] text-ink-mute">{f.moa}</span>}
              {f.label !== f.feature_id && (
                <span className="truncate font-mono text-[10px] text-ink-mute">{f.feature_id}</span>
              )}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

function AuxCorrelationTab({
  datasetId,
  layer,
  auxLayerKey,
  featurePlaceholder,
  featureLabel,
  defaultGene = "MYCN",
}: {
  datasetId: string;
  layer: AuxLayerSummary;
  auxLayerKey: "crispr_gene_effect" | "drug_sensitivity";
  featurePlaceholder: string;
  featureLabel: string;
  defaultGene?: string;
}) {
  const chartRef = useRef<SVGSVGElement>(null);
  const [gene, setGene, geneLink] = useLinkedGene(defaultGene);
  const [auxFeature, setAuxFeature] = useState(featurePlaceholder);
  const [method, setMethod] = useState<"pearson" | "spearman">("pearson");
  const [query, setQuery] = useState<{
    datasetId: string;
    gene: string;
    auxFeature: string;
    method: "pearson" | "spearman";
  } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["expression-vs-aux", auxLayerKey, query],
    queryFn: () => api.expressionVsAux(query!.datasetId, query!.gene, auxLayerKey, query!.auxFeature, query!.method),
    enabled: !!query,
  });

  const run = () => {
    const g = gene.trim().toUpperCase();
    const f = auxFeature.trim();
    if (g && f) setQuery({ datasetId, gene: g, auxFeature: f, method });
  };

  const title = data ? auxResultTitle(auxLayerKey, data.gene, data.aux_feature_label) : "";
  const caption = data
    ? auxCaption({
        layer: auxLayerKey,
        gene: data.gene,
        featureLabel: data.aux_feature_label,
        featureId: data.aux_feature,
        moa: data.aux_feature_moa,
        method: data.method,
        n: data.n,
        expressionUnit: data.expression_value_unit,
        valueDescription: data.value_description,
        nDatasetTotal: data.n_dataset_total,
        nExcluded: data.n_excluded,
        exclusionReason: data.exclusion_reason,
      })
    : "";

  return (
    <>
      <QueryPanel
        done={!!data && !!query}
        runKey={query}
        dirty={!!query && !!data && gene.trim().toUpperCase() !== query.gene}
        summary={
          query && (
            <SummaryParts
              parts={[
                <><Gene>{query.gene}</Gene> expression</>,
                <>vs {featureLabel.toLowerCase()} <span className="font-mono">{query.auxFeature}</span></>,
                <span className="capitalize">{query.method}</span>,
              ]}
            />
          )
        }
      >
        <div className="flex flex-col gap-3 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[12ch] flex-1 @sm:max-w-[30ch]">
            <FieldLabel link={geneLink}>Expression of</FieldLabel>
            <input
              aria-label="Expression of"
              value={gene}
              onChange={(e) => setGene(e.target.value)}
              onKeyDown={(e) => e.key === "Enter" && run()}
              className="wb-field wb-gene"
            />
          </div>
          <div className="min-w-[16ch] flex-1 @sm:max-w-[40ch]">
            <label className="wb-label">
              {featureLabel}
            </label>
            <AuxFeaturePicker
              datasetId={datasetId}
              auxLayerKey={auxLayerKey}
              value={auxFeature}
              onChange={setAuxFeature}
              onSubmit={run}
              placeholder={featurePlaceholder}
            />
          </div>
          <div className="min-w-[14ch]">
            <label className="wb-label">Method</label>
            <div className="flex flex-wrap overflow-hidden rounded-[3px] border border-field-edge bg-surface">
              {(["pearson", "spearman"] as const).map((m) => (
                <button
                  aria-pressed={method === m}
                  key={m}
                  type="button"
                  onClick={() => setMethod(m)}
                  className={`flex-1 whitespace-nowrap h-[28px] px-3 font-mono text-[11.5px] capitalize transition-colors ${
                    method === m ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
                  }`}
                >
                  {m}
                </button>
              ))}
            </div>
          </div>
          <PrimaryButton onClick={run} loading={isFetching} disabled={!gene.trim() || !auxFeature.trim()}>
            Correlate
          </PrimaryButton>
        </div>
        <div className="mt-3 border-t border-rule pt-2.5">
          <CoverageNote layer={layer} />
        </div>
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
          title="Nothing correlated yet"
          description={`Enter a gene and ${featureLabel.toLowerCase()} to see whether expression tracks with ${layer.value_label.toLowerCase()}.`}
        />
      )}

      {data && (
        <Panel
          title={title}
          action={
            <div className="flex items-center gap-2.5">
              <span className="font-mono text-[10.5px] text-ink-mute">
                n = {data.n}
                {data.n_excluded > 0 && ` of ${data.n_dataset_total}`}
              </span>
              <ExportButton
                svgRef={chartRef}
                filename={auxExportFilename({
                  datasetId,
                  layer: auxLayerKey,
                  gene: data.gene,
                  featureLabel: data.aux_feature_label,
                  method: data.method,
                })}
                title={title}
                subtitle={`n = ${data.n} · ${data.method}`}
                statLines={[
                  `${data.method} r = ${data.coefficient.toFixed(3)}`,
                  `p = ${data.p_value < 0.001 ? data.p_value.toExponential(2) : data.p_value.toFixed(4)}`,
                ]}
                caption={caption}
              />
            </div>
          }
        >
          {data.aux_feature_moa && (
            <p className="mb-2 text-[11.5px] text-ink-mute">Mechanism: {data.aux_feature_moa}</p>
          )}
          <div className="grid gap-x-5 @4xl:grid-cols-[minmax(0,1fr)_minmax(240px,300px)]">
          <div className="min-w-0">
          <ScatterPlot
            points={data.points}
            xLabel={`${data.gene} expression${data.expression_value_unit ? ` (${formatUnit(data.expression_value_unit)})` : ""}`}
            yLabel={auxAxisLabel(auxLayerKey, data.aux_feature_label)}
            svgRef={chartRef}
          />
          <FigureCaption text={caption} />
          </div>
          <div className="min-w-0 @4xl:border-l @4xl:border-rule @4xl:pl-5">
          <div className="mt-4 flex flex-wrap gap-x-8 gap-y-3 border-t border-rule pt-3 @4xl:mt-0 @4xl:flex-col @4xl:border-t-0 @4xl:pt-0">
            <Stat
              label={`${data.method[0].toUpperCase()}${data.method.slice(1)} r`}
              value={data.coefficient.toFixed(3)}
              tone={Math.abs(data.coefficient) > 0.4 ? "accent" : "neutral"}
            />
            <Stat
              label="p-value"
              value={data.p_value < 0.001 ? data.p_value.toExponential(2) : data.p_value.toFixed(4)}
              tone={data.p_value < 0.05 ? "hot" : "neutral"}
            />
          </div>
          {/* What the measurement means is in the caption under the chart. */}
          {data.n_excluded > 0 && (
            <p className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
              Showing {data.n} of {data.n_dataset_total} samples — {data.n_excluded} excluded ({data.exclusion_reason}).
            </p>
          )}
          </div>
          </div>
        </Panel>
      )}
    </>
  );
}

export function MultiOmicsPage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  const [datasetId, setDatasetId, datasetLink] = useLinkedDataset(datasets.map((d) => d.dataset_id));

  const { data: layersResult, isLoading: layersLoading } = useQuery({
    queryKey: ["aux-layers", datasetId],
    queryFn: () => api.auxLayers(datasetId),
    enabled: !!datasetId,
    retry: false,
  });
  const layers = useMemo(
    () => new Map((layersResult?.layers ?? []).map((l) => [l.layer, l])),
    [layersResult],
  );

  const [tab, setTab] = useState<Tab>("mutation");

  // Land the user on a tab this dataset can actually answer, rather than
  // showing an "unavailable" panel by default when another tab would work.
  useEffect(() => {
    if (!layersResult) return;
    if (!layers.has(TAB_LAYER[tab])) {
      const firstAvailable = (Object.keys(TAB_LAYER) as Tab[]).find((t) => layers.has(TAB_LAYER[t]));
      if (firstAvailable) setTab(firstAvailable);
    }
  }, [layersResult, layers, tab]);

  const activeLayer = layers.get(TAB_LAYER[tab]);

  return (
    <div className="mx-auto max-w-[1000px]">
      <PageHeader
        eyebrow="Multi-omics"
        title="Expression against mutations, dependencies, and drugs"
        description="Auxiliary measurements on the same samples — somatic mutation status, CRISPR knockout dependency, and drug sensitivity — joined to the expression matrix."
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
          <div className="flex flex-wrap gap-px overflow-hidden rounded-[3px] border border-field-edge bg-surface">
            {(Object.keys(TAB_LABELS) as Tab[]).map((t) => {
              const available = layers.has(TAB_LAYER[t]);
              return (
                <button
                  aria-pressed={tab === t && available}
                  key={t}
                  type="button"
                  onClick={() => available && setTab(t)}
                  disabled={!available}
                  // A bare `title` REPLACES the accessible name rather than
                  // adding to it, so every unavailable tab announced as the
                  // identical "Not available for this dataset" and a screen
                  // reader user couldn't tell which layer was missing --
                  // on a dataset with no layers at all, that was all three
                  // buttons reading the same (caught in QA). Name the layer
                  // first, then the reason.
                  title={available ? undefined : `${TAB_LABELS[t]} — not available for this dataset`}
                  className={`whitespace-nowrap h-[30px] px-3 text-[12.5px] font-medium transition-colors ${
                    // Guard on `available`: an unavailable tab that happened
                    // to be the selected one kept the active tint, so it read
                    // as both current and disabled at once.
                    tab === t && available
                      ? "bg-accent-soft text-accent-ink"
                      : "text-ink-mute hover:text-ink-soft"
                  } disabled:cursor-not-allowed disabled:opacity-40 disabled:hover:text-ink-mute`}
                >
                  {TAB_LABELS[t]}
                </button>
              );
            })}
          </div>
        </div>

        {layersLoading && <p className="text-[12.5px] text-ink-mute">Checking available layers…</p>}

        {!layersLoading && layers.size === 0 && (
          <Panel>
            <p className="text-[13px] text-ink-soft">
              This dataset carries no auxiliary measurement layers — only expression. Try DepMap (CRISPR and
              drug sensitivity) or TARGET ALL-P2 (somatic mutations).
            </p>
          </Panel>
        )}

        {activeLayer && tab === "mutation" && (
          <MutationTab key={`${datasetId}-mutation`} datasetId={datasetId} layer={activeLayer} />
        )}
        {activeLayer && tab === "crispr" && (
          <AuxCorrelationTab
            key={`${datasetId}-crispr`}
            datasetId={datasetId}
            layer={activeLayer}
            auxLayerKey="crispr_gene_effect"
            // MYB expression vs MYB dependency: screened on all 91 covered
            // lines, and the known MYB addiction (r = -0.65, METHODS.md 8).
            // RPL13A was screened on 13, so the first thing a user saw was
            // an underpowered result.
            defaultGene="MYB"
            featurePlaceholder="MYB"
            featureLabel="Dependency on gene"
          />
        )}
        {activeLayer && tab === "drug" && (
          <AuxCorrelationTab
            key={`${datasetId}-drug`}
            datasetId={datasetId}
            layer={activeLayer}
            auxLayerKey="drug_sensitivity"
            // Defaults to the single best-covered compound of the 6,790
            // (AZ-628, measured on all 79 covered lines). The previous
            // default was measured on 72, which worked but taught nothing
            // about how sparse the layer is -- and most neighbours of it
            // fail outright. "Browse" is the real answer here; this is
            // just a landing value that always runs.
            featurePlaceholder="BRD:BRD-K05804044-001-18-5"
            featureLabel="Compound"
          />
        )}
      </div>
    </div>
  );
}
