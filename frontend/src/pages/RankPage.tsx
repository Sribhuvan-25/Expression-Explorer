import { useState } from "react";
import { useContextHasGenes, useLinkedDataset, useLinkedGenes } from "../lib/workbench";
import { formatUnit } from "../lib/units";
import { useQuery } from "@tanstack/react-query";
import { api, type RankResult } from "../lib/api";
import {
  DataTable,
  EmptyState,
  ErrorState,
  GeneTagInput,
  PageHeader,
  Panel,
  PresetButton,
  PrimaryButton,
  QueryPanel,
  SummaryParts,
  Gene,
  FieldLabel,
} from "../components/ui";
import { exportRowsAsCsv } from "../lib/exportTable";

const PRESETS: Record<string, string[]> = {
  "MYCN only": ["MYCN"],
};

export function RankPage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  // Defaults to whatever dataset the registry reports first, same pattern
  // as the other pages — no dataset id hardcoded here.
  const [datasetId, setDatasetId, datasetLink] = useLinkedDataset(
    datasets.map((d) => d.dataset_id),
    "depmap",
  );

  const [activePreset, setActivePreset] = useState<string | null>(useContextHasGenes() ? null : "MYCN only");
  const [genes, setGenes, genesLink] = useLinkedGenes(["MYCN"], { onAdopt: () => setActivePreset(null) });
  const [query, setQuery] = useState<{ datasetId: string; genes: string[] } | null>(null);

  const isSingleGene = query ? query.genes.length === 1 : true;

  const { data, isFetching, error } = useQuery({
    queryKey: ["rank", query],
    // rankByGene and rankBySignature return different extra fields (gene
    // vs genes+method) on top of the shared RankResult shape -- normalized
    // to RankResult here so useQuery infers one concrete type instead of
    // a union it can't resolve (this broke `vite build`'s stricter
    // type-checking pass even though `tsc --noEmit` alone accepted it).
    queryFn: async (): Promise<RankResult> =>
      query!.genes.length === 1
        ? await api.rankByGene(query!.datasetId, query!.genes[0])
        : await api.rankBySignature(query!.datasetId, query!.genes),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) setQuery({ datasetId, genes });
  };

  // Say what the number is: one gene ranks by its expression in the
  // dataset's unit, several by AUCell score. The column read "Value".
  const queriedDataset = datasets.find((d) => d.dataset_id === query?.datasetId);
  const valueLabel = isSingleGene
    ? `${query?.genes[0] ?? ""} expression${queriedDataset?.expression_unit ? ` (${formatUnit(queriedDataset.expression_unit)})` : ""}`
    : "AUCell score (0–1)";

  const metadataKeys =
    data && data.rows.length > 0
      ? Object.keys(data.rows[0]).filter((k) => k !== "sample_id" && k !== "value")
      : [];

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Sample Ranking"
        title="Rank samples by gene or signature"
        description="Every sample sorted by expression, for picking candidates — which samples carry the highest expression of a gene or signature."
      />

      <div className="flex flex-col gap-3 px-5 py-4 max-sm:px-3">
        <QueryPanel
          done={!!data && !!query}
          runKey={query}
          dirty={!!query && !!data && (query.genes.join(",") !== genes.join(",") || query.datasetId !== datasetId)}
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
                  query.genes.length === 1 ? "by expression" : "by AUCell score",
                  datasets.find((d) => d.dataset_id === query.datasetId)?.display_name,
                ]}
              />
            )
          }
        >
          {datasets.length > 1 && (
            <div className="mb-4 max-w-[36ch]">
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

          <FieldLabel link={genesLink} noun="genes">Genes</FieldLabel>
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
            <PresetButton label="Custom" active={activePreset === null} onClick={() => setActivePreset(null)} />
          </div>

          <GeneTagInput
            genes={genes}
            onChange={(g) => {
              setGenes(g);
              setActivePreset(null);
            }}
          />

          <div className="mt-4 flex flex-wrap items-center justify-between gap-3 border-t border-rule pt-3">
            <p className="text-[12px] text-ink-mute">
              One gene ranks by raw expression; more than one ranks by signature score (AUCell).
            </p>
            <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length === 0}>
              Rank {genes.length} gene{genes.length === 1 ? "" : "s"}
            </PrimaryButton>
          </div>
        </QueryPanel>

        {error && <ErrorState message={(error as Error).message} />}

        {!query && !error && (
          <EmptyState
            icon={
              <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                <path d="M6 6h5M6 12h9M6 18h13" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
              </svg>
            }
            title="No ranking run yet"
            description="Choose a gene or gene set, then run to see every sample sorted high to low."
          />
        )}

        {data && (
          <Panel
            title={isSingleGene ? `Ranked by ${query!.genes[0]}` : `Ranked by signature (${query!.genes.join(", ")})`}
            action={
              <div className="flex items-center gap-2.5">
                <span className="font-mono text-[10.5px] text-ink-mute">n = {data.n}</span>
                <button
                  type="button"
                  onClick={() =>
                    exportRowsAsCsv(
                      `${query!.datasetId}-${isSingleGene ? query!.genes[0] : query!.genes.join("-")}-ranking`,
                      [
                        { key: "dataset_id", label: "dataset_id" },
                        { key: "rank", label: "#" },
                        { key: "sample_id", label: "sample" },
                        ...metadataKeys.map((k) => ({ key: k, label: k })),
                        { key: "value", label: valueLabel },
                      ],
                      data.rows.map((row, i) => ({ ...row, rank: i + 1, dataset_id: query!.datasetId })),
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
              </div>
            }
          >
            <DataTable
              columns={[
                { key: "rank", label: "#", align: "right" },
                { key: "sample", label: "Sample" },
                ...metadataKeys.map((k) => ({ key: k, label: k.replace(/_/g, " ") })),
                { key: "value", label: valueLabel, align: "right" },
              ]}
              rows={data.rows.map((row, i) => ({
                rank: i + 1,
                sample: <span className="font-mono">{row.sample_id}</span>,
                ...Object.fromEntries(metadataKeys.map((k) => [k, row[k] ?? "—"])),
                value: typeof row.value === "number" ? row.value.toFixed(3) : row.value,
              }))}
            />
          </Panel>
        )}
      </div>
    </div>
  );
}
