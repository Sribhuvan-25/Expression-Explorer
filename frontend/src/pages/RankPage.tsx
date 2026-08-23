import { useEffect, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import {
  DataTable,
  EmptyState,
  ErrorState,
  GeneTagInput,
  PageHeader,
  Panel,
  PresetButton,
  PrimaryButton,
} from "../components/ui";

const PRESETS: Record<string, string[]> = {
  "MYCN only": ["MYCN"],
  "PP2A subunits (unconfirmed)": ["PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B"],
};

export function RankPage() {
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];

  // Defaults to whatever dataset the registry reports first, same pattern
  // as the other pages — no dataset id hardcoded here.
  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && datasets.length > 0) {
      const depmap = datasets.find((d) => d.dataset_id === "depmap");
      setDatasetId((depmap ?? datasets[0]).dataset_id);
    }
  }, [datasets, datasetId]);

  const [genes, setGenes] = useState<string[]>(["MYCN"]);
  const [activePreset, setActivePreset] = useState<string | null>("MYCN only");
  const [query, setQuery] = useState<{ datasetId: string; genes: string[] } | null>(null);

  const isSingleGene = query ? query.genes.length === 1 : true;

  const { data, isFetching, error } = useQuery({
    queryKey: ["rank", query],
    queryFn: () =>
      query!.genes.length === 1
        ? api.rankByGene(query!.datasetId, query!.genes[0])
        : api.rankBySignature(query!.datasetId, query!.genes),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) setQuery({ datasetId, genes });
  };

  const metadataKeys =
    data && data.rows.length > 0
      ? Object.keys(data.rows[0]).filter((k) => k !== "sample_id" && k !== "value")
      : [];

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Sample Ranking"
        title="Rank samples by gene or signature"
        description="Every sample sorted by expression, for picking candidates — e.g. which DepMap cell lines have the highest MYCN expression, for deciding what to order."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <Panel title="Query">
          {datasets.length > 1 && (
            <div className="mb-4 max-w-[28ch]">
              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Dataset
              </label>
              <select
                value={datasetId}
                onChange={(e) => setDatasetId(e.target.value)}
                className="w-full rounded-md border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent"
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
            <PresetButton label="Custom" active={activePreset === null} onClick={() => setActivePreset(null)} />
          </div>

          <GeneTagInput
            genes={genes}
            onChange={(g) => {
              setGenes(g);
              setActivePreset(null);
            }}
          />

          <div className="mt-4 flex items-center justify-between">
            <p className="text-[12px] text-ink-mute">
              One gene ranks by raw expression; more than one ranks by signature score (AUCell).
            </p>
            <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length === 0}>
              Rank {genes.length} gene{genes.length === 1 ? "" : "s"}
            </PrimaryButton>
          </div>

          <p className="mt-2.5 text-[11.5px] text-warn">
            There is no confirmed canonical "PP2A signature" — the preset above is a guess at core catalytic and
            scaffold subunits, not a validated gene set.
          </p>
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
            description="Choose a gene or gene set, then run to see every sample sorted high to low."
          />
        )}

        {data && (
          <Panel
            title={isSingleGene ? `Ranked by ${query!.genes[0]}` : `Ranked by signature (${query!.genes.join(", ")})`}
            action={<span className="font-mono text-[10.5px] text-ink-mute">n = {data.n}</span>}
          >
            <DataTable
              columns={[
                { key: "rank", label: "#", align: "right" },
                { key: "sample", label: "Sample" },
                ...metadataKeys.map((k) => ({ key: k, label: k.replace(/_/g, " ") })),
                { key: "value", label: "Value", align: "right" },
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
