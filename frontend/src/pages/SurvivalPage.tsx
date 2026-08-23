import { useEffect, useRef, useState } from "react";
import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";
import { GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, Stat, EmptyState, ErrorState, DataTable } from "../components/ui";
import { SurvivalPlot } from "../components/SurvivalCurve";
import { ExportButton } from "../components/ExportButton";

const PRESETS: Record<string, string[]> = {
  "ETP-TF5 (paper)": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
  "ZMIZ1-5 (paper)": ["MEF2C", "BCL2", "MYB", "MYCN", "ZMIZ1"],
};

export function SurvivalPage() {
  const chartRef = useRef<SVGSVGElement>(null);
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];
  const survivalCapable = datasets.filter((d) => d.supports_survival);

  const [datasetId, setDatasetId] = useState<string>("");
  useEffect(() => {
    if (!datasetId && survivalCapable.length > 0) setDatasetId(survivalCapable[0].dataset_id);
  }, [survivalCapable, datasetId]);

  const [genes, setGenes] = useState<string[]>(PRESETS["ETP-TF5 (paper)"]);
  const [activePreset, setActivePreset] = useState<string | null>("ETP-TF5 (paper)");
  const [query, setQuery] = useState<{ datasetId: string; genes: string[] } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["survival", query],
    queryFn: () => api.survival(query!.datasetId, query!.genes),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) setQuery({ datasetId, genes });
  };

  return (
    <div className="mx-auto max-w-[900px]">
      <PageHeader
        eyebrow="Survival"
        title="Kaplan–Meier by signature score"
        description="Samples are binarized at the median signature score into LOW/HIGH groups; log-rank test and a Cox proportional hazards model follow the paper's approach."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        {datasetList && survivalCapable.length === 0 && (
          <Panel>
            <p className="text-[13px] text-ink-soft">
              No dataset currently registered supports survival analysis (none carry clinical follow-up data).
            </p>
          </Panel>
        )}

        {survivalCapable.length > 0 && (
          <>
            <Panel title="Query">
              {survivalCapable.length > 1 && (
                <div className="mb-4 max-w-[28ch]">
                  <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                    Dataset
                  </label>
                  <select
                    value={datasetId}
                    onChange={(e) => setDatasetId(e.target.value)}
                    className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent"
                  >
                    {survivalCapable.map((d) => (
                      <option key={d.dataset_id} value={d.dataset_id}>
                        {d.display_name}
                      </option>
                    ))}
                  </select>
                </div>
              )}

              <label className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
                Gene set
              </label>
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
                  Requires ≥10 patients with both a computable score and clinical follow-up data.
                </p>
                <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length === 0}>
                  Run survival analysis
                </PrimaryButton>
              </div>
            </Panel>

            {error && <ErrorState message={(error as Error).message} />}

            {!query && !error && (
              <EmptyState
                icon={
                  <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                    <path d="M4 18C10 18 10 6 20 6" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" fill="none" />
                  </svg>
                }
                title="No survival analysis run yet"
                description="Choose a gene set and run to see Kaplan–Meier curves split at the median signature score."
              />
            )}

            {data && (
              <>
                <Panel
                  title={`Survival — ${data.genes.join(", ")}`}
                  action={
                    <div className="flex items-center gap-2.5">
                      <span className="font-mono text-[10.5px] text-ink-mute">n = {data.n}</span>
                      <ExportButton
                        svgRef={chartRef}
                        filename={`survival-${data.genes.join("-")}`}
                        title={`Survival — ${data.genes.join(", ")}`}
                        subtitle={`n = ${data.n} · LOW/HIGH split at median signature score`}
                        statLines={[
                          `Log-rank p = ${data.logrank_p_value != null ? data.logrank_p_value.toExponential(2) : "—"}`,
                          ...(data.cox
                            ? [`Cox model p = ${data.cox.log_likelihood_p_value.toExponential(2)}`]
                            : []),
                        ]}
                      />
                    </div>
                  }
                >
                  <SurvivalPlot curves={data.curves} svgRef={chartRef} />
                </Panel>

                <Panel title="Statistics">
                  <div className="flex flex-wrap gap-8">
                    <Stat
                      label="Log-rank p"
                      value={data.logrank_p_value != null ? data.logrank_p_value.toExponential(2) : "—"}
                      tone={data.logrank_p_value != null && data.logrank_p_value < 0.05 ? "hot" : "neutral"}
                    />
                    {data.cox && (
                      <Stat
                        label="Cox model p"
                        value={data.cox.log_likelihood_p_value.toExponential(2)}
                        tone={data.cox.log_likelihood_p_value < 0.05 ? "hot" : "neutral"}
                      />
                    )}
                  </div>

                  {data.cox && (
                    <div className="mt-4">
                      <DataTable
                        columns={[
                          { key: "cov", label: "Covariate" },
                          { key: "coef", label: "Coef", align: "right" },
                          { key: "hr", label: "Hazard ratio", align: "right" },
                          { key: "p", label: "p", align: "right" },
                        ]}
                        rows={Object.entries(data.cox.coefficients).map(([name, c]) => ({
                          cov: <span className="font-mono">{name}</span>,
                          coef: c.coef.toFixed(3),
                          hr: c["exp(coef)"].toFixed(3),
                          p: c.p < 0.001 ? c.p.toExponential(2) : c.p.toFixed(4),
                        }))}
                      />
                    </div>
                  )}
                </Panel>
              </>
            )}
          </>
        )}
      </div>
    </div>
  );
}
