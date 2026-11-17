import { useRef, useState } from "react";
import { formatP } from "../lib/format";
import { useQuery } from "@tanstack/react-query";
import { api, type SurvivalResult, type TiePolicy } from "../lib/api";
import { GeneTagInput, PageHeader, Panel, PresetButton, PrimaryButton, Stat, EmptyState, ErrorState, DataTable, QueryPanel, SummaryParts, Gene, FieldLabel, FigureCaption } from "../components/ui";
import { useContextHasGenes, useLinkedDataset, useLinkedGenes } from "../lib/workbench";
import { SurvivalPlot, survivalColorKey } from "../components/SurvivalCurve";
import { ExportButton } from "../components/ExportButton";

const PRESETS: Record<string, string[]> = {
  "ETP-TF5": ["MEF2C", "LYL1", "HHEX", "LMO2", "MYCN"],
  "ZMIZ1-5": ["MEF2C", "BCL2", "MYB", "MYCN", "ZMIZ1"],
};

export function SurvivalPage() {
  const chartRef = useRef<SVGSVGElement>(null);
  const { data: datasetList } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = datasetList?.datasets ?? [];
  const survivalCapable = datasets.filter((d) => d.supports_survival);

  const [datasetId, setDatasetId, datasetLink] = useLinkedDataset(survivalCapable.map((d) => d.dataset_id));

  // The context gene list, once set, replaces the preset (see
  // lib/workbench.tsx); until then the ETP-TF5 preset is the default.
  const [activePreset, setActivePreset] = useState<string | null>(useContextHasGenes() ? null : "ETP-TF5");
  const [genes, setGenes, genesLink] = useLinkedGenes(PRESETS["ETP-TF5"], { onAdopt: () => setActivePreset(null) });
  const [cutoffMethod, setCutoffMethod] = useState<"median" | "quartile" | "custom">("median");
  const [cutoffHighPct, setCutoffHighPct] = useState(25);
  const [cutoffLowPct, setCutoffLowPct] = useState(25);
  const [tiePolicy, setTiePolicy] = useState<TiePolicy>("trim");
  const [query, setQuery] = useState<{
    datasetId: string;
    genes: string[];
    cutoffMethod: "median" | "quartile" | "custom";
    cutoffHighPct: number;
    cutoffLowPct: number;
    tiePolicy: TiePolicy;
  } | null>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["survival", query],
    queryFn: () =>
      api.survival(query!.datasetId, query!.genes, [], {
        method: query!.cutoffMethod,
        highPct: query!.cutoffHighPct,
        lowPct: query!.cutoffLowPct,
        tiePolicy: query!.tiePolicy,
      }),
    enabled: !!query,
  });

  const run = () => {
    if (genes.length > 0 && datasetId) {
      setQuery({ datasetId, genes, cutoffMethod, cutoffHighPct, cutoffLowPct, tiePolicy });
    }
  };

  // Hedged deliberately when scores tie at a cutoff: the old wording
  // asserted "middle 50% excluded" even when the tie policy had made the
  // arms something else entirely, so the caption contradicted the arm
  // sizes printed right beside it.
  const cutoffDescription =
    query?.cutoffMethod === "quartile"
      ? data?.tie_note
        ? "top vs. bottom quartile of scores, adjusted for tied scores at the cutoff"
        : "top vs. bottom quartile of scores (middle 50% excluded)"
      : query?.cutoffMethod === "custom"
        ? `top ${query.cutoffHighPct}% vs. bottom ${query.cutoffLowPct}% of scores`
        : "median signature score";

  return (
    <div className="mx-auto max-w-[1100px]">
      <PageHeader
        eyebrow="Survival"
        title="Kaplan–Meier by signature score"
        description="Samples are split into LOW and HIGH groups by signature score (at the median, or top vs. bottom quartile or percentile), then compared with a log-rank test and a Cox proportional hazards model."
      />

      <div className="flex flex-col gap-3 px-5 py-4 max-sm:px-3">
        {datasetList && survivalCapable.length === 0 && (
          <Panel>
            <p className="text-[13px] text-ink-soft">
              No dataset currently registered supports survival analysis (none carry clinical follow-up data).
            </p>
          </Panel>
        )}

        {survivalCapable.length > 0 && (
          <>
            <QueryPanel
              done={!!data && !!query}
              runKey={query}
              dirty={!!query && !!data && (query.genes.join(",") !== genes.join(",") || query.datasetId !== datasetId)}
              summary={
                query && (
                  <SummaryParts
                    parts={[
                      <>
                        {activePreset && query.genes.join(",") === PRESETS[activePreset]?.join(",") && (
                          <span className="mr-1.5 font-mono text-[12px] font-medium text-ink">{activePreset}</span>
                        )}
                        {query.genes.map((g, i) => (
                          <span key={g}>
                            {i > 0 && ", "}
                            <Gene>{g}</Gene>
                          </span>
                        ))}
                      </>,
                      query.cutoffMethod === "custom"
                        ? `top ${query.cutoffHighPct}% vs bottom ${query.cutoffLowPct}%`
                        : `${query.cutoffMethod} split`,
                      query.cutoffMethod !== "median" && `ties: ${query.tiePolicy}`,
                      datasets.find((d) => d.dataset_id === query.datasetId)?.display_name,
                    ]}
                  />
                )
              }
            >
              {survivalCapable.length > 1 && (
                <div className="mb-4 max-w-[28ch]">
                  <FieldLabel link={datasetLink} noun="dataset" display={datasetLink.contextValue ? survivalCapable.find((d) => d.dataset_id === datasetLink.contextValue)?.display_name : undefined}>
                    Dataset
                  </FieldLabel>
                  <select
                    aria-label="Dataset"
                    value={datasetId}
                    onChange={(e) => setDatasetId(e.target.value)}
                    className="wb-field"
                  >
                    {survivalCapable.map((d) => (
                      <option key={d.dataset_id} value={d.dataset_id}>
                        {d.display_name}
                      </option>
                    ))}
                  </select>
                </div>
              )}
              {/* With a single survival-capable dataset there is no picker,
                  so without this line the page never says which cohort it
                  runs on, nor that the cell-line datasets are not offered
                  (they have no clinical follow-up). */}
              {survivalCapable.length === 1 && (
                <p data-testid="survival-dataset-note" className="mb-4 text-[12px] text-ink-mute">
                  Dataset: <span className="text-ink-soft">{survivalCapable[0].display_name}</span> — the only
                  registered dataset with clinical follow-up.
                  {datasets.length > 1 &&
                    ` Not available for ${datasets
                      .filter((d) => !d.supports_survival)
                      .map((d) => d.display_name)
                      .join(", ")} (no survival data).`}
                </p>
              )}

              <FieldLabel link={genesLink} noun="genes">
                Gene set
              </FieldLabel>
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

              <div className="mt-4">
                <label className="wb-label">
                  Group cutoff
                </label>
                <div className="flex flex-wrap items-center gap-3">
                  <div className="flex overflow-hidden rounded-[3px] border border-field-edge bg-surface">
                    {(["median", "quartile", "custom"] as const).map((m) => (
                      <button
                        aria-pressed={cutoffMethod === m}
                        key={m}
                        type="button"
                        onClick={() => setCutoffMethod(m)}
                        className={`h-[26px] px-3 font-mono text-[11.5px] capitalize transition-colors ${
                          cutoffMethod === m ? "bg-accent-soft font-semibold text-accent-ink" : "text-ink-soft hover:bg-ground hover:text-ink"
                        }`}
                      >
                        {m}
                      </button>
                    ))}
                  </div>
                  {cutoffMethod === "custom" && (
                    <div className="flex items-center gap-2 text-[12px] text-ink-soft">
                      <span>High</span>
                      <input
                        aria-label="High cutoff percent"
                        type="number"
                        min={1}
                        max={99}
                        value={cutoffHighPct}
                        onChange={(e) => setCutoffHighPct(Number(e.target.value))}
                        className="wb-field wb-field-sm !w-14"
                      />
                      <span>% · Low</span>
                      <input
                        aria-label="Low cutoff percent"
                        type="number"
                        min={1}
                        max={99}
                        value={cutoffLowPct}
                        onChange={(e) => setCutoffLowPct(Number(e.target.value))}
                        className="wb-field wb-field-sm !w-14"
                      />
                      <span>%</span>
                    </div>
                  )}
                </div>
                {cutoffMethod !== "median" && (
                  <div className="mt-3">
                    <label className="wb-label">
                      Tied scores at the cutoff
                    </label>
                    <div className="flex overflow-hidden rounded-[3px] border border-field-edge bg-surface">
                      {(
                        [
                          ["trim", "Trim to size"],
                          ["exclude", "Exclude tied"],
                          ["inclusive", "Keep all tied"],
                        ] as const
                      ).map(([value, label]) => (
                        <button
                          aria-pressed={tiePolicy === value}
                          key={value}
                          type="button"
                          onClick={() => setTiePolicy(value)}
                          className={`h-[26px] px-3 font-mono text-[11.5px] transition-colors ${
                            tiePolicy === value
                              ? "bg-accent-soft text-accent-ink"
                              : "text-ink-mute hover:text-ink-soft"
                          }`}
                        >
                          {label}
                        </button>
                      ))}
                    </div>
                    <p className="mt-1.5 text-[11.5px] text-ink-mute">
                      {tiePolicy === "trim"
                        ? "Arms keep the requested size; surplus samples sharing the exact cutoff score are dropped in a fixed order."
                        : tiePolicy === "exclude"
                          ? "Every sample sharing the exact cutoff score is dropped from both arms — never splits equal scores, but can empty an arm."
                          : "Every sample sharing the exact cutoff score joins the arm, so an arm may be far larger than the requested percentile."}
                    </p>
                  </div>
                )}
                {cutoffMethod !== "median" && (
                  <p className="mt-1.5 text-[11.5px] text-ink-mute">
                    {cutoffMethod === "quartile"
                      ? "Only the top and bottom 25% of samples by score are compared — the middle 50% are excluded, plus any adjustment for tied scores above."
                      : `Only the top ${cutoffHighPct}% and bottom ${cutoffLowPct}% of samples by score are compared — the rest are excluded.`}
                  </p>
                )}
              </div>

              <div className="mt-4 flex flex-wrap items-center justify-between gap-3 border-t border-rule pt-3">
                <p className="text-[12px] text-ink-mute">
                  Requires ≥10 patients with both a computable score and clinical follow-up data.
                </p>
                <PrimaryButton onClick={run} loading={isFetching} disabled={genes.length === 0}>
                  Run survival analysis
                </PrimaryButton>
              </div>
            </QueryPanel>

            {error && <ErrorState message={(error as Error).message} />}

            {!query && !error && (
              <EmptyState
                icon={
                  <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                    <path d="M4 18C10 18 10 6 20 6" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" fill="none" />
                  </svg>
                }
                title="No survival analysis run yet"
                description="Choose a gene set and a cutoff, then run to see Kaplan–Meier curves for the LOW and HIGH score groups."
              />
            )}

            {data && (() => {
              // A tie policy (or a percentile cutoff) can empty one arm
              // entirely -- Cox still returns a fitted coefficient in that
              // case, but it's comparing a group against nothing, which
              // shows up as a hazard ratio and CI spanning many orders of
              // magnitude (e.g. HR ~18000, CI 0.000-6.9e14) that looks like
              // a real result rather than a degenerate one. Gate the table
              // on there being at least 2 arms, the same signal the plot
              // and tie_note already use to explain the emptied arm.
              const armCount = Object.keys(data.curves).length;
              const coxIsDegenerate = armCount < 2;
              const endpointName = data.endpoint?.name ?? "Survival";
              const caption = survivalCaption(data, cutoffDescription, { coxShown: !coxIsDegenerate });
              return (
              // minmax(0,1fr) even in one column: the implicit `auto` column
              // took the Cox table's min-width (28rem), so at phone width
              // both panels -- and the KM plot laid out inside -- ran 134px
              // past the pane's right edge. The table scrolls in its own region.
              <div className="grid grid-cols-[minmax(0,1fr)] items-start gap-3 @5xl:grid-cols-[minmax(0,1fr)_minmax(440px,500px)]">
                <Panel
                  title={`${endpointName} — ${data.genes.join(", ")}`}
                  action={
                    <div className="flex items-center gap-2.5">
                      <span className="font-mono text-[10.5px] text-ink-mute">
                        n = {data.n}
                        {data.n_excluded > 0 && ` of ${data.n_dataset_total}`}
                      </span>
                      <ExportButton
                        svgRef={chartRef}
                        filename={`${datasetId}-survival-${data.genes.join("-")}-${
                          query?.cutoffMethod === "custom"
                            ? `top${query.cutoffHighPct}-bottom${query.cutoffLowPct}`
                            : query?.cutoffMethod ?? "median"
                        }${query?.cutoffMethod && query.cutoffMethod !== "median" ? `-${query.tiePolicy}` : ""}`}
                        title={`${endpointName} — ${data.genes.join(", ")}${
                          datasets.find((d) => d.dataset_id === datasetId)?.display_name
                            ? ` — ${datasets.find((d) => d.dataset_id === datasetId)!.display_name}`
                            : ""
                        }`}
                        subtitle={`n = ${data.n} · LOW/HIGH split by ${cutoffDescription}`}
                        statLines={[
                          `Curves: ${survivalColorKey(data.curves)}`,
                          `Log-rank p = ${formatP(data.logrank_p_value)}`,
                          ...(data.cox && !coxIsDegenerate
                            ? [`Cox model p = ${formatP(data.cox.log_likelihood_p_value)}`]
                            : []),
                        ]}
                        caption={caption}
                      />
                    </div>
                  }
                >
                  <p className="mb-2 text-[11.5px] text-ink-mute">Groups split by {cutoffDescription}.</p>
                  {/* A percentile cutoff landing on a block of equal scores
                      is common here (a single gene scores exactly 0 for every
                      sample below the rank threshold), and it silently made
                      the arms something other than the quartiles the caption
                      promises. Surface it next to the plot with the arm sizes
                      visible, so the split can be reconciled on sight. */}
                  {data.tie_note && (
                    <p className="mb-2.5 rounded-[3px] border border-warn/30 bg-warn-soft px-2.5 py-2 text-[11.5px] text-warn">
                      {data.tie_note.message}{" "}
                      {Object.entries(data.curves)
                        .map(([label, c]) => `${label} n=${c.n}`)
                        .join(" · ")}
                      . Change <span className="font-mono">Tied scores at the cutoff</span> to compare policies.
                      {data.tie_note.sensitivity && (
                        <span data-testid="tie-sensitivity" className="mt-1 block">
                          Which tied samples are kept matters: across {data.tie_note.sensitivity.n_draws} random
                          choices, the log-rank p ranged {fmtP(data.tie_note.sensitivity.p_low)} to{" "}
                          {fmtP(data.tie_note.sensitivity.p_high)} (middle 95%; median{" "}
                          {fmtP(data.tie_note.sensitivity.p_median)}), against {fmtP(data.logrank_p_value)} shown.
                        </span>
                      )}
                    </p>
                  )}
                  <SurvivalPlot
                    curves={data.curves}
                    svgRef={chartRef}
                    endpointName={data.endpoint?.name}
                    timeUnit={data.endpoint?.time_unit}
                  />
                  <FigureCaption text={caption} />
                  {data.n_excluded > 0 && (
                    <p className="mt-3 border-t border-rule pt-2.5 text-[11.5px] text-ink-mute">
                      Showing {data.n} of {data.n_dataset_total} samples in this dataset — {data.n_excluded} excluded
                      {data.exclusion_breakdown && data.exclusion_breakdown.length === 1 ? (
                        <> {data.exclusion_breakdown[0].reason}.</>
                      ) : data.exclusion_breakdown && data.exclusion_breakdown.length > 1 ? (
                        // One line per reason (METHODS.md 6.1): a merged count
                        // put patients without follow-up and patients cut by
                        // the quartile split behind the same number.
                        <>
                          :
                          <span data-testid="survival-exclusions" className="mt-1 flex flex-col">
                            {data.exclusion_breakdown.map((e) => (
                              <span key={e.reason}>
                                <span className="font-mono">{e.n}</span> {e.reason}
                              </span>
                            ))}
                          </span>
                        </>
                      ) : (
                        <> ({data.exclusion_reason}).</>
                      )}
                    </p>
                  )}
                </Panel>

                <Panel title="Statistics">
                  <div className="flex flex-wrap gap-x-8 gap-y-3">
                    <Stat
                      label="Log-rank p"
                      value={formatP(data.logrank_p_value)}
                      tone={data.logrank_p_value != null && data.logrank_p_value < 0.05 ? "hot" : "neutral"}
                    />
                    {data.cox && !coxIsDegenerate && (
                      <Stat
                        label="Cox model p"
                        value={formatP(data.cox.log_likelihood_p_value)}
                        tone={data.cox.log_likelihood_p_value < 0.05 ? "hot" : "neutral"}
                      />
                    )}
                  </div>

                  {data.cox && coxIsDegenerate && (
                    <p className="mt-4 rounded-[3px] border border-warn/30 bg-warn-soft px-2.5 py-2 text-[11.5px] text-warn">
                      Cox model not shown — the tied-scores policy above left only one group
                      (n={Object.values(data.curves)[0]?.n ?? 0}), so there is nothing to compare it
                      against. A hazard ratio needs two groups; choose a different{" "}
                      <span className="font-mono">Tied scores at the cutoff</span> policy to restore
                      both arms.
                    </p>
                  )}

                  {data.cox && !coxIsDegenerate && (
                    <div className="mt-4">
                      <DataTable
                        columns={[
                          { key: "cov", label: "Covariate" },
                          { key: "coef", label: "Coef", align: "right" },
                          { key: "hr", label: "Hazard ratio", align: "right" },
                          { key: "ci", label: "95% CI", align: "right" },
                          { key: "p", label: "p", align: "right" },
                        ]}
                        rows={Object.entries(data.cox.coefficients).map(([name, c]) => ({
                          cov: <span className="font-mono">{name}</span>,
                          coef: c.coef.toFixed(3),
                          hr: c["exp(coef)"].toFixed(3),
                          ci: `${c["exp(coef) lower 95%"].toFixed(3)} – ${c["exp(coef) upper 95%"].toFixed(3)}`,
                          p: c.p < 0.001 ? c.p.toExponential(2) : c.p.toFixed(4),
                        }))}
                      />
                      {"signature_score" in data.cox.coefficients && (
                        // HR 23.5 on ETP-TF5 reads as an enormous per-unit
                        // effect unless the unit is stated: the score is an
                        // AUC bounded to [0, 1], so one unit is its full range.
                        <p data-testid="cox-unit-note" className="mt-2 text-[11.5px] text-ink-mute">
                          signature_score is the continuous AUC score (0–1); its hazard ratio compares a score of 1
                          with a score of 0 — the whole range, not a typical difference between two patients. It is
                          fit on the {data.n} patients shown. The HIGH/LOW curves and log-rank p above compare groups
                          split by {cutoffDescription}.
                        </p>
                      )}
                    </div>
                  )}
                </Panel>
              </div>
              );
            })()}
          </>
        )}
      </div>
    </div>
  );
}

/**
 * The Kaplan–Meier figure's caption (METHODS.md 5.11): endpoint, how the
 * groups were split, n per arm and exclusions, the tests, and what the
 * marks and the at-risk table mean. Shown under the plot and printed into
 * the exported PNG. Every number in it is one the analysis returned.
 * Exported for tests.
 */
export function survivalCaption(
  data: SurvivalResult,
  cutoffDescription: string,
  opts: { coxShown: boolean },
): string {
  const parts: string[] = [];
  const endpoint = data.endpoint;
  parts.push(
    endpoint
      ? `Kaplan–Meier estimates of ${endpoint.name.toLowerCase()}. Event: ${endpoint.event}; ${endpoint.censoring}.`
      : "Kaplan–Meier survival estimates.",
  );
  let split = `Patients are split into HIGH and LOW groups by ${cutoffDescription}`;
  if (cutoffDescription.startsWith("median")) split += " (a score equal to the median goes to LOW)";
  if (data.tie_note) {
    split += `; ${data.tie_note.n_tied_at_cutoff} samples tied at a cutoff score were resolved by the '${data.tie_note.tie_policy}' tie policy`;
  }
  parts.push(`${split}.`);
  const arms = Object.entries(data.curves)
    .map(([label, c]) => `${label} n = ${c.n}${c.n_events != null ? ` (${c.n_events} deaths)` : ""}`)
    .join(", ");
  const accounting =
    data.n_excluded > 0
      ? `${data.n} of ${data.n_dataset_total} patients analysed; ${data.n_excluded} excluded (${
          data.exclusion_breakdown?.length
            ? data.exclusion_breakdown.map((e) => `${e.n} ${e.reason}`).join("; ")
            : data.exclusion_reason
        })`
      : `all ${data.n_dataset_total} patients analysed`;
  parts.push(`${arms}; ${accounting}.`);
  parts.push(
    data.logrank_p_value != null
      ? `Groups compared with a two-sided log-rank test, p = ${formatP(data.logrank_p_value)}.`
      : "No log-rank test: it needs two groups of at least 2 patients.",
  );
  const score = data.cox?.coefficients.signature_score;
  if (data.cox && opts.coxShown && score) {
    parts.push(
      `Cox proportional-hazards model on the continuous signature score (AUC, 0–1) in the same ${data.cox.n} patients: ` +
        `HR ${score["exp(coef)"].toFixed(3)} (95% CI ${score["exp(coef) lower 95%"].toFixed(3)}–${score["exp(coef) upper 95%"].toFixed(3)}), ` +
        `likelihood-ratio p = ${formatP(data.cox.log_likelihood_p_value)}; the HR compares a score of 1 with a score of 0 — ` +
        "the whole range of the score, not the HIGH and LOW groups.",
    );
  }
  if (Object.values(data.curves).some((c) => c.censor_marks)) {
    parts.push("Tick marks on the curves are censored patients.");
  }
  if (Object.values(data.curves).every((c) => c.at_risk)) {
    parts.push(
      "The table under the axis gives the number at risk at each time shown: patients whose follow-up is at least that long (a death or censoring at exactly that time counts as at risk).",
    );
  }
  return parts.join(" ");
}

function fmtP(p: number | null | undefined): string {
  if (p == null) return "—";
  return p < 0.001 ? p.toExponential(2) : p.toFixed(3);
}
