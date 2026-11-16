import { useEffect, useMemo, useRef, useState } from "react";
import { formatUnit } from "../lib/units";
import { parseGeneInput, type ParsedGene } from "../lib/geneInput";
import { useQuery } from "@tanstack/react-query";
import {
  api,
  type ContrastDirection,
  type ContrastPanelCell,
  type ContrastPanelResult,
  type ContrastPanelRow,
  type ContrastSummary,
  type DatasetSummary,
} from "../lib/api";
import { EmptyState, ErrorState, PageHeader, Panel, PresetButton, PrimaryButton, Stat } from "../components/ui";
import { ExportButton } from "../components/ExportButton";
import { exportRowsAsCsv } from "../lib/exportTable";

// ---------------------------------------------------------------- helpers
// Every formatter accepts `unknown`: one malformed API field must render
// as "—", never throw and take the page down (CLAUDE.md, Frontend).

const isNum = (v: unknown): v is number => typeof v === "number" && Number.isFinite(v);

export { formatUnit };

export function fmtP(v: unknown): string {
  if (!isNum(v)) return "—";
  return v < 0.001 ? v.toExponential(2) : v.toFixed(4);
}

export function fmtExpr(v: unknown): string {
  if (!isNum(v)) return "—";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 100) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  return v.toFixed(3);
}

function fmtFc(v: unknown): string {
  if (!isNum(v)) return "—";
  return `${v > 0 ? "+" : ""}${v.toFixed(2)}`;
}

function CsvButton({ onClick }: { onClick: () => void }) {
  return (
    <button
      type="button"
      onClick={onClick}
      title="Download as CSV"
      className="flex items-center gap-1 rounded-[3px] border border-rule-firm bg-surface px-2 py-1 font-mono text-[10.5px] text-ink-mute transition-colors hover:border-accent hover:text-accent"
    >
      <svg width="10" height="10" viewBox="0 0 10 10" fill="none">
        <path d="M5 1v6M2.5 5L5 7.5L7.5 5M1.5 8.5h7" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" />
      </svg>
      CSV
    </button>
  );
}

// ------------------------------------------------------------ gene input

export const PP2A_PRESET_TEXT = [
  "MYC family: MYC, MYCL, MYCN",
  "PP2A inhibitors: CIP2A (KIAA1524), SET, ANP32A, ARPP19, PABIR1 (FAM122A)",
  "PP2A regulators: PPME1, LCMT1",
  "PP2A subunits and isoforms: PPP2R1A, PPP2R1B, PPP2CA, PPP2CB, PPP2R2A, PPP2R2B, PPP2R2C, PPP2R2D, PPP2R5A, PPP2R5B, PPP2R5C, PPP2R5D, PPP2R5E",
].join("\n");

export { parseGeneInput };
export type { ParsedGene };

// ------------------------------------------------- across both experiments

/**
 * Exported for tests: the "across both experiments" verdict for the
 * two-contrast case, in words that say WHICH experiment reaches q < 0.05.
 *
 * The earlier four-way label ("Opposite direction, not significant in
 * both") read as "significant in neither" for CIP2A, SET and PPP2R1B,
 * which are q < 0.002 in C4-2B -- and "Same direction in both" gave
 * PPP2R2C (q < 0.05 in both) the same words as ARPP19 (in neither).
 * Caught in user-standpoint testing.
 */
export function describeAcross(cells: (ContrastPanelCell | undefined)[], tags: string[]): string {
  if (cells.length !== 2 || tags.length !== 2) return "—";
  const [a, b] = cells;
  const usable = (c: ContrastPanelCell | undefined) => !!c?.tested && isNum(c.log2_fc);
  if (!usable(a) && !usable(b)) return "—";
  if (!usable(a)) return `Tested only in ${tags[1]}`;
  if (!usable(b)) return `Tested only in ${tags[0]}`;
  const prod = a!.log2_fc! * b!.log2_fc!;
  if (prod === 0) return "—";
  const sig = [a, b].map((c) => isNum(c!.q_value) && c!.q_value! < 0.05);
  const opposite = prod < 0;
  if (sig[0] && sig[1]) return opposite ? "Reciprocal (both q<0.05)" : "Same direction (both q<0.05)";
  const dir = opposite ? "Opposite direction" : "Same direction";
  if (sig[0]) return `${dir}; q<0.05 only in ${tags[0]}`;
  if (sig[1]) return `${dir}; q<0.05 only in ${tags[1]}`;
  return `${dir}; q≥0.05 in both`;
}

function perturbation(c: ContrastSummary): { gene: string; kind: "gain" | "loss" } | null {
  const m = (c.direction_note ?? "").match(/with (\S+) (overexpressed|knocked down|knocked out|depleted|silenced)/i);
  if (!m) return null;
  return { gene: m[1], kind: /overexpressed/i.test(m[2]) ? "gain" : "loss" };
}

function acrossLegend(contrasts: ContrastSummary[]): string {
  const [p, q] = contrasts.map(perturbation);
  if (p && q && p.gene === q.gene && p.kind !== q.kind) {
    return `${p.gene} is raised in one experiment and lowered in the other, so a reciprocal pattern (opposite log2FC, both q < 0.05) is what ${p.gene}-dependent regulation would produce. It is an association within these cell lines, not proof of direct regulation.`;
  }
  return "For an overexpression paired with a knockdown of the same gene, a reciprocal pattern (opposite log2FC, both q < 0.05) is what regulation dependent on that gene would produce. It is an association within these cell lines, not proof of direct regulation.";
}

function contrastTag(c: ContrastSummary, i: number): string {
  return c.stratum_value || `#${i + 1}`;
}

// ------------------------------------------------------ fold-change chart

interface ChartRow {
  heading?: string;
  row?: ContrastPanelRow;
}

function FoldChangeChart({
  items,
  contrasts,
  svgRef,
}: {
  items: ChartRow[];
  contrasts: ContrastSummary[];
  svgRef: React.RefObject<SVGSVGElement | null>;
}) {
  const W = 640;
  const plotX0 = 150;
  const plotX1 = 590;
  const subH = 12;
  const geneH = contrasts.length * subH + 8;
  const headH = 20;
  const top = 24;

  let maxAbs = 0;
  for (const it of items)
    for (const c of it.row?.contrasts ?? []) if (c?.tested && isNum(c.log2_fc)) maxAbs = Math.max(maxAbs, Math.abs(c.log2_fc));
  const L = maxAbs <= 1 ? 1 : Math.min(4, Math.ceil(maxAbs));
  const ticks = L === 1 ? [-1, -0.5, 0, 0.5, 1] : Array.from({ length: 2 * L + 1 }, (_, i) => i - L);
  const x = (v: number) => plotX0 + ((Math.max(-L, Math.min(L, v)) + L) / (2 * L)) * (plotX1 - plotX0);

  let y = top;
  const placed = items.map((it) => {
    const y0 = y;
    y += it.heading ? headH : geneH;
    return { ...it, y0 };
  });
  const H = y + 30;

  return (
    <svg
      ref={svgRef}
      viewBox={`0 0 ${W} ${H}`}
      width={W}
      height={H}
      className="block max-w-full"
      role="img"
      aria-label="log2 fold-change per contrast"
    >
      {ticks.map((t) => (
        <g key={t}>
          <line x1={x(t)} x2={x(t)} y1={top - 4} y2={y} stroke={t === 0 ? "var(--ink-mute)" : "var(--rule)"} strokeWidth={1} />
          <text x={x(t)} y={top - 8} textAnchor="middle" fontSize={9.5} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
            {t}
          </text>
        </g>
      ))}
      {placed.map((it, i) =>
        it.heading ? (
          <text key={`h${i}`} x={4} y={it.y0 + 14} fontSize={10} fontFamily="var(--f-mono)" fontWeight={600} className="fill-[var(--ink-soft)]">
            {it.heading.toUpperCase()}
          </text>
        ) : (
          <g key={`g${i}`}>
            <text x={4} y={it.y0 + 4 + subH / 2} dominantBaseline="middle" fontSize={11} fontStyle="italic" fill="var(--ink)">
              {it.row!.symbol}
            </text>
            {contrasts.map((c, ci) => {
              const cell = it.row!.contrasts?.find((x) => x?.contrast_id === c.contrast_id);
              const cy = it.y0 + 4 + ci * subH + subH / 2;
              const fc = cell?.tested && isNum(cell.log2_fc) ? cell.log2_fc : null;
              const sig = isNum(cell?.q_value) && cell!.q_value! < 0.05;
              const capped = fc != null && Math.abs(fc) > L;
              return (
                <g key={c.contrast_id}>
                  <text x={plotX0 - 6} y={cy} textAnchor="end" dominantBaseline="middle" fontSize={8.5} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
                    {contrastTag(c, ci)}
                  </text>
                  {fc == null ? (
                    <text x={x(0) + 4} y={cy} dominantBaseline="middle" fontSize={8.5} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
                      not tested
                    </text>
                  ) : (
                    <>
                      <rect
                        x={Math.min(x(0), x(fc))}
                        y={cy - 4}
                        width={Math.max(1, Math.abs(x(fc) - x(0)))}
                        height={8}
                        style={{ fill: fc > 0 ? "var(--hot)" : "var(--cool)" }}
                        fillOpacity={sig ? 1 : 0.35}
                      />
                      <text
                        x={capped ? (fc > 0 ? x(fc) - 3 : x(fc) + 3) : fc > 0 ? x(fc) + 3 : x(fc) - 3}
                        y={cy}
                        textAnchor={capped ? (fc > 0 ? "end" : "start") : fc > 0 ? "start" : "end"}
                        dominantBaseline="middle"
                        fontSize={8.5}
                        fontFamily="var(--f-mono)"
                        style={{ fill: capped ? "var(--surface)" : "var(--ink-soft)" }}
                      >
                        {fmtFc(fc)}
                        {capped ? " ›" : ""}
                      </text>
                    </>
                  )}
                </g>
              );
            })}
          </g>
        ),
      )}
      <text x={(plotX0 + plotX1) / 2} y={H - 8} textAnchor="middle" fontSize={10} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
        LOG2 FOLD-CHANGE (TREATMENT − REFERENCE) · SOLID = q &lt; 0.05 · BARS BEYOND ±{L} TRUNCATED
      </text>
    </svg>
  );
}

// ------------------------------------------------------------- gene panel

function MethodFootnote({ c }: { c: ContrastSummary }) {
  return (
    <li className="text-[11.5px] leading-snug text-ink-mute">
      <span className="font-medium text-ink-soft">{c.label}</span> — {c.method}. n = {c.n_reference} ({c.reference}:{" "}
      <span className="font-mono">{(c.reference_samples ?? []).join(", ")}</span>) vs {c.n_treatment} ({c.treatment}:{" "}
      <span className="font-mono">{(c.treatment_samples ?? []).join(", ")}</span>).{" "}
      {isNum(c.n_genes_tested) ? c.n_genes_tested.toLocaleString() : "?"} of{" "}
      {isNum(c.n_genes_input) ? c.n_genes_input.toLocaleString() : "?"} genes tested; q = Benjamini–Hochberg across all{" "}
      {isNum(c.n_genes_tested) ? c.n_genes_tested.toLocaleString() : "?"} tested genes in this experiment. {c.direction_note}.
      {isNum(c.n_dataset_total) && (
        <span data-testid="sample-accounting">
          {" "}
          Samples: {isNum(c.n_excluded) ? c.n_dataset_total - c.n_excluded : c.n_reference + c.n_treatment} of{" "}
          {c.n_dataset_total}
          {isNum(c.n_excluded) && c.n_excluded > 0 && c.exclusion_reason
            ? ` — ${c.exclusion_reason.replace(/\.\s*$/, "")}`
            : ""}
          .
        </span>
      )}
    </li>
  );
}

function orderRows(result: ContrastPanelResult, parsed: ParsedGene[]) {
  const byGene = new Map<string, ContrastPanelRow>();
  for (const r of result.rows ?? []) if (r?.gene) byGene.set(r.gene.toUpperCase(), r);
  const hasGroups = parsed.some((p) => p.group);
  const items: (ChartRow & { parsed?: ParsedGene })[] = [];
  let lastGroup: string | null | undefined = undefined;
  for (const p of parsed) {
    const row = byGene.get(p.gene.toUpperCase());
    if (!row) continue;
    if (hasGroups && p.group !== lastGroup) {
      items.push({ heading: p.group ?? "Other" });
      lastGroup = p.group;
    }
    items.push({ row, parsed: p });
  }
  // Rows the API returned that do not match an input token (should not
  // happen, but never silently drop data).
  const shown = new Set(items.map((i) => i.row).filter(Boolean));
  for (const r of result.rows ?? []) if (r && !shown.has(r)) items.push({ row: r });
  return items;
}

function CellGroup({
  cell,
  unit,
  reasonIndex,
}: {
  cell: ContrastPanelCell | undefined;
  unit: string;
  reasonIndex: (reason: string) => number;
}) {
  const td = "whitespace-nowrap border-b border-rule px-2.5 py-1.5 align-top tabular-nums";
  if (!cell) {
    return (
      <td colSpan={4} className={`${td} text-ink-mute`}>
        not in this contrast
      </td>
    );
  }
  const ref = cell.reference;
  const trt = cell.treatment;
  const reps = (vals: { value: number | null }[] | undefined) => (vals ?? []).map((v) => fmtExpr(v?.value)).join(", ");
  const noTpm = !isNum(ref?.mean_expression) && !isNum(trt?.mean_expression);
  const means = (
    <td className={`${td} text-right`}>
      {noTpm ? (
        <span className="text-ink-mute" title="No value for this gene in the dataset's expression matrix (see dataset notes); statistics come from the counts.">
          no {unit}
        </span>
      ) : (
        <>
          <div className="font-mono text-[12px] text-ink">
            {fmtExpr(ref?.mean_expression)} → {fmtExpr(trt?.mean_expression)}
          </div>
          <div className="font-mono text-[10px] text-ink-mute">
            {reps(ref?.values)} → {reps(trt?.values)}
          </div>
        </>
      )}
      {(isNum(ref?.mean_norm_cpm) || isNum(trt?.mean_norm_cpm)) && (
        <div data-testid="norm-cpm" className="mt-0.5" title="TMM-normalised CPM: the counts the test was fitted on">
          <div className="font-mono text-[11px] text-ink-soft">
            <span className="text-[9.5px] text-ink-mute">norm CPM </span>
            {fmtExpr(ref?.mean_norm_cpm)} → {fmtExpr(trt?.mean_norm_cpm)}
          </div>
          <div className="font-mono text-[10px] text-ink-mute">
            {reps(ref?.norm_cpm)} → {reps(trt?.norm_cpm)}
          </div>
        </div>
      )}
      {/* null = cannot judge (untested / no expression value): no warning. */}
      {cell.display_direction_agrees === false && (
        <div
          data-testid="direction-disagrees"
          className="mt-1 max-w-[26ch] whitespace-normal text-left text-[10.5px] leading-snug text-warn"
          title="The displayed expression means come from a different quantification than the counts the test used. The fold-change and its statistics follow the normalised counts shown above."
        >
          ⚠ {unit} change points the other way from the count-based fold-change
        </div>
      )}
    </td>
  );
  if (!cell.tested || !isNum(cell.log2_fc)) {
    const reason = cell.not_tested_reason || "Not tested (no reason given by the API).";
    return (
      <>
        {means}
        <td colSpan={3} className={`${td} text-ink-mute`} title={reason}>
          not tested<sup className="ml-0.5 font-mono text-[9px]">[{reasonIndex(reason)}]</sup>
        </td>
      </>
    );
  }
  const fc = cell.log2_fc;
  const sig = isNum(cell.q_value) && cell.q_value < 0.05;
  return (
    <>
      {means}
      <td className={`${td} text-right font-mono font-semibold ${fc > 0 ? "text-hot" : fc < 0 ? "text-cool" : "text-ink"}`}>
        {fmtFc(fc)} {fc > 0 ? "↑" : fc < 0 ? "↓" : ""}
      </td>
      <td className={`${td} text-right font-mono text-[12px]`}>{fmtP(cell.p_value)}</td>
      <td className={`${td} text-right font-mono text-[12px]`}>
        {sig ? (
          <span className="rounded-[2px] bg-accent-soft px-1 font-semibold text-accent-ink" title="q < 0.05">
            {fmtP(cell.q_value)}
          </span>
        ) : (
          fmtP(cell.q_value)
        )}
      </td>
    </>
  );
}

export function GenePanelTab({ dataset, contrasts }: { dataset: DatasetSummary; contrasts: ContrastSummary[] }) {
  const [text, setText] = useState(PP2A_PRESET_TEXT);
  const [query, setQuery] = useState<{ datasetId: string; text: string; parsed: ParsedGene[] } | null>(null);
  const chartRef = useRef<SVGSVGElement>(null);

  const { data, isFetching, error } = useQuery({
    queryKey: ["contrast-panel", query?.datasetId, query?.parsed.map((p) => p.gene)],
    queryFn: () => api.contrastPanel(query!.datasetId, query!.parsed.map((p) => p.gene)),
    enabled: !!query,
    retry: false,
  });

  const parsedNow = useMemo(() => parseGeneInput(text), [text]);
  const run = () => {
    if (parsedNow.length > 0) setQuery({ datasetId: dataset.dataset_id, text, parsed: parsedNow });
  };
  const stale = !!data && !!query && (query.text !== text || query.datasetId !== dataset.dataset_id);

  const shownContrasts = data?.contrasts?.length ? data.contrasts : contrasts;
  const unit = formatUnit(data?.expression_unit ?? dataset.expression_unit);
  const items = useMemo(() => (data && query ? orderRows(data, query.parsed) : []), [data, query]);
  const twoWay = shownContrasts.length === 2;

  // Footnote numbering for not-tested reasons. Computed up front, not
  // while child cells render, so the footnote list below sees all of them.
  const reasons = useMemo(() => {
    const out: string[] = [];
    for (const it of items)
      for (const c of it.row?.contrasts ?? [])
        if (c && (!c.tested || !isNum(c.log2_fc))) {
          const r = c.not_tested_reason || "Not tested (no reason given by the API).";
          if (!out.includes(r)) out.push(r);
        }
    return out;
  }, [items]);
  const reasonIndex = (r: string) => reasons.indexOf(r) + 1;
  const anyDisagree = items.some((it) =>
    (it.row?.contrasts ?? []).some((c) => c?.display_direction_agrees === false),
  );

  const exportCsv = () => {
    if (!data || !query) return;
    const maxRep = Math.max(1, ...shownContrasts.map((c) => Math.max(c.n_reference ?? 0, c.n_treatment ?? 0)));
    const repCols: { key: string; label: string }[] = [];
    for (const arm of ["reference", "treatment"])
      for (let i = 1; i <= maxRep; i++)
        repCols.push(
          { key: `${arm}_sample_${i}`, label: `${arm}_sample_${i}` },
          { key: `${arm}_value_${i}`, label: `${arm}_${unit}_${i}` },
          { key: `${arm}_norm_cpm_${i}`, label: `${arm}_norm_cpm_${i}` },
        );
    const rows: Record<string, unknown>[] = [];
    let group: string | null = null;
    for (const it of items) {
      if (it.heading) {
        group = it.heading;
        continue;
      }
      const r = it.row!;
      const across = twoWay ? describeAcross(shownContrasts.map((c) => r.contrasts?.find((x) => x?.contrast_id === c.contrast_id)), shownContrasts.map(contrastTag)) : "";
      for (const c of shownContrasts) {
        const cell = r.contrasts?.find((x) => x?.contrast_id === c.contrast_id);
        const row: Record<string, unknown> = {
          dataset_id: dataset.dataset_id,
          category: group ?? "",
          query: r.gene,
          symbol: r.symbol,
          feature_id: r.feature_id,
          alias_note: r.alias_note ?? "",
          gene_name: r.gene_name ?? "",
          in_expression_matrix: r.in_expression_matrix ?? "",
          n_loci_with_name: r.n_loci_with_name ?? "",
          identity_notes: (r.identity_notes ?? []).join(" | "),
          contrast_id: c.contrast_id,
          contrast_label: c.label,
          method: c.method,
          q_scope: `BH across all ${c.n_genes_tested} tested genes in this experiment`,
          direction_note: c.direction_note,
          reference: c.reference,
          treatment: c.treatment,
          n_reference: c.n_reference,
          n_treatment: c.n_treatment,
          unit,
          mean_reference: cell?.reference?.mean_expression,
          mean_treatment: cell?.treatment?.mean_expression,
          mean_norm_cpm_reference: cell?.reference?.mean_norm_cpm,
          mean_norm_cpm_treatment: cell?.treatment?.mean_norm_cpm,
          tested: cell?.tested ?? false,
          not_tested_reason: cell && !cell.tested ? cell.not_tested_reason ?? "" : cell ? "" : "not in this contrast",
          log2_fc: cell?.log2_fc,
          ave_log_cpm: cell?.ave_log_cpm,
          t: cell?.t,
          p_value: cell?.p_value,
          q_value: cell?.q_value,
          display_direction_agrees: cell?.display_direction_agrees,
          across_both_experiments: across,
        };
        for (const arm of ["reference", "treatment"] as const) {
          // norm CPM matched by sample id, never by position, so a value
          // cannot land beside another sample's expression value.
          const cpm = new Map((cell?.[arm]?.norm_cpm ?? []).map((v) => [v?.sample_id, v?.value]));
          (cell?.[arm]?.values ?? []).forEach((v, i) => {
            row[`${arm}_sample_${i + 1}`] = v?.sample_id;
            row[`${arm}_value_${i + 1}`] = v?.value;
            row[`${arm}_norm_cpm_${i + 1}`] = cpm.get(v?.sample_id) ?? null;
          });
        }
        rows.push(row);
      }
    }
    const base = Object.keys(rows[0] ?? { dataset_id: 0 }).filter((k) => !/_(sample|value|norm_cpm)_\d+$/.test(k));
    exportRowsAsCsv(
      `${dataset.dataset_id}-contrast-panel`,
      [...base.map((k) => ({ key: k, label: k })), ...repCols],
      rows,
    );
  };

  const thSub = "border-b border-rule bg-surface-sunk px-2.5 py-1 font-mono text-[9.5px] uppercase tracking-wider text-ink-mute";
  const totalCols = 1 + shownContrasts.length * 4 + (twoWay ? 1 : 0);

  return (
    <>
      <Panel title="Genes">
        <div className="mb-3 flex flex-wrap gap-1.5">
          <PresetButton label="PP2A signalling panel" active={text === PP2A_PRESET_TEXT} onClick={() => setText(PP2A_PRESET_TEXT)} />
        </div>
        <label htmlFor="contrast-genes" className="mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute">
          Gene symbols — comma, space or newline separated; "Heading: genes" lines group the results
        </label>
        <textarea
          id="contrast-genes"
          value={text}
          onChange={(e) => setText(e.target.value)}
          rows={5}
          className="w-full rounded-[3px] border border-rule bg-ground px-3 py-2 font-mono text-[12.5px] text-ink outline-none focus:border-accent"
        />
        <div className="mt-3 flex items-center justify-between gap-3">
          <span className="text-[11.5px] text-ink-mute">
            {parsedNow.length} gene{parsedNow.length === 1 ? "" : "s"} · names in parentheses are shown as notes, not sent
          </span>
          <PrimaryButton onClick={run} loading={isFetching} disabled={parsedNow.length === 0}>
            Run panel
          </PrimaryButton>
        </div>
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}

      {!query && !error && (
        <EmptyState
          icon={
            <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
              <path d="M12 3V21M6 7H12M12 12H19M8 17H12" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
            </svg>
          }
          title="No panel run yet"
          description="Enter genes or use the preset, then run to see each gene's fold-change in every designed experiment."
        />
      )}

      {data && (
        <>
          {stale && (
            <p className="rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink">
              Genes or dataset changed since this ran — run again to refresh.
            </p>
          )}
          <div className={stale ? "flex flex-col gap-5 opacity-50" : "flex flex-col gap-5"}>
            {(data.missing ?? []).length > 0 && (
              <div role="note" className="rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink">
                <p className="font-medium">Not in the results ({data.missing.length}):</p>
                <ul className="mt-1 list-disc pl-5">
                  {data.missing.map((m, i) => (
                    <li key={i}>
                      <span className="font-mono">{m?.gene}</span> — {m?.reason}
                    </li>
                  ))}
                </ul>
              </div>
            )}

            <Panel
              title={`Fold-change by experiment — ${dataset.display_name}`}
              action={
                <ExportButton
                  svgRef={chartRef}
                  filename={`${dataset.dataset_id}-contrast-panel-log2fc`}
                  title={`log2 fold-change per contrast — ${dataset.display_name}`}
                  subtitle={shownContrasts[0]?.method ?? ""}
                  statLines={[
                    ...shownContrasts.map(
                      (c, i) =>
                        `${contrastTag(c, i)} = ${c.label} (n = ${c.n_reference} vs ${c.n_treatment}; q = BH across ${c.n_genes_tested} genes)`,
                    ),
                    `Dataset ${dataset.dataset_id}; log2FC = treatment − reference; solid bar = q < 0.05`,
                  ]}
                />
              }
            >
              <div className="overflow-x-auto">
                <FoldChangeChart items={items} contrasts={shownContrasts} svgRef={chartRef} />
              </div>
            </Panel>

            <Panel title="Per-gene statistics" action={<CsvButton onClick={exportCsv} />}>
              <div className="max-h-[70vh] overflow-auto rounded-[3px] border border-rule">
                <table className="w-full min-w-[56rem] border-collapse text-[12.5px]">
                  <thead>
                    <tr>
                      <th rowSpan={2} className={`${thSub} text-left align-bottom`}>Gene</th>
                      {shownContrasts.map((c) => (
                        <th key={c.contrast_id} colSpan={4} className="border-b border-l border-rule bg-surface-sunk px-2.5 py-1.5 text-left text-[11.5px] font-semibold text-ink">
                          {c.label}
                        </th>
                      ))}
                      {twoWay && (
                        <th rowSpan={2} className={`${thSub} border-l text-left align-bottom`}>Across both experiments</th>
                      )}
                    </tr>
                    <tr>
                      {shownContrasts.map((c) => (
                        <FragmentHeaders key={c.contrast_id} unit={unit} cls={thSub} />
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {items.map((it, i) =>
                      it.heading ? (
                        <tr key={`h${i}`} data-testid="category-heading">
                          <td colSpan={totalCols} className="border-b border-rule bg-surface-sunk px-2.5 py-1 font-mono text-[10px] font-semibold uppercase tracking-wider text-ink-soft">
                            {it.heading}
                          </td>
                        </tr>
                      ) : (
                        <tr key={`r${i}`} data-testid="gene-row">
                          <td className="min-w-[26ch] border-b border-rule px-2.5 py-1.5 align-top">
                            <span className="font-mono italic text-ink">{it.row!.symbol}</span>
                            {it.parsed?.note && <span className="ml-1 font-mono text-[10.5px] text-ink-mute">({it.parsed.note})</span>}
                            {it.row!.alias_note ? (
                              <div className="text-[10.5px] text-ink-mute">{it.row!.alias_note}</div>
                            ) : (
                              it.row!.gene &&
                              it.row!.symbol &&
                              it.row!.gene.toUpperCase() !== it.row!.symbol.toUpperCase() && (
                                <div className="text-[10.5px] text-ink-mute">queried as {it.row!.gene}</div>
                              )
                            )}
                            {/* Count-only genes and names shared by several
                                loci (e.g. PAR X/Y copies) -- the API says
                                why the row needs care; show it on the row. */}
                            {(it.row!.identity_notes ?? []).map((note, ni) => (
                              <div
                                key={ni}
                                data-testid="identity-note"
                                className="mt-0.5 max-w-[34ch] whitespace-normal text-[10.5px] leading-snug text-ink-mute"
                              >
                                {note}
                              </div>
                            ))}
                          </td>
                          {shownContrasts.map((c) => (
                            <CellGroup
                              key={c.contrast_id}
                              cell={it.row!.contrasts?.find((x) => x?.contrast_id === c.contrast_id)}
                              unit={unit}
                              reasonIndex={reasonIndex}
                            />
                          ))}
                          {twoWay && (
                            <td className="min-w-[18ch] border-b border-l border-rule px-2.5 py-1.5 align-top text-[12px] text-ink-soft">
                              {describeAcross(shownContrasts.map((c) => it.row!.contrasts?.find((x) => x?.contrast_id === c.contrast_id)), shownContrasts.map(contrastTag))}
                            </td>
                          )}
                        </tr>
                      ),
                    )}
                  </tbody>
                </table>
              </div>
              {twoWay && <p className="mt-2.5 text-[11.5px] text-ink-soft">Across both experiments: {acrossLegend(shownContrasts)}</p>}
              <p className="mt-1.5 text-[11.5px] text-ink-mute">
                log2FC, p and q are computed on TMM-normalised counts (norm CPM, shown under each {unit} mean with its
                replicates). The {unit} means (replicates below, small) are the dataset's own values, shown for
                orientation; they can come from a different quantification, so log2FC need not equal log2 of the {unit}{" "}
                ratio.
              </p>
              {anyDisagree && (
                <p data-testid="direction-disagrees-legend" className="mt-1 text-[11.5px] text-warn">
                  ⚠ marks a cell where the mean {unit} moves the opposite way from the count-based fold-change. The
                  fold-change and its statistics follow the normalised counts, which are shown.
                </p>
              )}
              {reasons.length > 0 && (
                <ol className="mt-2 flex flex-col gap-0.5 text-[11.5px] text-ink-mute">
                  {reasons.map((r, i) => (
                    <li key={i}>
                      [{i + 1}] {r}
                    </li>
                  ))}
                </ol>
              )}
              <ul className="mt-3 flex flex-col gap-1.5 border-t border-rule pt-2.5">
                {shownContrasts.map((c) => (
                  <MethodFootnote key={c.contrast_id} c={c} />
                ))}
              </ul>
            </Panel>
          </div>
        </>
      )}
    </>
  );
}

function FragmentHeaders({ unit, cls }: { unit: string; cls: string }) {
  return (
    <>
      <th className={`${cls} border-l text-right`}>Mean {unit}: ref → trt</th>
      <th className={`${cls} text-right`}>log2FC</th>
      <th className={`${cls} text-right`}>p</th>
      <th className={`${cls} text-right`}>q</th>
    </>
  );
}

// -------------------------------------------------------------- top genes

function TopGenesTab({ dataset, contrasts }: { dataset: DatasetSummary; contrasts: ContrastSummary[] }) {
  const [contrastId, setContrastId] = useState("");
  const [direction, setDirection] = useState<ContrastDirection>("both");
  const [n, setN] = useState(50);
  useEffect(() => {
    if (!contrasts.some((c) => c.contrast_id === contrastId)) setContrastId(contrasts[0]?.contrast_id ?? "");
  }, [contrasts, contrastId]);

  const { data, isFetching, error } = useQuery({
    queryKey: ["contrast-top", dataset.dataset_id, contrastId, n, direction],
    queryFn: () => api.contrastTop(dataset.dataset_id, contrastId, n, direction),
    enabled: !!contrastId,
    retry: false,
  });
  const unit = formatUnit(data?.expression_unit ?? dataset.expression_unit);
  const genes = data?.genes ?? [];
  const c = data?.contrast;
  const anyNoTpm = genes.some((g) => g?.in_expression_matrix === false);
  const sel = "w-full rounded-[3px] border border-rule bg-ground px-3 py-2 text-[13px] text-ink outline-none focus:border-accent";
  const lbl = "mb-1 block font-mono text-[10px] uppercase tracking-wider text-ink-mute";

  return (
    <>
      <Panel title="Query">
        <div className="flex flex-col gap-4 @sm:flex-row @sm:flex-wrap @sm:items-end">
          <div className="min-w-[24ch] flex-[2]">
            <label className={lbl}>Contrast</label>
            <select value={contrastId} onChange={(e) => setContrastId(e.target.value)} className={sel}>
              {contrasts.map((x) => (
                <option key={x.contrast_id} value={x.contrast_id}>
                  {x.label}
                </option>
              ))}
            </select>
          </div>
          <div className="min-w-[12ch] flex-1">
            <label className={lbl}>Direction</label>
            <select value={direction} onChange={(e) => setDirection(e.target.value as ContrastDirection)} className={sel}>
              <option value="both">Both</option>
              <option value="up">Up (log2FC &gt; 0)</option>
              <option value="down">Down (log2FC &lt; 0)</option>
            </select>
          </div>
          <div className="min-w-[10ch] flex-1">
            <label className={lbl}>Show top</label>
            <select value={n} onChange={(e) => setN(Number(e.target.value))} className={sel}>
              {[25, 50, 100, 250, 500].map((v) => (
                <option key={v} value={v}>
                  {v}
                </option>
              ))}
            </select>
          </div>
        </div>
      </Panel>

      {error && <ErrorState message={(error as Error).message} />}
      {isFetching && !data && <p className="text-[12.5px] text-ink-mute">Loading…</p>}

      {data && c && (
        <Panel
          title={`Top ${genes.length} ${data.direction === "up" ? "up-regulated " : data.direction === "down" ? "down-regulated " : ""}genes — ${c.label}`}
          action={
            <CsvButton
              onClick={() =>
                exportRowsAsCsv(
                  `${dataset.dataset_id}-${c.contrast_id}-top-${direction}`,
                  [
                    { key: "dataset_id", label: "dataset_id" },
                    { key: "contrast_id", label: "contrast_id" },
                    { key: "method", label: "method" },
                    { key: "rank", label: "#" },
                    { key: "symbol", label: "gene" },
                    { key: "feature_id", label: "feature_id" },
                    { key: "gene_name", label: "gene_name" },
                    { key: "n_loci_with_name", label: "n_loci_with_name" },
                    { key: "log2_fc", label: "log2FC (treatment - reference)" },
                    { key: "ave_log_cpm", label: "ave_log_cpm" },
                    { key: "mean_reference", label: `mean ${unit} ${c.reference}` },
                    { key: "mean_treatment", label: `mean ${unit} ${c.treatment}` },
                    { key: "mean_norm_cpm_reference", label: "mean norm CPM reference" },
                    { key: "mean_norm_cpm_treatment", label: "mean norm CPM treatment" },
                    { key: "p_value", label: "p-value" },
                    { key: "q_value", label: `q (BH across ${c.n_genes_tested} tested genes)` },
                    { key: "in_expression_matrix", label: `has ${unit}` },
                    { key: "display_direction_agrees", label: `${unit} direction agrees with log2FC (blank = cannot judge)` },
                  ],
                  genes.map((g, i) => ({ ...g, rank: i + 1, dataset_id: dataset.dataset_id, contrast_id: c.contrast_id, method: c.method })),
                )
              }
            />
          }
        >
          <div className="mb-3 flex gap-8">
            <Stat label={data.direction === "up" ? "Up genes with q < 0.05" : data.direction === "down" ? "Down genes with q < 0.05" : "Genes with q < 0.05"} value={isNum(data.n_significant_q05) ? data.n_significant_q05.toLocaleString() : "—"} tone="accent" />
            <Stat label="Genes tested" value={isNum(c.n_genes_tested) ? c.n_genes_tested.toLocaleString() : "—"} />
          </div>
          <div className="max-h-[70vh] overflow-auto rounded-[3px] border border-rule">
            <table className="w-full min-w-[40rem] border-collapse text-[12.5px]">
              <thead>
                <tr>
                  {["#", "Gene", "log2FC", `Mean ${unit}: ref → trt`, "p", "q"].map((h, i) => (
                    <th key={h} className={`sticky top-0 z-10 border-b border-rule bg-surface-sunk px-3 py-1.5 font-mono text-[10px] uppercase tracking-wider text-ink-mute ${i === 1 ? "text-left" : "text-right"}`}>
                      {h}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {genes.map((g, i) => (
                  <tr key={g?.feature_id ?? i} className="border-b border-rule last:border-none">
                    <td className="px-3 py-1.5 text-right tabular-nums">{i + 1}</td>
                    <td className="px-3 py-1.5 font-mono">
                      {g?.symbol}
                      {g?.gene_name && g.gene_name !== g.symbol && (
                        <span className="ml-1.5 text-[11px] italic text-ink-soft">{g.gene_name}</span>
                      )}
                      {g?.in_expression_matrix === false && <span className="ml-1 text-warn" title={`No ${unit} for this gene`}>*</span>}
                      {isNum(g?.n_loci_with_name) && g.n_loci_with_name > 1 && (
                        <span
                          data-testid="shared-name"
                          className="ml-1.5 whitespace-nowrap rounded-[2px] border border-warn/40 px-1 text-[10px] not-italic text-warn"
                          title={`The name ${g.gene_name ?? g.symbol} is shared by ${g.n_loci_with_name} loci (e.g. PAR X/Y copies); reads can be split between copies, so read this fold-change with care.`}
                        >
                          shared name
                        </span>
                      )}
                    </td>
                    <td className={`px-3 py-1.5 text-right font-mono font-semibold ${isNum(g?.log2_fc) && g.log2_fc > 0 ? "text-hot" : "text-cool"}`}>
                      {fmtFc(g?.log2_fc)}
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono text-[12px]">
                      {g?.in_expression_matrix === false ? <span className="text-ink-mute">no {unit}*</span> : `${fmtExpr(g?.mean_reference)} → ${fmtExpr(g?.mean_treatment)}`}
                      {(isNum(g?.mean_norm_cpm_reference) || isNum(g?.mean_norm_cpm_treatment)) && (
                        <div className="text-[10px] text-ink-mute" title="TMM-normalised CPM: the counts the test was fitted on">
                          norm CPM {fmtExpr(g?.mean_norm_cpm_reference)} → {fmtExpr(g?.mean_norm_cpm_treatment)}
                        </div>
                      )}
                      {g?.display_direction_agrees === false && (
                        <div
                          data-testid="direction-disagrees"
                          className="text-[10px] text-warn"
                          title="Statistics are computed on normalised counts; the expression means come from a different quantification."
                        >
                          ⚠ {unit} change points the other way from the count-based fold-change
                        </div>
                      )}
                    </td>
                    <td className="px-3 py-1.5 text-right font-mono text-[12px]">{fmtP(g?.p_value)}</td>
                    <td className="px-3 py-1.5 text-right font-mono text-[12px]">
                      {isNum(g?.q_value) && g.q_value < 0.05 ? (
                        <span className="rounded-[2px] bg-accent-soft px-1 font-semibold text-accent-ink">{fmtP(g.q_value)}</span>
                      ) : (
                        fmtP(g?.q_value)
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          {genes.some((g) => isNum(g?.n_loci_with_name) && g.n_loci_with_name > 1) && (
            <p className="mt-2 text-[11.5px] text-ink-mute">
              "shared name": several loci carry this gene name (e.g. PAR X/Y copies). Reads can be split between copies,
              so read the fold-change with care; for a secondary locus the Ensembl id is shown with the name beside it.
            </p>
          )}
          {genes.some((g) => g?.display_direction_agrees === false) && (
            <p data-testid="top-direction-legend" className="mt-2 text-[11.5px] text-warn">
              ⚠ marks genes whose mean {unit} moves the opposite way from the count-based fold-change. log2FC, p and q
              are computed on TMM-normalised counts (norm CPM, shown); trust them for direction.
            </p>
          )}
          {anyNoTpm && (
            <p className="mt-2 text-[11.5px] text-ink-mute">
              * No {unit} for this gene: the submitter's {unit} file is truncated, so some genes have no expression value
              here. Their statistics come from the complete counts file and are valid.
            </p>
          )}
          <ul className="mt-3 border-t border-rule pt-2.5">
            <MethodFootnote c={c} />
          </ul>
        </Panel>
      )}
    </>
  );
}

// ------------------------------------------------------------------ page

type Tab = "panel" | "top";

export function ContrastsPage() {
  const { data: datasetList, error: listError } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = useMemo(() => (datasetList?.datasets ?? []).filter((d) => (d.n_contrasts ?? 0) > 0), [datasetList]);
  const [datasetId, setDatasetId] = useState("");
  useEffect(() => {
    if (!datasets.some((d) => d.dataset_id === datasetId) && datasets.length > 0) setDatasetId(datasets[0].dataset_id);
  }, [datasets, datasetId]);
  const dataset = datasets.find((d) => d.dataset_id === datasetId);
  const [tab, setTab] = useState<Tab>("panel");

  const { data: contrastList, error: contrastError } = useQuery({
    queryKey: ["contrasts", datasetId],
    queryFn: () => api.contrasts(datasetId),
    enabled: !!datasetId,
    retry: false,
  });
  const contrasts = contrastList?.contrasts ?? [];

  return (
    <div className="mx-auto max-w-[1200px]">
      <PageHeader
        eyebrow="Perturbation Contrasts"
        title="Designed experiments, gene by gene"
        description="Precomputed count-based differential tests (limma-trend moderated t) for each within-experiment contrast — fold-change, p and genome-wide q per gene, side by side."
      />
      <div className="flex flex-col gap-5 px-8 py-6">
        {listError && <ErrorState message={(listError as Error).message} />}
        {datasetList && datasets.length === 0 && (
          <EmptyState
            icon={
              <svg width="28" height="28" viewBox="0 0 24 24" fill="none">
                <path d="M12 3V21M6 7H12M12 12H19M8 17H12" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
              </svg>
            }
            title="No dataset declares design contrasts"
            description="Perturbation contrasts are only available for datasets with designed experiments. Use Expression Compare for the others."
          />
        )}
        {dataset && (
          <Panel title="Dataset">
            <div className="max-w-[48ch]">
              <select
                aria-label="Dataset"
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
            <p className="mt-2 font-mono text-[10.5px] text-ink-mute">
              {dataset.accession}
              {dataset.n_samples != null && ` · ${dataset.n_samples} samples`}
              {dataset.assay_type && ` · ${dataset.assay_type.replace(/_/g, "-")}`}
              {` · ${formatUnit(dataset.expression_unit)}`}
              {` · ${dataset.n_contrasts} contrasts`}
            </p>
            {dataset.notes && <p className="mt-2 text-[11.5px] leading-relaxed text-ink-soft">{dataset.notes}</p>}
          </Panel>
        )}
        {contrastError && <ErrorState message={(contrastError as Error).message} />}

        {dataset && (
          <>
            <div className="flex self-start overflow-hidden rounded-[3px] border border-rule">
              {([
                { key: "panel" as const, label: "Gene panel" },
                { key: "top" as const, label: "Top genes" },
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
            {/* Both tabs stay mounted; CSS toggles them (METHODS.md 7.6). */}
            <div className={tab === "panel" ? "contents" : "hidden"}>
              <GenePanelTab dataset={dataset} contrasts={contrasts} />
            </div>
            <div className={tab === "top" ? "contents" : "hidden"}>
              <TopGenesTab dataset={dataset} contrasts={contrasts} />
            </div>
          </>
        )}
      </div>
    </div>
  );
}
