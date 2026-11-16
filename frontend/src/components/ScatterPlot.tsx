import { useEffect, useMemo, useRef, useState } from "react";
import { niceTicks, tickDigits } from "../lib/axis";

// Ordered so the first groups are clearly distinct. --cool used to sit
// third, right after --accent: two near-identical blues (#2455A4 vs
// #2C6E9B), so GSE66586's non-amplified and non-malignant lines were
// indistinguishable in PCA. Names are used for the text colour key in
// exported PNGs, which do not include the HTML legend.
export const GROUP_COLORS = ["var(--accent)", "var(--hot)", "var(--warn)", "#8B6DBF", "#5C8A6B", "var(--cool)"];
const GROUP_COLOR_NAMES = ["blue", "red", "amber", "purple", "green", "teal"];

/** Groups in plotting order (sorted), as ScatterPlot assigns colours. */
export function scatterGroups(points: { group?: string | null }[]): string[] {
  const seen = new Set<string>();
  for (const p of points) if (p.group != null) seen.add(p.group);
  return Array.from(seen).sort();
}

// Samples with no value in the colouring column, when colouring is on.
// They used to fall back to --accent -- the first group's blue -- so on
// TARGET coloured by etp_status the 279 samples with no ETP call read as
// ETP (19 samples). A neutral grey, drawn underneath, named in the key.
export const NO_GROUP_COLOR = "#A3A3A3";

/** "blue = A (n=3), red = B (n=5), grey = no value (n=2)" -- the legend,
 * as text, for exported images. */
export function scatterColorKey(points: { group?: string | null }[]): string | null {
  const groups = scatterGroups(points);
  if (groups.length === 0) return null;
  const count = (g: string) => points.filter((p) => p.group === g).length;
  const nNone = points.filter((p) => p.group == null).length;
  return [
    ...groups.map((g, i) => `${GROUP_COLOR_NAMES[i % GROUP_COLOR_NAMES.length]} = ${g} (n=${count(g)})`),
    ...(nNone > 0 ? [`grey = no value (n=${nNone})`] : []),
  ].join(", ");
}

interface PlottablePoint {
  sample_id: string;
  x: number;
  y: number;
  // Optional group label for color-coding (e.g. PCA points colored by
  // etp_status) -- plain correlation points (gene A vs gene B) have no
  // natural grouping and are left uncolored (a single accent color).
  group?: string | null;
}

// Same responsive-width approach as SurvivalCurve.tsx: measure the real
// container and lay the chart out at that width, rather than drawing at a
// fixed size and letting viewBox scale it up (which magnifies type/stroke
// widths along with the chart).
export function ScatterPlot({
  points,
  xLabel,
  yLabel,
  svgRef,
}: {
  points: PlottablePoint[];
  xLabel: string;
  yLabel: string;
  svgRef?: React.RefObject<SVGSVGElement | null>;
}) {
  const [hover, setHover] = useState<{ x: number; y: number; sampleId: string; xv: number; yv: number } | null>(null);

  const groups = useMemo(() => scatterGroups(points), [points]);
  const colorFor = (group: string | null | undefined) =>
    group == null
      ? groups.length > 0
        ? NO_GROUP_COLOR
        : "var(--accent)"
      : GROUP_COLORS[groups.indexOf(group) % GROUP_COLORS.length];
  const nUngrouped = groups.length > 0 ? points.filter((p) => p.group == null).length : 0;
  // Ungrouped points first, so grey never sits on top of a coloured one.
  const drawOrder = nUngrouped > 0 ? [...points.filter((p) => p.group == null), ...points.filter((p) => p.group != null)] : points;

  const containerRef = useRef<HTMLDivElement>(null);
  const [containerWidth, setContainerWidth] = useState(0);
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const measure = () => setContainerWidth(el.clientWidth);
    measure();
    const observer = new ResizeObserver(measure);
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const width = Math.max(containerWidth || 640, 220); // lay out at a phone-width container rather than scaling down
  const margin = { top: 16, right: 20, bottom: 46, left: 56 };
  const plotW = Math.max(width - margin.left - margin.right, 130);
  const plotH = Math.round(Math.min(400, Math.max(220, plotW * 0.62)));
  const height = plotH + margin.top + margin.bottom;

  // Ticks on round values with the domain widened to the outer ticks.
  // The evenly divided padded range printed "6.7, 8.6, 10.6, 12.5" and,
  // on a narrow range, repeated labels at one decimal.
  const { xTicks, yTicks } = useMemo(() => {
    const xs = points.map((p) => p.x);
    const ys = points.map((p) => p.y);
    // Pad slightly so points at the extremes aren't drawn flush against
    // the axis line.
    const xPad = (Math.max(...xs) - Math.min(...xs)) * 0.06 || 1;
    const yPad = (Math.max(...ys) - Math.min(...ys)) * 0.06 || 1;
    // Non-negative data (log2(TPM+1), TPM) keeps its axis at 0: padding
    // put DepMap MYCN's axis at -2, a value no sample can have.
    const low = (vals: number[], pad: number) => (Math.min(...vals) >= 0 ? Math.max(0, Math.min(...vals) - pad) : Math.min(...vals) - pad);
    return {
      xTicks: niceTicks(low(xs, xPad), Math.max(...xs) + xPad, 5),
      yTicks: niceTicks(low(ys, yPad), Math.max(...ys) + yPad, 5),
    };
  }, [points]);
  const xMin = xTicks[0];
  const xMax = xTicks[xTicks.length - 1];
  const yMin = yTicks[0];
  const yMax = yTicks[yTicks.length - 1];
  const xDigits = tickDigits(xTicks[1] - xTicks[0]);
  const yDigits = tickDigits(yTicks[1] - yTicks[0]);

  const x = (v: number) => margin.left + ((v - xMin) / (xMax - xMin)) * plotW;
  const y = (v: number) => margin.top + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

  return (
    <div ref={containerRef} className="relative w-full">
      <svg
        ref={svgRef}
        viewBox={`0 0 ${width} ${height}`}
        width={width}
        height={height}
        className="block max-w-full"
        role="img"
        aria-label={`Scatter plot of ${xLabel} vs ${yLabel}`}
      >
        {yTicks.map((t, i) => (
          <g key={i}>
            <line x1={margin.left} x2={margin.left + plotW} y1={y(t)} y2={y(t)} stroke="var(--rule)" strokeWidth={1} />
            <text
              x={margin.left - 8}
              y={y(t)}
              textAnchor="end"
              dominantBaseline="middle"
              fontSize={10}
              fontFamily="var(--f-mono)"
              className="fill-[var(--ink-mute)]"
            >
              {t.toFixed(yDigits)}
            </text>
          </g>
        ))}
        {xTicks.map((t, i) => (
          <text
            key={i}
            x={x(t)}
            y={margin.top + plotH + 18}
            textAnchor="middle"
            fontSize={10}
            fontFamily="var(--f-mono)"
            className="fill-[var(--ink-mute)]"
          >
            {t.toFixed(xDigits)}
          </text>
        ))}
        <text
          x={margin.left + plotW / 2}
          y={height - 6}
          textAnchor="middle"
          fontSize={10.5}
          fontFamily="var(--f-mono)"
          className="fill-[var(--ink-mute)]"
          style={{ letterSpacing: "0.04em" }}
        >
          {xLabel.toUpperCase()}
        </text>
        <text
          transform={`translate(14, ${margin.top + plotH / 2}) rotate(-90)`}
          textAnchor="middle"
          fontSize={10.5}
          fontFamily="var(--f-mono)"
          className="fill-[var(--ink-mute)]"
          style={{ letterSpacing: "0.04em" }}
        >
          {yLabel.toUpperCase()}
        </text>

        <line x1={margin.left} x2={margin.left} y1={margin.top} y2={margin.top + plotH} stroke="var(--rule-firm)" />
        <line
          x1={margin.left}
          x2={margin.left + plotW}
          y1={margin.top + plotH}
          y2={margin.top + plotH}
          stroke="var(--rule-firm)"
        />

        {drawOrder.map((p) => {
          const color = colorFor(p.group);
          return (
            <circle
              key={p.sample_id}
              cx={x(p.x)}
              cy={y(p.y)}
              r={3.5}
              fill={color}
              fillOpacity={0.55}
              stroke={color}
              strokeWidth={1}
              className="cursor-pointer"
              onMouseEnter={(e) => {
                const rect = (e.target as SVGCircleElement).ownerSVGElement!.getBoundingClientRect();
                setHover({ x: e.clientX - rect.left, y: e.clientY - rect.top, sampleId: p.sample_id, xv: p.x, yv: p.y });
              }}
              onMouseLeave={() => setHover(null)}
            />
          );
        })}
      </svg>

      {(groups.length > 1 || (groups.length === 1 && nUngrouped > 0)) && (
        <div className="mt-1 flex flex-wrap items-center justify-center gap-4">
          {groups.map((g) => (
            <span key={g} className="flex items-center gap-1.5 text-[12px] text-ink-soft">
              <span className="h-2.5 w-2.5 rounded-sm" style={{ background: colorFor(g) }} />
              {g} <span className="font-mono text-[10.5px] text-ink-mute">{points.filter((p) => p.group === g).length}</span>
            </span>
          ))}
          {nUngrouped > 0 && (
            <span className="flex items-center gap-1.5 text-[12px] text-ink-soft">
              <span className="h-2.5 w-2.5 rounded-sm" style={{ background: NO_GROUP_COLOR }} />
              no value <span className="font-mono text-[10.5px] text-ink-mute">{nUngrouped}</span>
            </span>
          )}
        </div>
      )}

      {hover && (
        <div
          className="pointer-events-none absolute z-10 rounded-[3px] border border-rule bg-surface px-2.5 py-1.5 shadow-lg"
          style={{ left: hover.x + 12, top: hover.y - 10 }}
        >
          <p className="font-mono text-[10.5px] text-ink-mute">{hover.sampleId}</p>
          <p className="font-mono text-[12px] font-semibold text-ink">
            {hover.xv.toFixed(3)}, {hover.yv.toFixed(3)}
          </p>
        </div>
      )}
    </div>
  );
}
