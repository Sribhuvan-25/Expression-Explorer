import { useEffect, useMemo, useRef, useState } from "react";
import type { ComparePoint } from "../lib/api";

const GROUP_COLORS = ["var(--cool)", "var(--hot)", "var(--accent)", "var(--warn)", "#8B6DBF", "#5C8A6B"];

// Past this many groups a vertical box plot stops working at any width:
// 35 groups need ~3000px to keep labels legible, which means the reader
// scrolls sideways and never sees the whole chart. Flipping the axes puts
// group names in a left column where they can be read in full, and lets
// the chart grow downward -- vertical scrolling is normal on a page,
// horizontal scrolling inside a pane is not.
const HORIZONTAL_THRESHOLD = 7;
// Reserved width for the per-row sample count in the horizontal layout,
// between the group name and the value axis.
const N_GUTTER = 34;

function quantile(sorted: number[], q: number): number {
  const pos = (sorted.length - 1) * q;
  const base = Math.floor(pos);
  const rest = pos - base;
  if (sorted[base + 1] !== undefined) {
    return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
  }
  return sorted[base];
}

interface GroupStats {
  group: string;
  values: number[];
  min: number;
  q1: number;
  median: number;
  q3: number;
  max: number;
  color: string;
}

export function BoxPlot({
  points,
  valueLabel = "Expression",
  svgRef,
}: {
  points: ComparePoint[];
  valueLabel?: string;
  svgRef?: React.RefObject<SVGSVGElement | null>;
}) {
  const [hover, setHover] = useState<{ x: number; y: number; sample: string; value: number } | null>(null);

  // Expression data is routinely right-skewed -- a single high-expressing
  // outlier on a linear axis compresses every other group into a sliver
  // near zero, even though the underlying stats (computed on raw values,
  // not the chart) are unaffected. log(x+1) rather than a bare log: real
  // zero values are common in this data and log(0) is undefined, but
  // log1p maps 0 -> 0 cleanly while still compressing the long tail.
  // `null` means "not chosen yet, use the automatic pick"; clicking either
  // button pins it explicitly.
  const [scaleChoice, setScaleChoice] = useState<"linear" | "log" | null>(null);
  const skewed = useMemo(() => {
    const vals = points.map((p) => p.value).filter((v) => Number.isFinite(v));
    if (vals.length < 4) return false;
    const sorted = [...vals].sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)];
    const max = sorted[sorted.length - 1];
    return max > 0 && max / Math.max(median, 1e-6) > 10;
  }, [points]);
  const logScale = scaleChoice === null ? skewed : scaleChoice === "log";
  const setLogScale = (v: boolean) => setScaleChoice(v ? "log" : "linear");
  const transform = (v: number) => (logScale ? Math.log1p(Math.max(0, v)) : v);
  const untransform = (v: number) => (logScale ? Math.expm1(v) : v);

  // The chart sizes itself to whatever space it's actually given rather
  // than to a hardcoded width.
  // Read clientWidth off the element rather than the observer's
  // contentRect: the SVG is a child of the measured div, so once the SVG
  // is sized from the measurement, contentRect reports the SVG's own width
  // and the chart latches at whatever it first rendered at. clientWidth
  // reflects the width the *parent layout* grants the div, which is the
  // number we actually want and doesn't feed back.
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

  const groups = useMemo<GroupStats[]>(() => {
    const byGroup = new Map<string, number[]>();
    for (const p of points) {
      if (!byGroup.has(p.group)) byGroup.set(p.group, []);
      byGroup.get(p.group)!.push(p.value);
    }
    // Sort by median so the chart reads as a ranking rather than an
    // arbitrary order -- with many groups this is the difference between
    // a wall of boxes and a visible trend.
    return Array.from(byGroup.entries())
      .map(([group, values]) => {
        const sorted = [...values].sort((a, b) => a - b);
        return {
          group,
          values,
          min: sorted[0],
          q1: quantile(sorted, 0.25),
          median: quantile(sorted, 0.5),
          q3: quantile(sorted, 0.75),
          max: sorted[sorted.length - 1],
          color: "",
        };
      })
      .sort((a, b) => b.median - a.median)
      .map((g, i) => ({ ...g, color: GROUP_COLORS[i % GROUP_COLORS.length] }));
  }, [points]);

  const horizontal = groups.length > HORIZONTAL_THRESHOLD;

  const availableWidth = Math.max(containerWidth || 720, 320);

  // Domain is computed in transformed space so log mode actually
  // compresses the long tail instead of relabeling a linear axis. Tick
  // values are generated evenly in transformed space then untransformed
  // for their labels, so log ticks land at genuinely log-spaced values.
  //
  // The headroom above the data is added in RAW space, before
  // transforming. Scaling the transformed max instead inflates the axis
  // exponentially once it is untransformed for the tick labels: with a
  // real max of 6.40, a 4% bump in log space printed a top tick of 7.02 --
  // a 10% overstatement, and worse for larger values. The axis would then
  // claim an expression level that does not occur in the data.
  const rawValues = points.map((p) => p.value).filter((v) => Number.isFinite(v));
  const rawMax = Math.max(...rawValues, 1) * 1.04;
  const tMax = transform(rawMax);
  const tMin = logScale ? 0 : Math.min(0, Math.min(...rawValues.map(transform)));
  const span = tMax - tMin || 1;
  const tickCount = 5;

  // Digits shown on the value axis: log mode lands on values like 0.6 or
  // 148, so a fixed 1-decimal format either loses precision at the low end
  // or adds noise at the high end.
  const fmt = (v: number) => (Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 10 ? v.toFixed(1) : v.toFixed(2));

  if (horizontal) {
    // ---- Horizontal layout: one row per group, names in a left column ----
    // Give names a real column, proportional to the longest name but
    // capped so the plot keeps most of the width. No truncation to a fixed
    // character count here -- that's what made four distinct subtypes all
    // render as "B-Lymphoblastic Leukemia/Lymp…".
    const longest = groups.reduce((m, g) => Math.max(m, g.group.length), 0);
    const labelW = Math.round(
      Math.min(Math.max(longest * 6.1 + 16 + N_GUTTER, 150), Math.max(availableWidth * 0.44, 190)),
    );
    const margin = { top: 14, right: 26, bottom: 46, left: 12 };
    const rowH = 30;
    const plotH = groups.length * rowH;
    const height = plotH + margin.top + margin.bottom;
    const width = availableWidth;
    const plotLeft = margin.left + labelW;
    const plotW = Math.max(width - plotLeft - margin.right, 120);
    const x = (v: number) => plotLeft + ((transform(v) - tMin) / span) * plotW;
    const boxH = Math.min(17, rowH * 0.6);
    const tickValues = Array.from({ length: tickCount + 1 }, (_, i) => untransform(tMin + (span * i) / tickCount));

    return (
      <div ref={containerRef} className="w-full">
        <ScaleToggle logScale={logScale} setLogScale={setLogScale} auto={scaleChoice === null && skewed} />
        <div className="relative">
          <svg
            ref={svgRef}
            viewBox={`0 0 ${width} ${height}`}
            width={width}
            height={height}
            className="block max-w-full"
            role="img"
            aria-label={`Box plot of ${valueLabel} by group`}
          >
            {tickValues.map((t, i) => (
              <g key={i}>
                <line x1={x(t)} x2={x(t)} y1={margin.top} y2={margin.top + plotH} stroke="var(--rule)" strokeWidth={1} />
                <text
                  x={x(t)}
                  y={margin.top + plotH + 16}
                  textAnchor="middle"
                  fontSize={10}
                  fontFamily="var(--f-mono)"
                  className="fill-[var(--ink-mute)]"
                >
                  {fmt(t)}
                </text>
              </g>
            ))}
            <text
              x={plotLeft + plotW / 2}
              y={height - 6}
              textAnchor="middle"
              fontSize={10.5}
              fontFamily="var(--f-mono)"
              className="fill-[var(--ink-mute)]"
              style={{ letterSpacing: "0.04em" }}
            >
              {valueLabel.toUpperCase()}
            </text>
            <line
              x1={plotLeft}
              x2={plotLeft}
              y1={margin.top}
              y2={margin.top + plotH}
              stroke="var(--rule-firm)"
              strokeWidth={1}
            />

            {groups.map((g, i) => {
              const cy = margin.top + rowH * i + rowH / 2;
              const single = g.values.length === 1;
              return (
                <g key={g.group}>
                  {i % 2 === 1 && (
                    <rect
                      x={plotLeft}
                      y={margin.top + rowH * i}
                      width={plotW}
                      height={rowH}
                      fill="var(--ink)"
                      fillOpacity={0.022}
                    />
                  )}
                  {/* n= gets its own reserved gutter between the name and
                      the axis, so it neither collides with a long group
                      name nor overlaps the boxes inside the plot. */}
                  <text
                    x={plotLeft - N_GUTTER}
                    y={cy}
                    textAnchor="end"
                    dominantBaseline="middle"
                    fontSize={11}
                    fill="var(--ink)"
                  >
                    <title>{`${g.group} (n=${g.values.length})`}</title>
                    {fitLabel(g.group, labelW - N_GUTTER - 12)}
                  </text>
                  <text
                    x={plotLeft - 8}
                    y={cy}
                    textAnchor="end"
                    dominantBaseline="middle"
                    fontSize={9}
                    fontFamily="var(--f-mono)"
                    className="fill-[var(--ink-mute)]"
                  >
                    {g.values.length}
                  </text>

                  {single ? (
                    // A single sample has no quartile range -- drawing a
                    // 1px "box" reads as a rendering bug. Show the point.
                    <circle cx={x(g.median)} cy={cy} r={4} fill={g.color} stroke="var(--surface)" strokeWidth={1} />
                  ) : (
                    <>
                      <line x1={x(g.min)} x2={x(g.q1)} y1={cy} y2={cy} stroke="var(--ink-mute)" strokeWidth={1.2} />
                      <line x1={x(g.q3)} x2={x(g.max)} y1={cy} y2={cy} stroke="var(--ink-mute)" strokeWidth={1.2} />
                      <line
                        x1={x(g.min)}
                        x2={x(g.min)}
                        y1={cy - boxH / 4}
                        y2={cy + boxH / 4}
                        stroke="var(--ink-mute)"
                        strokeWidth={1.2}
                      />
                      <line
                        x1={x(g.max)}
                        x2={x(g.max)}
                        y1={cy - boxH / 4}
                        y2={cy + boxH / 4}
                        stroke="var(--ink-mute)"
                        strokeWidth={1.2}
                      />
                      <rect
                        x={x(g.q1)}
                        y={cy - boxH / 2}
                        width={Math.max(1.5, x(g.q3) - x(g.q1))}
                        height={boxH}
                        fill={g.color}
                        fillOpacity={0.14}
                        stroke={g.color}
                        strokeWidth={1.3}
                        rx={2}
                      />
                      <line
                        x1={x(g.median)}
                        x2={x(g.median)}
                        y1={cy - boxH / 2}
                        y2={cy + boxH / 2}
                        stroke={g.color}
                        strokeWidth={2}
                      />
                    </>
                  )}

                  {points
                    .filter((p) => p.group === g.group)
                    .map((p) => (
                      <circle
                        key={p.sample_id}
                        cx={x(p.value)}
                        cy={cy + jitter(p.sample_id, boxH * 0.8)}
                        r={2.1}
                        fill={g.color}
                        fillOpacity={0.5}
                        stroke="var(--surface)"
                        strokeWidth={0.4}
                        onMouseEnter={(e) => {
                          const rect = (e.target as SVGCircleElement).ownerSVGElement!.getBoundingClientRect();
                          setHover({
                            x: e.clientX - rect.left,
                            y: e.clientY - rect.top,
                            sample: p.sample_id,
                            value: p.value,
                          });
                        }}
                        onMouseLeave={() => setHover(null)}
                        className="cursor-pointer"
                      />
                    ))}
                </g>
              );
            })}
          </svg>
          <HoverCard hover={hover} />
        </div>
      </div>
    );
  }

  // ---- Vertical layout: few enough groups that columns read well ----
  const longestLabel = groups.reduce((m, g) => Math.max(m, g.group.length), 0);
  const margin = { top: 20, right: 24, bottom: 44, left: 56 };
  const innerWidth = Math.max(availableWidth - margin.left - margin.right, 200);
  const bandW = innerWidth / Math.max(groups.length, 1);
  // With <= 7 groups there's room for horizontal labels unless the names
  // are genuinely long.
  const rotateLabels = longestLabel * 6.4 > bandW;
  const bottomMargin = rotateLabels ? Math.min(170, 54 + longestLabel * 4.4) : margin.bottom;
  const plotH = Math.round(Math.min(420, Math.max(280, innerWidth * 0.36)));
  const height = plotH + margin.top + bottomMargin;
  const width = availableWidth;
  const y = (v: number) => margin.top + plotH - ((transform(v) - tMin) / span) * plotH;
  const boxW = Math.min(120, Math.max(20, bandW * 0.42));
  const tickValues = Array.from({ length: tickCount + 1 }, (_, i) => untransform(tMin + (span * i) / tickCount));

  return (
    <div ref={containerRef} className="w-full">
      <ScaleToggle logScale={logScale} setLogScale={setLogScale} auto={scaleChoice === null && skewed} />
      <div className="relative">
        <svg
          ref={svgRef}
          viewBox={`0 0 ${width} ${height}`}
          width={width}
          height={height}
          className="block max-w-full"
          role="img"
          aria-label={`Box plot of ${valueLabel} by group`}
        >
          {tickValues.map((t, i) => (
            <g key={i}>
              <line
                x1={margin.left}
                x2={margin.left + innerWidth}
                y1={y(t)}
                y2={y(t)}
                stroke="var(--rule)"
                strokeWidth={1}
              />
              <text
                x={margin.left - 8}
                y={y(t)}
                textAnchor="end"
                dominantBaseline="middle"
                className="fill-[var(--ink-mute)]"
                fontSize={10}
                fontFamily="var(--f-mono)"
              >
                {fmt(t)}
              </text>
            </g>
          ))}
          <line
            x1={margin.left}
            x2={margin.left}
            y1={margin.top}
            y2={margin.top + plotH}
            stroke="var(--rule-firm)"
            strokeWidth={1}
          />
          <line
            x1={margin.left}
            x2={margin.left + innerWidth}
            y1={margin.top + plotH}
            y2={margin.top + plotH}
            stroke="var(--rule-firm)"
            strokeWidth={1}
          />
          <text
            transform={`translate(${14}, ${margin.top + plotH / 2}) rotate(-90)`}
            textAnchor="middle"
            fontSize={10.5}
            fontFamily="var(--f-mono)"
            className="fill-[var(--ink-mute)]"
            style={{ letterSpacing: "0.04em" }}
          >
            {valueLabel.toUpperCase()}
          </text>

          {groups.map((g, i) => {
            const cx = margin.left + bandW * i + bandW / 2;
            const single = g.values.length === 1;
            return (
              <g key={g.group}>
                {single ? (
                  <circle cx={cx} cy={y(g.median)} r={4.5} fill={g.color} stroke="var(--surface)" strokeWidth={1} />
                ) : (
                  <>
                    <line x1={cx} x2={cx} y1={y(g.max)} y2={y(g.q3)} stroke="var(--ink-mute)" strokeWidth={1.2} />
                    <line x1={cx} x2={cx} y1={y(g.min)} y2={y(g.q1)} stroke="var(--ink-mute)" strokeWidth={1.2} />
                    <line
                      x1={cx - boxW / 4}
                      x2={cx + boxW / 4}
                      y1={y(g.max)}
                      y2={y(g.max)}
                      stroke="var(--ink-mute)"
                      strokeWidth={1.2}
                    />
                    <line
                      x1={cx - boxW / 4}
                      x2={cx + boxW / 4}
                      y1={y(g.min)}
                      y2={y(g.min)}
                      stroke="var(--ink-mute)"
                      strokeWidth={1.2}
                    />
                    <rect
                      x={cx - boxW / 2}
                      y={y(g.q3)}
                      width={boxW}
                      height={Math.max(1.5, y(g.q1) - y(g.q3))}
                      fill={g.color}
                      fillOpacity={0.12}
                      stroke={g.color}
                      strokeWidth={1.4}
                      rx={2}
                    />
                    <line
                      x1={cx - boxW / 2}
                      x2={cx + boxW / 2}
                      y1={y(g.median)}
                      y2={y(g.median)}
                      stroke={g.color}
                      strokeWidth={2}
                    />
                  </>
                )}

                {points
                  .filter((p) => p.group === g.group)
                  .map((p) => (
                    <circle
                      key={p.sample_id}
                      cx={cx + jitter(p.sample_id, boxW * 0.85)}
                      cy={y(p.value)}
                      r={2.6}
                      fill={g.color}
                      fillOpacity={0.55}
                      stroke="var(--surface)"
                      strokeWidth={0.5}
                      onMouseEnter={(e) => {
                        const rect = (e.target as SVGCircleElement).ownerSVGElement!.getBoundingClientRect();
                        setHover({
                          x: e.clientX - rect.left,
                          y: e.clientY - rect.top,
                          sample: p.sample_id,
                          value: p.value,
                        });
                      }}
                      onMouseLeave={() => setHover(null)}
                      className="cursor-pointer"
                    />
                  ))}

                {rotateLabels ? (
                  <text
                    x={cx}
                    y={margin.top + plotH + 12}
                    textAnchor="end"
                    fontSize={11}
                    fontWeight={600}
                    fill="var(--ink)"
                    transform={`rotate(-45 ${cx} ${margin.top + plotH + 12})`}
                  >
                    <title>{`${g.group} (n=${g.values.length})`}</title>
                    {`${g.group} (${g.values.length})`}
                  </text>
                ) : (
                  <>
                    <text
                      x={cx}
                      y={margin.top + plotH + 18}
                      textAnchor="middle"
                      fontSize={11}
                      fontWeight={600}
                      fill="var(--ink)"
                    >
                      <title>{g.group}</title>
                      {fitLabel(g.group, bandW - 8)}
                    </text>
                    <text
                      x={cx}
                      y={margin.top + plotH + 32}
                      textAnchor="middle"
                      fontSize={9.5}
                      fontFamily="var(--f-mono)"
                      className="fill-[var(--ink-mute)]"
                    >
                      n={g.values.length}
                    </text>
                  </>
                )}
              </g>
            );
          })}
        </svg>
        <HoverCard hover={hover} />
      </div>
    </div>
  );
}

function jitter(s: string, amplitude: number): number {
  let h = 0;
  for (let c = 0; c < s.length; c++) h = (h * 31 + s.charCodeAt(c)) % 1000;
  return (h / 1000 - 0.5) * amplitude;
}

// Truncate to the pixel budget actually available rather than a fixed
// character count, and keep the *end* of a name when the start is a shared
// prefix -- several DepMap subtypes differ only after 40 identical
// characters, so head-truncation renders them all identical.
function fitLabel(label: string, pxBudget: number): string {
  const perChar = 6.1;
  const maxChars = Math.max(6, Math.floor(pxBudget / perChar));
  if (label.length <= maxChars) return label;
  const keepEnd = Math.floor(maxChars * 0.6);
  const keepStart = maxChars - keepEnd - 1;
  return `${label.slice(0, keepStart)}…${label.slice(label.length - keepEnd)}`;
}

function ScaleToggle({
  logScale,
  setLogScale,
  auto,
}: {
  logScale: boolean;
  setLogScale: (v: boolean) => void;
  auto: boolean;
}) {
  return (
    <div className="mb-2 flex items-center justify-end gap-1.5">
      {auto && (
        <span className="mr-1 font-mono text-[10px] text-ink-mute" title="Log scale was selected automatically because this data is right-skewed">
          auto
        </span>
      )}
      <button
        type="button"
        onClick={() => setLogScale(false)}
        className={`rounded-[3px] px-2 py-1 font-mono text-[10.5px] transition-colors ${
          !logScale ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
        }`}
      >
        Linear
      </button>
      <button
        type="button"
        onClick={() => setLogScale(true)}
        className={`rounded-[3px] px-2 py-1 font-mono text-[10.5px] transition-colors ${
          logScale ? "bg-accent-soft text-accent-ink" : "text-ink-mute hover:text-ink-soft"
        }`}
      >
        Log scale
      </button>
    </div>
  );
}

function HoverCard({
  hover,
}: {
  hover: { x: number; y: number; sample: string; value: number } | null;
}) {
  if (!hover) return null;
  return (
    <div
      className="pointer-events-none absolute z-10 rounded-[3px] border border-rule bg-surface px-2.5 py-1.5 shadow-lg"
      style={{ left: hover.x + 12, top: hover.y - 10 }}
    >
      <p className="font-mono text-[10.5px] text-ink-mute">{hover.sample}</p>
      <p className="font-mono text-[12px] font-semibold text-ink">{hover.value.toFixed(3)}</p>
    </div>
  );
}
