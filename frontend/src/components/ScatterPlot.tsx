import { useEffect, useMemo, useRef, useState } from "react";

const GROUP_COLORS = ["var(--accent)", "var(--hot)", "var(--cool)", "var(--warn)", "#8B6DBF", "#5C8A6B"];

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

  const groups = useMemo(() => {
    const seen = new Set<string>();
    for (const p of points) if (p.group != null) seen.add(p.group);
    return Array.from(seen).sort();
  }, [points]);
  const colorFor = (group: string | null | undefined) =>
    group == null ? "var(--accent)" : GROUP_COLORS[groups.indexOf(group) % GROUP_COLORS.length];

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

  const width = Math.max(containerWidth || 640, 320);
  const margin = { top: 16, right: 20, bottom: 46, left: 56 };
  const plotW = Math.max(width - margin.left - margin.right, 160);
  const plotH = Math.round(Math.min(400, Math.max(220, plotW * 0.62)));
  const height = plotH + margin.top + margin.bottom;

  const { xMin, xMax, yMin, yMax } = useMemo(() => {
    const xs = points.map((p) => p.x);
    const ys = points.map((p) => p.y);
    // Pad the domain slightly so points at the extremes aren't drawn
    // flush against the axis line.
    const xPad = (Math.max(...xs) - Math.min(...xs)) * 0.06 || 1;
    const yPad = (Math.max(...ys) - Math.min(...ys)) * 0.06 || 1;
    return {
      xMin: Math.min(...xs) - xPad,
      xMax: Math.max(...xs) + xPad,
      yMin: Math.min(...ys) - yPad,
      yMax: Math.max(...ys) + yPad,
    };
  }, [points]);

  const x = (v: number) => margin.left + ((v - xMin) / (xMax - xMin)) * plotW;
  const y = (v: number) => margin.top + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

  const xTicks = Array.from({ length: 5 }, (_, i) => xMin + ((xMax - xMin) * i) / 4);
  const yTicks = Array.from({ length: 5 }, (_, i) => yMin + ((yMax - yMin) * i) / 4);

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
              {t.toFixed(1)}
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
            {t.toFixed(1)}
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

        {points.map((p) => {
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

      {groups.length > 1 && (
        <div className="mt-1 flex flex-wrap items-center justify-center gap-4">
          {groups.map((g) => (
            <span key={g} className="flex items-center gap-1.5 text-[12px] text-ink-soft">
              <span className="h-2.5 w-2.5 rounded-sm" style={{ background: colorFor(g) }} />
              {g}
            </span>
          ))}
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
