import { useMemo, useState } from "react";
import type { ComparePoint } from "../lib/api";

const GROUP_COLORS = ["var(--cool)", "var(--hot)", "var(--accent)", "var(--warn)", "#8B6DBF", "#5C8A6B"];

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

  const groups = useMemo<GroupStats[]>(() => {
    const byGroup = new Map<string, number[]>();
    for (const p of points) {
      if (!byGroup.has(p.group)) byGroup.set(p.group, []);
      byGroup.get(p.group)!.push(p.value);
    }
    return Array.from(byGroup.entries()).map(([group, values], i) => {
      const sorted = [...values].sort((a, b) => a - b);
      return {
        group,
        values,
        min: sorted[0],
        q1: quantile(sorted, 0.25),
        median: quantile(sorted, 0.5),
        q3: quantile(sorted, 0.75),
        max: sorted[sorted.length - 1],
        color: GROUP_COLORS[i % GROUP_COLORS.length],
      };
    });
  }, [points]);

  const width = 640;
  const height = 340;
  const margin = { top: 20, right: 24, bottom: 44, left: 52 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  const allValues = points.map((p) => p.value);
  const yMax = Math.max(...allValues, 1) * 1.08;
  const yMin = Math.min(0, Math.min(...allValues));
  const y = (v: number) => margin.top + plotH - ((v - yMin) / (yMax - yMin)) * plotH;

  const bandW = plotW / groups.length;
  const boxW = Math.min(64, bandW * 0.38);

  const yTicks = 5;
  const tickValues = Array.from({ length: yTicks + 1 }, (_, i) => yMin + ((yMax - yMin) * i) / yTicks);

  return (
    <div className="relative">
      <svg
        ref={svgRef}
        viewBox={`0 0 ${width} ${height}`}
        className="w-full"
        role="img"
        aria-label={`Box plot of ${valueLabel} by group`}
      >
        {tickValues.map((t, i) => (
          <g key={i}>
            <line x1={margin.left} x2={width - margin.right} y1={y(t)} y2={y(t)} stroke="var(--rule)" strokeWidth={1} />
            <text x={margin.left - 8} y={y(t)} textAnchor="end" dominantBaseline="middle" className="fill-[var(--ink-mute)]" fontSize={10} fontFamily="var(--f-mono)">
              {t.toFixed(1)}
            </text>
          </g>
        ))}
        <line x1={margin.left} x2={margin.left} y1={margin.top} y2={margin.top + plotH} stroke="var(--rule-firm)" strokeWidth={1} />
        <line x1={margin.left} x2={width - margin.right} y1={margin.top + plotH} y2={margin.top + plotH} stroke="var(--rule-firm)" strokeWidth={1} />

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
          const jitterSeed = (s: string) => {
            let h = 0;
            for (let c = 0; c < s.length; c++) h = (h * 31 + s.charCodeAt(c)) % 1000;
            return (h / 1000 - 0.5) * boxW * 0.85;
          };
          return (
            <g key={g.group}>
              <line x1={cx} x2={cx} y1={y(g.max)} y2={y(g.q3)} stroke="var(--ink-mute)" strokeWidth={1.2} />
              <line x1={cx} x2={cx} y1={y(g.min)} y2={y(g.q1)} stroke="var(--ink-mute)" strokeWidth={1.2} />
              <line x1={cx - boxW / 4} x2={cx + boxW / 4} y1={y(g.max)} y2={y(g.max)} stroke="var(--ink-mute)" strokeWidth={1.2} />
              <line x1={cx - boxW / 4} x2={cx + boxW / 4} y1={y(g.min)} y2={y(g.min)} stroke="var(--ink-mute)" strokeWidth={1.2} />

              <rect
                x={cx - boxW / 2}
                y={y(g.q3)}
                width={boxW}
                height={Math.max(1, y(g.q1) - y(g.q3))}
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

              {points
                .filter((p) => p.group === g.group)
                .map((p) => (
                  <circle
                    key={p.sample_id}
                    cx={cx + jitterSeed(p.sample_id)}
                    cy={y(p.value)}
                    r={2.6}
                    fill={g.color}
                    fillOpacity={0.55}
                    stroke="var(--surface)"
                    strokeWidth={0.5}
                    onMouseEnter={(e) => {
                      const rect = (e.target as SVGCircleElement).ownerSVGElement!.getBoundingClientRect();
                      setHover({ x: e.clientX - rect.left, y: e.clientY - rect.top, sample: p.sample_id, value: p.value });
                    }}
                    onMouseLeave={() => setHover(null)}
                    className="cursor-pointer"
                  />
                ))}

              <text
                x={cx}
                y={margin.top + plotH + 18}
                textAnchor="middle"
                fontSize={11}
                fontWeight={600}
                fill="var(--ink)"
              >
                {g.group}
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
            </g>
          );
        })}
      </svg>
      {hover && (
        <div
          className="pointer-events-none absolute z-10 rounded-[3px] border border-rule bg-surface px-2.5 py-1.5 shadow-lg"
          style={{ left: hover.x + 12, top: hover.y - 10 }}
        >
          <p className="font-mono text-[10.5px] text-ink-mute">{hover.sample}</p>
          <p className="font-mono text-[12px] font-semibold text-ink">{hover.value.toFixed(3)}</p>
        </div>
      )}
    </div>
  );
}
