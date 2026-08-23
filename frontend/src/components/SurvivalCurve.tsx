import { useMemo, useState } from "react";
import type { SurvivalCurve } from "../lib/api";

const CURVE_COLORS: Record<string, string> = {
  HIGH: "var(--hot)",
  LOW: "var(--cool)",
};
const FALLBACK_COLORS = ["var(--accent)", "var(--warn)", "#8B6DBF"];

export function SurvivalPlot({ curves }: { curves: Record<string, SurvivalCurve> }) {
  const [hover, setHover] = useState<{ x: number; y: number; label: string; days: number; prob: number } | null>(null);

  const width = 640;
  const height = 340;
  const margin = { top: 20, right: 24, bottom: 44, left: 52 };
  const plotW = width - margin.left - margin.right;
  const plotH = height - margin.top - margin.bottom;

  const maxDays = useMemo(() => {
    let m = 0;
    for (const c of Object.values(curves)) {
      for (const p of c.points) m = Math.max(m, p.days);
    }
    return m || 1;
  }, [curves]);

  const x = (days: number) => margin.left + (days / maxDays) * plotW;
  const y = (prob: number) => margin.top + plotH - prob * plotH;

  const yTicks = [0, 0.25, 0.5, 0.75, 1.0];
  const xTicks = Array.from({ length: 5 }, (_, i) => Math.round((maxDays * i) / 4));

  const entries = Object.entries(curves);

  return (
    <div className="relative">
      <svg viewBox={`0 0 ${width} ${height}`} className="w-full" role="img" aria-label="Kaplan-Meier survival curves">
        {yTicks.map((t) => (
          <g key={t}>
            <line x1={margin.left} x2={width - margin.right} y1={y(t)} y2={y(t)} stroke="var(--rule)" strokeWidth={1} />
            <text x={margin.left - 8} y={y(t)} textAnchor="end" dominantBaseline="middle" fontSize={10} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
              {t.toFixed(2)}
            </text>
          </g>
        ))}
        {xTicks.map((t) => (
          <text key={t} x={x(t)} y={margin.top + plotH + 18} textAnchor="middle" fontSize={10} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]">
            {t}
          </text>
        ))}
        <text x={margin.left + plotW / 2} y={height - 4} textAnchor="middle" fontSize={10.5} fontFamily="var(--f-mono)" className="fill-[var(--ink-mute)]" style={{ letterSpacing: "0.04em" }}>
          DAYS
        </text>
        <text
          transform={`translate(14, ${margin.top + plotH / 2}) rotate(-90)`}
          textAnchor="middle"
          fontSize={10.5}
          fontFamily="var(--f-mono)"
          className="fill-[var(--ink-mute)]"
          style={{ letterSpacing: "0.04em" }}
        >
          SURVIVAL PROBABILITY
        </text>

        <line x1={margin.left} x2={margin.left} y1={margin.top} y2={margin.top + plotH} stroke="var(--rule-firm)" />
        <line x1={margin.left} x2={width - margin.right} y1={margin.top + plotH} y2={margin.top + plotH} stroke="var(--rule-firm)" />

        {entries.map(([label, curve], i) => {
          const color = CURVE_COLORS[label] ?? FALLBACK_COLORS[i % FALLBACK_COLORS.length];
          const sorted = [...curve.points].sort((a, b) => a.days - b.days);
          let path = "";
          sorted.forEach((p, idx) => {
            const px = x(p.days);
            const py = y(p.survival_probability);
            if (idx === 0) {
              path += `M ${x(0)} ${y(1)} L ${px} ${y(1)} L ${px} ${py}`;
            } else {
              const prevY = y(sorted[idx - 1].survival_probability);
              path += ` L ${px} ${prevY} L ${px} ${py}`;
            }
          });
          return (
            <g key={label}>
              <path d={path} fill="none" stroke={color} strokeWidth={2} strokeLinejoin="round" />
              {sorted.map((p, idx) => (
                <circle
                  key={idx}
                  cx={x(p.days)}
                  cy={y(p.survival_probability)}
                  r={3}
                  fill={color}
                  opacity={0}
                  onMouseEnter={(e) => {
                    const rect = (e.target as SVGCircleElement).ownerSVGElement!.getBoundingClientRect();
                    setHover({ x: e.clientX - rect.left, y: e.clientY - rect.top, label, days: p.days, prob: p.survival_probability });
                  }}
                  onMouseLeave={() => setHover(null)}
                  className="cursor-pointer"
                />
              ))}
            </g>
          );
        })}
      </svg>

      <div className="mt-1 flex items-center justify-center gap-5">
        {entries.map(([label, curve], i) => {
          const color = CURVE_COLORS[label] ?? FALLBACK_COLORS[i % FALLBACK_COLORS.length];
          return (
            <span key={label} className="flex items-center gap-1.5 text-[12px] text-ink-soft">
              <span className="h-2.5 w-2.5 rounded-sm" style={{ background: color }} />
              {label}
              <span className="font-mono text-[10.5px] text-ink-mute">n={curve.n}</span>
            </span>
          );
        })}
      </div>

      {hover && (
        <div
          className="pointer-events-none absolute z-10 rounded-[3px] border border-rule bg-surface px-2.5 py-1.5 shadow-lg"
          style={{ left: hover.x + 12, top: hover.y - 10 }}
        >
          <p className="font-mono text-[10.5px] text-ink-mute">{hover.label} &middot; day {hover.days}</p>
          <p className="font-mono text-[12px] font-semibold text-ink">{(hover.prob * 100).toFixed(1)}%</p>
        </div>
      )}
    </div>
  );
}
