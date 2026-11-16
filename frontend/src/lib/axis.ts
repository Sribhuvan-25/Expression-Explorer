/** Shared axis helpers for the SVG charts. */

/**
 * Linear-axis ticks on round values (1, 2, 2.5 or 5 x 10^k), with the
 * domain widened to the outermost ticks so gridlines meet the plot edges.
 * Exported for tests.
 */
export function niceTicks(lo: number, hi: number, count = 5): number[] {
  if (!(hi > lo)) {
    const pad = Math.abs(lo) * 0.05 || 1;
    lo -= pad;
    hi += pad;
  }
  const raw = (hi - lo) / count;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((st) => st >= raw * (1 - 1e-9))!;
  const start = Math.floor(lo / step + 1e-9) * step;
  const n = Math.ceil(hi / step - 1e-9) - Math.round(start / step);
  // Round away float noise (0.30000000000000004) at the step's precision.
  const digits = Math.max(0, -Math.floor(Math.log10(step)) + 2);
  return Array.from({ length: n + 1 }, (_, i) => Number((start + i * step).toFixed(digits)));
}

/**
 * Decimals needed to print ticks spaced `step` apart without repeats or
 * noise: 1 -> 0, 0.5 -> 1, 2.5 -> 1, 0.25 -> 2, 20 -> 0.
 */
export function tickDigits(step: number): number {
  const e = Math.floor(Math.log10(step) + 1e-9);
  const mantissa = Math.round((step / 10 ** e) * 100) / 100;
  return Math.min(4, Math.max(0, -e + (mantissa === 2.5 ? 1 : 0)));
}
