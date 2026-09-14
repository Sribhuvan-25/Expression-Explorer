/**
 * The binning rule decides whether the score histogram is readable.
 *
 * A fixed 36 bins left 39% of bins empty on the 52-sample GDS4299 cohort
 * -- a gap-toothed comb that cannot answer the question the chart exists
 * to answer. These pin the behaviour against the real cohort sizes this
 * app ships with, so a future "simplification" back to a constant fails
 * loudly. See METHODS.md 5.6.
 */
import { describe, expect, it } from "vitest";

import { histogramBinCount } from "../SignaturePage";

/** Bin the scores the way ScoreDistribution does, and report the gaps. */
function emptyBinFraction(scores: number[], bins: number): number {
  const min = Math.min(...scores);
  const max = Math.max(...scores);
  const counts = new Array(bins).fill(0);
  for (const s of scores) {
    counts[Math.min(bins - 1, Math.floor(((s - min) / (max - min)) * bins))] += 1;
  }
  return counts.filter((c) => c === 0).length / bins;
}

/** Deterministic spread-out scores — no RNG, so failures are reproducible. */
function evenScores(n: number): number[] {
  return Array.from({ length: n }, (_, i) => i / (n - 1));
}

/** The zero-inflated shape a single-gene signature produces (METHODS.md 3.1a). */
function zeroInflated(n: number, zeroFraction: number): number[] {
  const zeros = Math.round(n * zeroFraction);
  return [
    ...Array.from({ length: zeros }, () => 0),
    ...Array.from({ length: n - zeros }, (_, i) => (i + 1) / (n - zeros)),
  ];
}

describe("histogramBinCount", () => {
  it("stays within the legible clamp for every cohort size we ship", () => {
    for (const n of [4, 12, 52, 186, 469, 5000]) {
      const bins = histogramBinCount(evenScores(n));
      expect(bins, `n=${n}`).toBeGreaterThanOrEqual(6);
      expect(bins, `n=${n}`).toBeLessThanOrEqual(36);
      expect(Number.isInteger(bins), `n=${n}`).toBe(true);
    }
  });

  it("leaves far fewer empty bins than a fixed 36 on the small cohort", () => {
    // GDS4299: 52 samples. Real signature scores clump -- they are not
    // evenly spread -- which is what leaves a fixed 36 bins gap-toothed.
    // Perfectly even scores would fill all 36 and hide the problem, so
    // this uses a clumped shape: most scores in a narrow band, a few
    // stragglers pulling the range wide.
    const clumped = [
      ...Array.from({ length: 46 }, (_, i) => 0.30 + (i / 45) * 0.12),
      0.02, 0.05, 0.71, 0.83, 0.94, 1.0,
    ];
    expect(clumped).toHaveLength(52);

    const adaptive = histogramBinCount(clumped);
    const adaptiveGaps = emptyBinFraction(clumped, adaptive);
    const fixedGaps = emptyBinFraction(clumped, 36);

    expect(fixedGaps).toBeGreaterThan(0.3); // the comb this replaced
    expect(adaptiveGaps).toBeLessThan(fixedGaps);
  });

  it("keeps meaningful resolution on the large cohort", () => {
    // TARGET: 469 samples, which rendered perfectly well at 36 bins.
    // Freedman-Diaconis alone collapsed this to 9 -- fixing the small
    // cohort by flattening the large one would trade one defect for
    // another, which is what the sqrt(n) floor prevents.
    expect(histogramBinCount(evenScores(469))).toBeGreaterThanOrEqual(20);
  });

  it("gives a larger cohort at least as many bins as a smaller one", () => {
    const sizes = [12, 52, 186, 469];
    const bins = sizes.map((n) => histogramBinCount(evenScores(n)));
    for (let i = 1; i < bins.length; i += 1) {
      expect(bins[i], `n=${sizes[i]} vs n=${sizes[i - 1]}`).toBeGreaterThanOrEqual(bins[i - 1]);
    }
  });

  it("survives a zero-inflated distribution without degenerating", () => {
    // Half the cohort tied at exactly 0 is the normal shape for a
    // single-gene signature, and it is where a naive IQR rule breaks.
    const bins = histogramBinCount(zeroInflated(466, 0.47));
    expect(bins).toBeGreaterThanOrEqual(6);
    expect(bins).toBeLessThanOrEqual(36);
  });

  it("falls back cleanly when the IQR is zero", () => {
    // Every score identical except one: IQR collapses to 0 and the
    // Freedman-Diaconis width would divide by zero.
    const scores = [...Array.from({ length: 40 }, () => 0.5), 1];
    const bins = histogramBinCount(scores);
    expect(Number.isFinite(bins)).toBe(true);
    expect(bins).toBeGreaterThanOrEqual(6);
    expect(bins).toBeLessThanOrEqual(36);
  });
});
