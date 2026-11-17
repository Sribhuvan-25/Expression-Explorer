/**
 * Exported PNGs carry a text colour key ("teal = ETP (n=19), ..."), since
 * the HTML legend is not part of the image. The names are hard-coded next
 * to the palette, so when the Workbench redesign changed --accent from
 * blue (#2455A4) to teal (#0B6B66) the key kept saying "blue" for teal
 * points -- and "teal" for the actually blue --cool. This resolves every
 * palette entry to the colour index.css really renders and checks its
 * hue against the name the key prints.
 */
// @ts-expect-error -- node built-in; the frontend tsconfig carries no node types,
// and vitest's css: false makes "?raw" imports of CSS come back empty.
import { readFileSync } from "node:fs";
import { describe, expect, it } from "vitest";

import { GROUP_COLORS, GROUP_COLOR_NAMES } from "../ScatterPlot";
import { CURVE_COLORS, CURVE_COLOR_NAMES, FALLBACK_COLORS, FALLBACK_COLOR_NAMES } from "../SurvivalCurve";

const css: string = readFileSync(`${(globalThis as unknown as { process: { cwd(): string } }).process.cwd()}/src/index.css`, "utf8");
const rootBlock = css.slice(css.indexOf(":root"), css.indexOf("}", css.indexOf(":root")));

function hex(color: string): string {
  const m = /^var\((--[\w-]+)\)$/.exec(color);
  if (!m) return color;
  const def = new RegExp(`${m[1]}:\\s*(#[0-9A-Fa-f]{6})`).exec(rootBlock);
  if (!def) throw new Error(`${m[1]} not defined in :root`);
  return def[1];
}

function hue(h: string): number {
  const [r, g, b] = [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
  const max = Math.max(r, g, b);
  const d = max - Math.min(r, g, b);
  if (d === 0) return 0;
  const raw = max === r ? ((g - b) / d) % 6 : max === g ? (b - r) / d + 2 : (r - g) / d + 4;
  return (raw * 60 + 360) % 360;
}

// Hue ranges (degrees) each printed name may cover.
const RANGES: Record<string, [number, number][]> = {
  red: [[0, 22], [340, 360]],
  amber: [[30, 55]],
  green: [[75, 165]],
  teal: [[165, 192]],
  blue: [[192, 245]],
  purple: [[245, 300]],
};
const named = (name: string, color: string) =>
  (RANGES[name] ?? []).some(([lo, hi]) => hue(hex(color)) >= lo && hue(hex(color)) <= hi);

describe("colour names printed in exported keys", () => {
  it("scatter/PCA: each palette colour is called what it looks like", () => {
    GROUP_COLORS.forEach((c, i) => expect(`${GROUP_COLOR_NAMES[i]}: ${named(GROUP_COLOR_NAMES[i], c)}`).toBe(`${GROUP_COLOR_NAMES[i]}: true`));
  });

  it("survival: HIGH/LOW and fallback curves are called what they look like", () => {
    for (const [label, c] of Object.entries(CURVE_COLORS)) expect(named(CURVE_COLOR_NAMES[label], c)).toBe(true);
    FALLBACK_COLORS.forEach((c, i) => expect(`${FALLBACK_COLOR_NAMES[i]}: ${named(FALLBACK_COLOR_NAMES[i], c)}`).toBe(`${FALLBACK_COLOR_NAMES[i]}: true`));
  });
});
