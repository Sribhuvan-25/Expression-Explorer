import { GROUP_COLORS } from "../components/ScatterPlot";

type Point = { group: string; value: number };

// The same interpolation as BoxPlot's quantile(sorted, 0.5), so the two
// orders agree to the last bit.
function median(values: number[]): number {
  const s = [...values].sort((a, b) => a - b);
  const pos = (s.length - 1) * 0.5;
  const base = Math.floor(pos);
  return s[base + 1] !== undefined ? s[base] + (pos - base) * (s[base + 1] - s[base]) : s[base];
}

/**
 * Group labels in the order the box plot draws them: median descending,
 * ties in order of first appearance (Array.prototype.sort is stable).
 */
export function labelsByMedian(points: Point[]): string[] {
  const byGroup = new Map<string, number[]>();
  for (const p of points) {
    if (!byGroup.has(p.group)) byGroup.set(p.group, []);
    byGroup.get(p.group)!.push(p.value);
  }
  return Array.from(byGroup.entries())
    .map(([group, values]) => ({ group, m: median(values) }))
    .sort((a, b) => b.m - a.m)
    .map((g) => g.group);
}

/**
 * One colour per group label for a whole multi-facet result, so a label
 * has the same colour in every facet. Colouring by position within each
 * facet painted "non-ETP" ochre in TARGET (third by median) and red in
 * GDS4299 (second), in the same result.
 *
 * Labels are taken facet by facet in display order, each facet's in its
 * own median order, and a label keeps the slot of its first appearance.
 * A single facet therefore gets exactly the colours it had before. Pass
 * every facet's full points, not the filtered ones, or hiding a group
 * would recolour the others.
 */
export function groupColorMap(facets: Point[][]): Map<string, string> {
  const map = new Map<string, string>();
  for (const points of facets)
    for (const label of labelsByMedian(points))
      if (!map.has(label)) map.set(label, GROUP_COLORS[map.size % GROUP_COLORS.length]);
  return map;
}
