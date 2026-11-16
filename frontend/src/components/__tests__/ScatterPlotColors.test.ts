/**
 * PCA colours groups with GROUP_COLORS in sorted-group order. --accent
 * (#2455A4) and --cool (#2C6E9B) are both blue; with --cool third,
 * GSE66586's three MYCN groups showed two of them in near-identical blue,
 * and the exported PNG had no legend at all.
 */
import { describe, expect, it } from "vitest";

import { GROUP_COLORS, scatterColorKey } from "../ScatterPlot";

describe("scatter group colours", () => {
  it("never puts the two blues among the first three groups", () => {
    const firstThree = GROUP_COLORS.slice(0, 3);
    expect(firstThree.includes("var(--accent)") && firstThree.includes("var(--cool)")).toBe(false);
    expect(new Set(firstThree).size).toBe(3);
  });

  it("produces a text colour key in the same order the plot assigns colours", () => {
    const pts = ["non-malignant", "MYCN-amplified", "MYCN non-amplified", "MYCN-amplified"].map((group) => ({ group }));
    expect(scatterColorKey(pts)).toBe("blue = MYCN non-amplified, red = MYCN-amplified, amber = non-malignant");
    expect(scatterColorKey([{ group: null }])).toBeNull();
  });
});
