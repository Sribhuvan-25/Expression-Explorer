/**
 * PCA colours groups with GROUP_COLORS in sorted-group order. --accent
 * (#2455A4) and --cool (#2C6E9B) are both blue; with --cool third,
 * GSE66586's three MYCN groups showed two of them in near-identical blue,
 * and the exported PNG had no legend at all.
 */
import { describe, expect, it } from "vitest";

import { GROUP_COLORS, NO_GROUP_COLOR, scatterColorKey } from "../ScatterPlot";

describe("scatter group colours", () => {
  it("never puts the two blues among the first three groups", () => {
    const firstThree = GROUP_COLORS.slice(0, 3);
    expect(firstThree.includes("var(--accent)") && firstThree.includes("var(--cool)")).toBe(false);
    expect(new Set(firstThree).size).toBe(3);
  });

  it("produces a text colour key in the same order the plot assigns colours", () => {
    const pts = ["non-malignant", "MYCN-amplified", "MYCN non-amplified", "MYCN-amplified"].map((group) => ({ group }));
    expect(scatterColorKey(pts)).toBe(
      "blue = MYCN non-amplified (n=1), red = MYCN-amplified (n=2), amber = non-malignant (n=1)",
    );
    expect(scatterColorKey([{ group: null }])).toBeNull();
  });

  it("names samples with no value in the key instead of folding them into a group", () => {
    // TARGET by etp_status: most samples have no ETP call.
    const pts: { group: string | null }[] = [...[...Array(3)].map(() => ({ group: "ETP" })), ...[...Array(5)].map(() => ({ group: null }))];
    expect(scatterColorKey(pts)).toBe("blue = ETP (n=3), grey = no value (n=5)");
    expect(GROUP_COLORS).not.toContain(NO_GROUP_COLOR);
  });
});

describe("correlation colouring default", () => {
  it("picks a grouping with 2-8 groups among the plotted samples", async () => {
    const { pickColorColumn } = await import("../../pages/CorrelationPage");
    const pts = Array.from({ length: 12 }, (_, i) => ({
      group_columns: { lineage: "Lymphoid", subtype: ["A", "B", "C"][i % 3], cell_line_name: `L${i}` },
    }));
    expect(pickColorColumn(["lineage", "subtype", "cell_line_name"], pts)).toBe("subtype");
    expect(pickColorColumn(["lineage", "subtype"], pts, "etp_status")).toBe("subtype");
    expect(pickColorColumn(["lineage", "cell_line_name"], pts)).toBe("");
  });
});
