/**
 * Coloured scatter with samples that have no value in the colouring
 * column. They used to be drawn in --accent, the first group's blue, so on
 * TARGET PCA coloured by etp_status the 279 uncalled samples read as ETP.
 */
import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { GROUP_COLORS, NO_GROUP_COLOR, ScatterPlot } from "../ScatterPlot";

const pts = [
  { sample_id: "a", x: 1, y: 1, group: "ETP" },
  { sample_id: "b", x: 2, y: 2, group: "non-ETP" },
  { sample_id: "c", x: 3, y: 3, group: null },
  { sample_id: "d", x: 4, y: 4, group: null },
];

describe("scatter points with no group value", () => {
  it("draws them grey, underneath, and lists them in the legend", () => {
    const { container } = render(<ScatterPlot points={pts} xLabel="PC1" yLabel="PC2" />);
    const circles = [...container.querySelectorAll("svg circle")];
    const fills = circles.map((c) => c.getAttribute("fill"));
    expect(fills.slice(0, 2)).toEqual([NO_GROUP_COLOR, NO_GROUP_COLOR]);
    expect(fills.slice(2)).toEqual([GROUP_COLORS[0], GROUP_COLORS[1]]);
    expect(screen.getByText("no value")).toBeInTheDocument();
  });

  it("keeps the single accent colour when nothing is grouped", () => {
    const { container } = render(
      <ScatterPlot points={pts.map((p) => ({ ...p, group: null }))} xLabel="A" yLabel="B" />,
    );
    const fills = new Set([...container.querySelectorAll("svg circle")].map((c) => c.getAttribute("fill")));
    expect([...fills]).toEqual(["var(--accent)"]);
    expect(screen.queryByText("no value")).toBeNull();
  });
});
