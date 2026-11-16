/**
 * With n = 2 vs 2 an exact Mann-Whitney cannot give p below 0.333, so
 * "p = 0.33" on GSE324984 read as a negative result. The notice must
 * appear exactly when min_p_achievable > 0.05, and an older response
 * without the field must neither crash nor show it.
 */
import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { DatasetPanel } from "../ComparePage";
import { MinPNotice } from "../../components/MinPNotice";
import type { CompareMultiDatasetResult, PairwiseTest } from "../../lib/api";

function result(tests: PairwiseTest[], perGroup: number): CompareMultiDatasetResult {
  const points = ["A", "B"].flatMap((g) =>
    Array.from({ length: perGroup }, (_, i) => ({ sample_id: `${g}${i}`, group: g, value: i })),
  );
  return {
    dataset_id: "d",
    display_name: "D",
    skipped: false,
    assay_type: "rna_seq",
    expression_unit: "tpm",
    n_dataset_total: points.length,
    n_excluded: 0,
    exclusion_reason: null,
    points,
    pairwise_tests: tests,
    kruskal_wallis: { h_stat: 1, p_value: 0.3 },
  };
}
const test = (n: number, minP?: number): PairwiseTest => ({
  group_a: "A",
  group_b: "B",
  n_a: n,
  n_b: n,
  u_stat: 4,
  p_value: 0.3333,
  q_value: 0.4,
  ...(minP === undefined ? {} : { min_p_achievable: minP }),
});

describe("min-p notice on Expression Compare", () => {
  it("appears for 2 vs 2", () => {
    render(<DatasetPanel result={result([test(2, 1 / 3)], 2)} gene="MYCL" />);
    expect(screen.getByText(/smallest p-value this test can produce is 0\.333/)).toBeInTheDocument();
    expect(screen.getByText(/n = 2 vs 2/)).toBeInTheDocument();
  });

  it("does not appear for large n", () => {
    render(<DatasetPanel result={result([test(40, 1e-20)], 40)} gene="MYCL" />);
    expect(screen.queryByText(/smallest p-value/)).not.toBeInTheDocument();
  });

  it("does not crash or appear when the field is absent", () => {
    render(<DatasetPanel result={result([test(2)], 2)} gene="MYCL" />);
    expect(screen.queryByText(/smallest p-value/)).not.toBeInTheDocument();
  });
});

/**
 * GSE324984 grouped by `arm` yields 6 pairs, 4 of which compare different
 * cell lines -- differences that mostly reflect the cell line, not MYCL.
 * Those rows must be visibly tagged and the caveat shown; a response
 * without the fields (older backend) must render unchanged.
 */
describe("stratum-crossing pairs on Expression Compare", () => {
  const crossing = { column: "cell_line", group_a: "C4-2B", group_b: "NCI-H660" };
  it("tags only the crossing pair and shows the caveat", () => {
    const r = result(
      [
        { ...test(2, 1 / 3), group_a: "C4-2B OE", group_b: "C4-2B EV", crosses_stratum: null },
        { ...test(2, 1 / 3), group_a: "C4-2B OE", group_b: "H660 shCtrl", crosses_stratum: crossing },
      ],
      2,
    );
    r.stratum_caveat = "Some pairs compare groups from different cell_line.";
    render(<DatasetPanel result={r} gene="MYCL" />);
    const tagged = screen.getAllByTestId("crosses-stratum");
    expect(tagged).toHaveLength(1);
    expect(tagged[0].textContent).toContain("H660 shCtrl");
    expect(tagged[0].textContent).toMatch(/different cell_line/);
    expect(screen.getByTestId("stratum-caveat").textContent).toMatch(/different cell_line/);
  });

  it("renders without the fields (older responses)", () => {
    render(<DatasetPanel result={result([test(2, 1 / 3)], 2)} gene="MYCL" />);
    expect(screen.queryByTestId("crosses-stratum")).not.toBeInTheDocument();
    expect(screen.queryByTestId("stratum-caveat")).not.toBeInTheDocument();
  });
});

/**
 * The genome-wide differential ranking reports p = 0.22 on a 2-vs-2 design
 * (and 0.079 on 5 vs 2), below the exact floor of 0.333 / 0.095 -- an
 * approximation. The notice must not state a floor the table then
 * contradicts without saying why.
 */
describe("p-values shown below the exact floor", () => {
  const sizes = [{ n_a: 2, n_b: 2, min_p_achievable: 1 / 3 }];
  it("explains them when a shown p is below the floor", () => {
    render(<MinPNotice tests={sizes} smallestShownP={0.2207} />);
    expect(screen.getByTestId("below-floor").textContent).toMatch(/normal approximation/);
  });
  it("says nothing extra when shown p-values respect the floor", () => {
    render(<MinPNotice tests={sizes} smallestShownP={1 / 3} />);
    expect(screen.getByText(/smallest p-value/)).toBeInTheDocument();
    expect(screen.queryByTestId("below-floor")).not.toBeInTheDocument();
  });
});

describe("underpoweredSizes", () => {
  it("treats mirrored group sizes as one statement, largest floor first", async () => {
    const { underpoweredSizes } = await import("../../components/MinPNotice");
    const sizes = underpoweredSizes([
      { n_a: 4, n_b: 3, min_p_achievable: 2 / 35 },
      { n_a: 3, n_b: 4, min_p_achievable: 2 / 35 },
      { n_a: 2, n_b: 2, min_p_achievable: 1 / 3 },
      { n_a: 30, n_b: 29, min_p_achievable: null },
    ]);
    expect(sizes.map((s) => s.label)).toEqual(["2 vs 2", "4 vs 3"]);
  });
});
