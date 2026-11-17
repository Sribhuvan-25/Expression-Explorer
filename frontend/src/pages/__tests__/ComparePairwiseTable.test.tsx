/**
 * At 390 px the pairwise table cut off its p and FDR q columns: the shared
 * DataTable had a 28rem minimum width, wider than a phone's panel, so a
 * phone user saw comparisons with no statistics. Values also wrapped mid-
 * number ("2.29e-" / "7"). jsdom does no layout, so this pins the two
 * properties that caused it; the 390 px check itself is done in a browser.
 */
import { render, screen, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { DatasetPanel } from "../ComparePage";
import type { CompareMultiDatasetResult } from "../../lib/api";

const result: CompareMultiDatasetResult = {
  dataset_id: "target_all_p2",
  display_name: "TARGET ALL-P2",
  skipped: false,
  assay_type: "rna_seq",
  expression_unit: "tpm",
  n_dataset_total: 6,
  n_excluded: 0,
  exclusion_reason: null,
  points: [
    { sample_id: "S1", group: "ETP", value: 5 },
    { sample_id: "S2", group: "ETP", value: 6 },
    { sample_id: "S3", group: "non-ETP", value: 1 },
    { sample_id: "S4", group: "non-ETP", value: 2 },
  ],
  pairwise_tests: [
    { group_a: "ETP", group_b: "non-ETP", n_a: 19, n_b: 146, u_stat: 2401, p_value: 2.29e-7, q_value: 6.88e-7 },
  ] as CompareMultiDatasetResult["pairwise_tests"],
  kruskal_wallis: { h_stat: 34.7, p_value: 2.89e-8 },
};

describe("Compare pairwise table on a phone", () => {
  it("fits its container instead of imposing a minimum width", () => {
    render(<DatasetPanel result={result} gene="MYCN" />);
    const table = within(screen.getByRole("region", { name: /Pairwise tests/ })).getByRole("table");
    expect(table.className).not.toMatch(/min-w-/);
  });

  it("never wraps a statistic mid-value", () => {
    render(<DatasetPanel result={result} gene="MYCN" />);
    const region = screen.getByRole("region", { name: /Pairwise tests/ });
    for (const value of ["19 / 146", "2.29e-7", "6.88e-7"]) {
      expect(within(region).getByRole("cell", { name: value })).toHaveClass("whitespace-nowrap");
    }
  });
});
