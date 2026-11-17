/**
 * A group label is one colour in every facet of a Compare result.
 *
 * BoxPlot coloured groups by their median rank within each plot, so with
 * MYCN by etp_status across TARGET-ALL-P2 (ETP > near-ETP > non-ETP) and
 * GDS4299 (ETP > non-ETP) "non-ETP" was the third colour (ochre) in one
 * facet and the second (red) in the other -- the reader matched groups
 * across datasets by colour and got them wrong. Rendered through the real
 * page so the shared map has to reach every facet.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { ComparePage } from "../ComparePage";
import { api, type CompareMultiDatasetResult } from "../../lib/api";
import { WorkbenchProvider } from "../../lib/workbench";

function facet(id: string, groups: [string, number[]][]): CompareMultiDatasetResult {
  return {
    dataset_id: id,
    display_name: id,
    skipped: false,
    assay_type: "rna_seq",
    expression_unit: "log2_tpm",
    n_dataset_total: 0,
    n_excluded: 0,
    exclusion_reason: null,
    points: groups.flatMap(([group, values]) => values.map((value, i) => ({ sample_id: `${id}-${group}-${i}`, group, value }))),
    pairwise_tests: [],
    kruskal_wallis: { h_stat: 1, p_value: 0.01 },
  };
}

function boxFill(panel: Element, group: string): string | null {
  return panel.querySelector(`g[data-group="${group}"] rect[fill-opacity]`)?.getAttribute("fill") ?? null;
}

describe("Expression Compare group colours across facets", () => {
  it("gives a label the same colour in every dataset's box plot", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({
      datasets: [
        { dataset_id: "target_all_p2", display_name: "TARGET", group_columns: ["etp_status"], n_samples: 469 },
        { dataset_id: "gds4299", display_name: "GDS4299", group_columns: ["etp_status"], n_samples: 52 },
      ],
    } as never);
    vi.spyOn(api, "geneInfo").mockResolvedValue(null as never);
    vi.spyOn(api, "compareMulti").mockResolvedValue({
      gene: "MYCN",
      group_column: "etp_status",
      datasets: [
        facet("target_all_p2", [
          ["ETP", [8, 9, 10]],
          ["near-ETP", [5, 6, 7]],
          ["non-ETP", [1, 2, 3]],
        ]),
        facet("gds4299", [
          ["ETP", [7, 8, 9]],
          ["non-ETP", [3, 4, 5]],
        ]),
      ],
    });

    render(
      <QueryClientProvider client={new QueryClient({ defaultOptions: { queries: { retry: false } } })}>
        <WorkbenchProvider>
          <ComparePage />
        </WorkbenchProvider>
      </QueryClientProvider>,
    );
    const run = screen.getByRole("button", { name: /run comparison/i });
    await waitFor(() => expect(run).toBeEnabled());
    fireEvent.click(run);

    await waitFor(() => expect(document.querySelectorAll("svg[role=img]").length).toBe(2));
    const [target, gds] = [...document.querySelectorAll("svg[role=img]")];
    for (const label of ["ETP", "non-ETP"]) {
      expect(boxFill(target, label)).not.toBeNull();
      expect(boxFill(gds, label)).toBe(boxFill(target, label));
    }
    // Still three distinct colours where there are three groups.
    expect(new Set(["ETP", "near-ETP", "non-ETP"].map((g) => boxFill(target, g))).size).toBe(3);
  });
});
