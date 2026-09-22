/**
 * Guards the "single-valued grouping column" gap on the Top Differential
 * Genes tab.
 *
 * When a grouping column has fewer than 2 distinct values -- e.g.
 * DepMap's own `lineage` is filtered to "Lymphoid" only upstream (see
 * app/ingest/depmap.py) -- the "Rank genes" button and both group
 * selects go disabled (groupOptions.length < 2), but groupA/groupB never
 * get populated (the auto-select effect requires >= 2 options), so the
 * pre-existing "Group A and Group B must be different" hint never fires
 * either (it requires groupA to be truthy first). The result: a user
 * lands on the tab with a disabled button and zero text anywhere on the
 * page saying why. Found during regression testing, 2026-09-21.
 *
 * Distinct from the §7.10 "unavailable" case (a column that's declared
 * but has zero CURRENT values, e.g. a degraded upstream fetch): this is
 * the structural case where the column genuinely has exactly one value,
 * every time, for this dataset.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { DifferentialGenesTab } from "../ComparePage";
import { api } from "../../lib/api";
import type { DatasetSummary } from "../../lib/api";

const SINGLE_VALUED_DATASET: DatasetSummary = {
  dataset_id: "depmap",
  display_name: "DepMap cell line RNA-seq",
  group_columns: ["lineage"],
  supports_survival: false,
  n_samples: 186,
  assay_type: "rna_seq",
};

const TWO_VALUED_DATASET: DatasetSummary = {
  dataset_id: "target_all_p2",
  display_name: "TARGET ALL-P2 (paediatric T-ALL)",
  group_columns: ["vital_status"],
  supports_survival: true,
  n_samples: 469,
  assay_type: "rna_seq",
};

function renderTab() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false, gcTime: 0 } } });
  return render(
    <QueryClientProvider client={client}>
      <DifferentialGenesTab />
    </QueryClientProvider>,
  );
}

describe("DifferentialGenesTab single-value grouping column", () => {
  it("explains why Rank genes is disabled when the column has only one value", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [SINGLE_VALUED_DATASET] });
    vi.spyOn(api, "groupValues").mockResolvedValue({
      group_column: "lineage",
      values: [{ value: "Lymphoid", n: 186 }],
      unavailable: false,
    });

    renderTab();

    // The exact defect: previously nothing on the page said why.
    const message = await screen.findByText(/only one value on this dataset/i);
    expect(message).toBeInTheDocument();
    expect(message.textContent).toContain("Lymphoid");
    expect(message.textContent).toContain("lineage");

    await waitFor(() => {
      expect(screen.getByRole("button", { name: /rank genes/i })).toBeDisabled();
    });

    // Must not be conflated with the "must be different" case -- that
    // message requires groupA to be set, which it never is here.
    expect(screen.queryByText(/must be different/i)).not.toBeInTheDocument();
    // Must not be conflated with the §7.10 "unavailable" case either --
    // this column genuinely has one real value, not zero.
    expect(screen.queryByText(/may be temporarily unavailable/i)).not.toBeInTheDocument();
  });

  it("shows no message and enables Rank genes once the column has two values", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [TWO_VALUED_DATASET] });
    vi.spyOn(api, "groupValues").mockResolvedValue({
      group_column: "vital_status",
      values: [
        { value: "Alive", n: 358 },
        { value: "Dead", n: 106 },
      ],
      unavailable: false,
    });

    renderTab();

    await waitFor(() => {
      expect(screen.getByRole("button", { name: /rank genes/i })).toBeEnabled();
    });
    expect(screen.queryByText(/only one value on this dataset/i)).not.toBeInTheDocument();
  });

  it("still shows the §7.10 unavailable message for a declared-but-empty column, not the single-value one", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [SINGLE_VALUED_DATASET] });
    vi.spyOn(api, "groupValues").mockResolvedValue({
      group_column: "lineage",
      values: [],
      unavailable: true,
    });

    renderTab();

    expect(await screen.findByText(/may be temporarily unavailable/i)).toBeInTheDocument();
    expect(screen.queryByText(/only one value on this dataset/i)).not.toBeInTheDocument();
  });
});
