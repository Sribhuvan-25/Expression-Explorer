/**
 * A dataset can now be "real but not queryable yet" -- the backend warms
 * its cache on a background thread so the deploy healthcheck never waits
 * on a multi-minute download (METHODS.md 7.9).
 *
 * That state has to be *said*. A warming dataset has no sample count, so
 * without explicit handling the sidebar shows a name above a blank line,
 * which reads as a broken dataset rather than a pending one.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { DatasetStatus } from "../DatasetStatus";
import { api } from "../../lib/api";
import type { DatasetSummary } from "../../lib/api";

function renderWithQuery() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false, gcTime: 0 } } });
  return render(
    <QueryClientProvider client={client}>
      <DatasetStatus />
    </QueryClientProvider>,
  );
}

const READY: DatasetSummary = {
  dataset_id: "depmap",
  display_name: "DepMap cell line RNA-seq",
  group_columns: ["lineage"],
  supports_survival: false,
  n_samples: 186,
  assay_type: "rna_seq",
  warming: false,
};

const WARMING: DatasetSummary = {
  dataset_id: "target_all_p2",
  display_name: "TARGET ALL-P2",
  group_columns: ["etp_status"],
  supports_survival: true,
  warming: true,
};

describe("DatasetStatus warming state", () => {
  it("says a warming dataset is preparing rather than showing a blank line", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [WARMING] });
    renderWithQuery();
    expect(await screen.findByText(/preparing/i)).toBeInTheDocument();
  });

  it("shows provenance for a ready dataset and does not call it preparing", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [READY] });
    renderWithQuery();
    expect(await screen.findByText(/186 samples/)).toBeInTheDocument();
    expect(screen.queryByText(/preparing/i)).not.toBeInTheDocument();
  });

  it("shows both states at once during a partial warm-up", async () => {
    // The realistic mid-deploy picture: one dataset usable, another still
    // downloading. Both facts have to be visible simultaneously.
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [READY, WARMING] });
    renderWithQuery();
    expect(await screen.findByText(/186 samples/)).toBeInTheDocument();
    expect(screen.getByText(/preparing/i)).toBeInTheDocument();
    expect(screen.getByText("DepMap cell line RNA-seq")).toBeInTheDocument();
    expect(screen.getByText("TARGET ALL-P2")).toBeInTheDocument();
  });
});
