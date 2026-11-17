/**
 * Top differential genes adopts the context bar's dataset. Its query panel
 * collapses after a run; when the context then changed the dataset, the
 * panel stayed collapsed over inputs that no longer matched the result.
 * It must reopen, as every other analysis does.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { act, fireEvent, render, screen, waitFor } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { DifferentialGenesTab } from "../ComparePage";
import { api } from "../../lib/api";
import { useWorkbench, WorkbenchProvider } from "../../lib/workbench";

let setContext: ((n: { datasetId?: string; genes?: string[] }) => void) | undefined;
function Grab() {
  setContext = useWorkbench()?.setContext;
  return null;
}

afterEach(() => localStorage.clear());

// The page loads datasets and group values before the button enables;
// under a heavily loaded machine that took longer than waitFor's 1 s
// default and the test failed intermittently, with nothing wrong.
const SLOW = { timeout: 5000 };

describe("differential genes with the context bar", () => {
  it("reopens the collapsed query when the context dataset changes", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({
      datasets: [
        { dataset_id: "target_all_p2", display_name: "TARGET", group_columns: ["etp_status"], n_samples: 469 },
        { dataset_id: "gds4299", display_name: "GDS4299", group_columns: ["etp_status"], n_samples: 52 },
      ],
    } as never);
    vi.spyOn(api, "groupValues").mockResolvedValue({
      group_column: "etp_status",
      values: [{ value: "ETP", n: 19 }, { value: "non-ETP", n: 146 }],
      unavailable: false,
    });
    vi.spyOn(api, "differential").mockResolvedValue({
      group_column: "etp_status", group_a: "ETP", group_b: "non-ETP", n_a: 19, n_b: 146,
      n_genes_tested: 100, top_n: 50, genes: [],
    } as never);
    render(
      <QueryClientProvider client={new QueryClient({ defaultOptions: { queries: { retry: false } } })}>
        <WorkbenchProvider>
          <Grab />
          <DifferentialGenesTab />
        </WorkbenchProvider>
      </QueryClientProvider>,
    );
    await waitFor(() => expect(screen.getByRole("button", { name: /rank genes/i })).toBeEnabled(), SLOW);
    fireEvent.click(screen.getByRole("button", { name: /rank genes/i }));
    expect(await screen.findByTestId("query-summary", {}, SLOW)).toHaveTextContent("TARGET");

    await act(async () => setContext!({ datasetId: "gds4299" }));
    await waitFor(() => expect(screen.queryByTestId("query-summary")).toBeNull(), SLOW);
  });
});
