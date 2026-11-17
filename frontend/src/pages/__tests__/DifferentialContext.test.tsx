/**
 * Top differential genes adopts the context bar's dataset. Its query panel
 * collapses after a run; when the context then changed the dataset, the
 * panel stayed collapsed over inputs that no longer matched the result.
 * It must reopen, as every other analysis does.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { act, fireEvent, render, screen } from "@testing-library/react";
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

// Deterministic waiting. This test used to poll with waitFor against a
// wall-clock budget and failed intermittently under parallel test load:
// the page's own chain (datasets query, dataset and column effects,
// group-values query, group effect) takes ~100 ms, but waitFor's first
// getByRole(..., { name }) costs ~240 ms on a cold jsdom and blocks the
// event loop, so on a starved CPU the 1 s budget ran out (reproduced 2 of
// 8 concurrent runs). Now each wait advances the event loop one task at a
// time inside act() until the condition holds, bounded by a number of
// turns, not milliseconds -- the outcome no longer depends on CPU speed --
// and the conditions use cheap queries.
async function settleUntil(check: () => boolean, maxTurns = 500): Promise<void> {
  for (let i = 0; i < maxTurns && !check(); i++) {
    await act(async () => {
      await new Promise((r) => setTimeout(r, 0));
    });
  }
}
const rankButton = () => screen.queryByText(/^rank genes$/i)?.closest("button") ?? null;

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
    await settleUntil(() => rankButton()?.disabled === false);
    expect(rankButton()).toBeEnabled();
    fireEvent.click(rankButton()!);
    await settleUntil(() => screen.queryByTestId("query-summary") !== null);
    expect(screen.getByTestId("query-summary")).toHaveTextContent("TARGET");

    await act(async () => setContext!({ datasetId: "gds4299" }));
    await settleUntil(() => screen.queryByTestId("query-summary") === null);
    expect(screen.queryByTestId("query-summary")).toBeNull();
  });
});
