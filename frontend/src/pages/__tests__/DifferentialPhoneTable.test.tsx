/**
 * At 390 px the Top differential genes table scrolled sideways and hid
 * its Δ median column behind the shared table's 28rem minimum width; an
 * earlier attempt to fit it let gene symbols break mid-name ("SHAN" /
 * "K3"). It now fits its container and keeps symbols and numbers whole.
 * jsdom does no layout, so this pins the properties; the 390 px
 * measurement is done in a browser (METHODS.md 5.12).
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { DifferentialGenesTab } from "../ComparePage";
import { api } from "../../lib/api";

const SLOW = { timeout: 5000 };

describe("Top differential genes table on a phone", () => {
  it("fits its container and never breaks a gene symbol or a number", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({
      datasets: [{ dataset_id: "target_all_p2", display_name: "TARGET", group_columns: ["etp_status"], n_samples: 469 }],
    } as never);
    vi.spyOn(api, "groupValues").mockResolvedValue({
      group_column: "etp_status",
      values: [{ value: "non-ETP", n: 146 }, { value: "near-ETP", n: 25 }],
      unavailable: false,
    });
    vi.spyOn(api, "differential").mockResolvedValue({
      group_column: "etp_status", group_a: "non-ETP", group_b: "near-ETP", n_a: 146, n_b: 25,
      n_genes_tested: 60660, top_n: 50,
      genes: [{ symbol: "SHANK3", p_value: 3.73e-9, q_value: 8.22e-5, median_diff: -6.68 }],
    } as never);
    render(
      <QueryClientProvider client={new QueryClient({ defaultOptions: { queries: { retry: false } } })}>
        <DifferentialGenesTab />
      </QueryClientProvider>,
    );
    await waitFor(() => expect(screen.getByRole("button", { name: /rank genes/i })).toBeEnabled(), SLOW);
    fireEvent.click(screen.getByRole("button", { name: /rank genes/i }));
    const region = await screen.findByRole("region", { name: /Top differential genes/ }, SLOW);
    expect(within(region).getByRole("table").className).not.toMatch(/min-w-/);
    expect(within(region).getByRole("cell", { name: "SHANK3" })).toHaveClass("whitespace-nowrap");
    for (const value of ["3.73e-9", "8.22e-5", "-6.680"]) {
      expect(within(region).getByRole("cell", { name: value })).toHaveClass("whitespace-nowrap");
    }
  });
});
