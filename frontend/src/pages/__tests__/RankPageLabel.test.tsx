import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { RankPage } from "../RankPage";
import { api } from "../../lib/api";

describe("Sample Ranking value column", () => {
  it("names the gene and unit instead of 'Value'", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({
      datasets: [
        { dataset_id: "depmap", display_name: "DepMap", group_columns: [], expression_unit: "log2_tpm" },
        { dataset_id: "target_all_p2", display_name: "TARGET", group_columns: [], expression_unit: "tpm" },
      ],
    } as never);
    vi.spyOn(api, "rankByGene").mockResolvedValue({
      gene: "MYCN",
      n: 1,
      rows: [{ sample_id: "ACH-1", value: 6.4 }],
    } as never);
    render(
      <QueryClientProvider client={new QueryClient()}>
        <RankPage />
      </QueryClientProvider>,
    );
    await waitFor(() => expect(screen.getByLabelText("Dataset")).toHaveValue("depmap"));
    fireEvent.click(screen.getByRole("button", { name: /rank 1 gene/i }));
    expect(await screen.findByRole("columnheader", { name: "MYCN expression (log2(TPM+1))" })).toBeInTheDocument();
  });
});
