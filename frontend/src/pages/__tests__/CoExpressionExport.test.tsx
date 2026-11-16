/**
 * The co-expression CSV was named "<gene>-<direction>-network.csv" with no
 * dataset anywhere: one gene's network exported from two datasets
 * overwrote one file, and the file did not say where it came from
 * (METHODS.md 5.7).
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";

import { NetworkTab } from "../CorrelationPage";
import * as exportTable from "../../lib/exportTable";
import { api } from "../../lib/api";

describe("co-expression CSV identity", () => {
  it("names the dataset in the file and in every row", async () => {
    const spy = vi.spyOn(exportTable, "exportRowsAsCsv").mockImplementation(() => {});
    vi.spyOn(api, "coExpression").mockResolvedValue({
      gene: "MYCN", method: "pearson", direction: "positive",
      network: [{ feature_id: "PPA1", symbol: "PPA1", coefficient: 0.97 }],
    } as unknown as Awaited<ReturnType<typeof api.coExpression>>);
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <NetworkTab datasetId="gse66586" />
      </QueryClientProvider>,
    );
    await userEvent.click(screen.getByRole("button", { name: /build network/i }));
    await userEvent.click(await screen.findByTitle("Download as CSV"));
    const [filename, , rows] = spy.mock.calls.at(-1)!;
    expect(filename).toContain("gse66586");
    expect((rows as Record<string, unknown>[])[0].dataset_id).toBe("gse66586");
  });
});
