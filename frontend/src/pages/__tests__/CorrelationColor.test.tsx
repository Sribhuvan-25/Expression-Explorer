/**
 * GSE66586 MYCN vs PPP2R5B: r = -0.85 over 10 lines, mostly because the 5
 * MYCN-amplified lines sit apart from the rest. The scatter used to be one
 * colour with no units and no hint that r pools the groups.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";

import { CorrelationTab } from "../CorrelationPage";
import { api } from "../../lib/api";

const groups = ["non-malignant", "MYCN-amplified", "MYCN-amplified", "MYCN non-amplified"];

function setup() {
  vi.spyOn(api, "listDatasets").mockResolvedValue({
    datasets: [
      {
        dataset_id: "gse66586",
        display_name: "GSE66586",
        group_columns: ["mycn_status", "cell_line"],
        expression_unit: "log2_intensity",
      },
    ],
  } as unknown as Awaited<ReturnType<typeof api.listDatasets>>);
  vi.spyOn(api, "correlation").mockResolvedValue({
    gene_a: "MYCN",
    gene_b: "PPP2R5B",
    method: "pearson",
    n: 4,
    coefficient: -0.85,
    p_value: 0.002,
    points: groups.map((g, i) => ({
      sample_id: `S${i}`,
      x: 7 + i * 2,
      y: 10 - i * 0.6,
      group_columns: { mycn_status: g, cell_line: `L${i}` },
    })),
  } as unknown as Awaited<ReturnType<typeof api.correlation>>);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <CorrelationTab datasetId="gse66586" datasetName="GSE66586" />
    </QueryClientProvider>,
  );
}

describe("correlation scatter shows what r pools", () => {
  it("colours by the first grouping by default, names units, and says r spans the groups", async () => {
    const { container } = setup();
    expect(await screen.findByDisplayValue("mycn_status")).toBeInTheDocument();
    await userEvent.click(screen.getByRole("button", { name: /run correlation/i }));
    expect(await screen.findByTestId("pools-groups-note")).toHaveTextContent(
      "r is computed across all 4 samples together, spanning 3 mycn_status groups",
    );
    const fills = new Set([...container.querySelectorAll("svg circle")].map((c) => c.getAttribute("fill")));
    expect(fills.size).toBe(3);
    expect(screen.getByText("MYCN (LOG2 INTENSITY)")).toBeInTheDocument();
  });

  it("drops the note when colouring is turned off", async () => {
    setup();
    await userEvent.selectOptions(await screen.findByDisplayValue("mycn_status"), "");
    await userEvent.click(screen.getByRole("button", { name: /run correlation/i }));
    await screen.findByText("-0.850");
    expect(screen.queryByTestId("pools-groups-note")).toBeNull();
  });
});

describe("correlation axes name the scale the points are on", () => {
  it("labels TPM data as log2(TPM+1) when the response says so", async () => {
    setup();
    vi.mocked(api.correlation).mockResolvedValue({
      ...(await api.correlation("gse66586", "MYCN", "PPP2R5B", "pearson")),
      value_unit: "log2_tpm",
    } as unknown as Awaited<ReturnType<typeof api.correlation>>);
    await screen.findByDisplayValue("mycn_status");
    await userEvent.click(screen.getByRole("button", { name: /run correlation/i }));
    expect(await screen.findByText("MYCN (LOG2(TPM+1))")).toBeInTheDocument();
  });
});

describe("correlation PNG names the dataset", () => {
  it("so the same pair from two datasets does not overwrite", async () => {
    const exportFigure = await import("../../lib/exportFigure");
    const spy = vi.spyOn(exportFigure, "exportSvgAsPng").mockResolvedValue(undefined);
    setup();
    await screen.findByDisplayValue("mycn_status");
    await userEvent.click(screen.getByRole("button", { name: /run correlation/i }));
    await userEvent.click(await screen.findByRole("button", { name: /png/i }));
    expect(spy.mock.calls.at(-1)![1]).toBe("gse66586-MYCN-vs-PPP2R5B-correlation");
  });
});
