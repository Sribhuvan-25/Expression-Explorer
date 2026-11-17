/**
 * The CRISPR result was titled "MYB vs MYB" -- MYB expression against
 * MYB's CRISPR gene effect -- on screen and in the exported PNG, and the
 * figure had no caption. Drives the real Multi-omics page and asserts on
 * what it hands the exporter.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { MultiOmicsPage } from "../MultiOmicsPage";
import { api } from "../../lib/api";
import * as exportFigure from "../../lib/exportFigure";

const CRISPR_DESCRIPTION =
  "CRISPR knockout gene-effect score. 0 = no effect; -1 is about the median of known common-essential genes, so more negative means the line depends more strongly on that gene.";

beforeEach(() => {
  vi.restoreAllMocks();
  vi.spyOn(api, "listDatasets").mockResolvedValue({
    datasets: [{ dataset_id: "depmap", display_name: "DepMap", group_columns: [], expression_unit: "log2_tpm" }],
  } as unknown as Awaited<ReturnType<typeof api.listDatasets>>);
  vi.spyOn(api, "auxLayers").mockResolvedValue({
    layers: [
      {
        layer: "crispr_gene_effect",
        value_label: "Gene effect (Chronos)",
        value_description: CRISPR_DESCRIPTION,
        source_note: "",
        n_features: 10,
        n_usable_features: 10,
        n_samples_covered: 4,
        n_dataset_total: 6,
      },
    ],
  } as unknown as Awaited<ReturnType<typeof api.auxLayers>>);
  vi.spyOn(api, "expressionVsAux").mockResolvedValue({
    gene: "MYB",
    expression_value_unit: "log2_tpm",
    aux_feature: "MYB",
    aux_feature_label: "MYB",
    aux_feature_moa: null,
    layer: "crispr_gene_effect",
    value_label: "Gene effect (Chronos)",
    value_description: CRISPR_DESCRIPTION,
    method: "pearson",
    n: 4,
    n_dataset_total: 6,
    n_excluded: 2,
    n_missing_layer: 2,
    n_missing_feature: 0,
    exclusion_reason: "2 have no crispr gene effect data",
    coefficient: -0.652,
    p_value: 2.5e-12,
    points: [1, 2, 3, 4].map((v) => ({ sample_id: `L${v}`, x: v, y: -v / 2 })),
  } as unknown as Awaited<ReturnType<typeof api.expressionVsAux>>);
});

async function runCrispr() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  const utils = render(
    <QueryClientProvider client={qc}>
      <MultiOmicsPage />
    </QueryClientProvider>,
  );
  // The page lands on the CRISPR tab once the layer list arrives.
  await userEvent.click(await screen.findByRole("button", { name: /^correlate$/i }, { timeout: 5000 }));
  await screen.findByText("-0.652", undefined, { timeout: 5000 });
  return utils;
}

describe("expression vs CRISPR names both sides", () => {
  it("titles the result with what each axis measures, not 'MYB vs MYB'", async () => {
    await runCrispr();
    expect(screen.queryByText("MYB vs MYB")).toBeNull();
    expect(screen.getByText("MYB expression vs MYB CRISPR gene effect")).toBeInTheDocument();
    expect(screen.getByText("MYB CRISPR gene effect (Chronos)")).toBeInTheDocument();
  });

  it("exports the same title, an unambiguous file name, and the on-screen caption", async () => {
    const spy = vi.spyOn(exportFigure, "exportSvgAsPng").mockResolvedValue(undefined);
    const { container } = await runCrispr();
    await userEvent.click(screen.getByRole("button", { name: /png/i }));
    await waitFor(() => expect(spy).toHaveBeenCalled());
    const [, filename, options] = spy.mock.calls.at(-1)!;
    expect(filename).toBe("depmap-MYB-expression-vs-MYB-crispr-gene-effect-pearson");
    expect(options.title).toBe("MYB expression vs MYB CRISPR gene effect");

    const caption = container.querySelector("[data-figure-caption]");
    expect(caption).not.toBeNull();
    expect(options.caption).toBe(caption!.textContent);
    // Says what was computed, on what, and who is missing and why.
    expect(options.caption).toContain("Each point is one of 4 samples: MYB expression (log2(TPM+1)) (x)");
    expect(options.caption).toContain("CRISPR knockout gene effect of MYB (DepMap Chronos)");
    expect(options.caption).toContain("2 of 6 samples excluded: 2 have no crispr gene effect data.");
  });
});
