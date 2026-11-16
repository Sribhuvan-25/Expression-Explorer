/**
 * Guards the export-identity fix at its real call site.
 *
 * Expression Compare wrote every PNG to "<dataset>-by-group.png" with the
 * dataset name as the only title, so exporting two genes silently
 * overwrote the first file and the gene appeared nowhere in the image.
 * See METHODS.md 5.7.
 *
 * Testing ExportButton alone is not enough -- it only proves the button
 * forwards whatever props it gets. The defect was the *call site* not
 * passing the gene, so this renders the real DatasetPanel and asserts on
 * what reaches the exporter.
 */
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { DatasetPanel } from "../ComparePage";
import * as exportFigure from "../../lib/exportFigure";
import type { CompareMultiDatasetResult } from "../../lib/api";

function result(datasetId: string, displayName: string): CompareMultiDatasetResult {
  return {
    dataset_id: datasetId,
    display_name: displayName,
    skipped: false,
    assay_type: "rna_seq",
    expression_unit: "tpm",
    n_dataset_total: 6,
    n_excluded: 0,
    exclusion_reason: null,
    points: [
      { sample_id: "S1", group: "ETP", value: 5 },
      { sample_id: "S2", group: "ETP", value: 6 },
      { sample_id: "S3", group: "ETP", value: 7 },
      { sample_id: "S4", group: "non-ETP", value: 1 },
      { sample_id: "S5", group: "non-ETP", value: 2 },
      { sample_id: "S6", group: "non-ETP", value: 3 },
    ],
    pairwise_tests: [],
    kruskal_wallis: { h_stat: 1, p_value: 0.5 },
  };
}

let spy: ReturnType<typeof vi.spyOn>;

beforeEach(() => {
  vi.restoreAllMocks();
  spy = vi.spyOn(exportFigure, "exportSvgAsPng").mockResolvedValue(undefined);
});

async function exportFrom(gene: string, datasetId: string, displayName: string) {
  const { unmount } = render(<DatasetPanel result={result(datasetId, displayName)} gene={gene} />);
  await userEvent.click(screen.getByRole("button", { name: /png/i }));
  await waitFor(() => expect(spy).toHaveBeenCalled());
  const [, filename, options] = spy.mock.calls.at(-1)!;
  unmount();
  return { filename: filename as string, title: (options as { title: string }).title };
}

describe("Expression Compare export identity", () => {
  it("exporting two genes from one dataset produces two distinct files", async () => {
    // The exact defect: the second export silently overwrote the first.
    const a = await exportFrom("MEF2C", "depmap", "DepMap cell line RNA-seq");
    spy.mockClear();
    const b = await exportFrom("DTX1", "depmap", "DepMap cell line RNA-seq");

    expect(a.filename).not.toBe(b.filename);
    expect(a.filename).toContain("MEF2C");
    expect(b.filename).toContain("DTX1");
  });

  it("names the gene inside the image too, not only in the filename", async () => {
    const { title } = await exportFrom("MEF2C", "depmap", "DepMap cell line RNA-seq");
    expect(title).toContain("MEF2C");
    expect(title).toContain("DepMap cell line RNA-seq");
  });

  it("keeps one gene distinguishable across two datasets", async () => {
    const a = await exportFrom("MEF2C", "depmap", "DepMap");
    spy.mockClear();
    const b = await exportFrom("MEF2C", "gds4299", "GDS4299");
    expect(a.filename).not.toBe(b.filename);
  });
});

describe("Expression Compare export carries the small-sample floor", () => {
  it("states the smallest possible p in the image when it exceeds 0.05", async () => {
    // The exported PNG travels without the on-page notice: "p = 0.0952"
    // alone in a figure reads as a negative result (GSE66586, 5 vs 2).
    const r = result("gse66586", "GSE66586");
    r.pairwise_tests = [
      { group_a: "A", group_b: "B", n_a: 5, n_b: 2, u_stat: 10, p_value: 0.0952, min_p_achievable: 2 / 21 },
      { group_a: "C", group_b: "D", n_a: 30, n_b: 40, u_stat: 10, p_value: 0.01, min_p_achievable: null },
    ] as CompareMultiDatasetResult["pairwise_tests"];
    render(<DatasetPanel result={r} gene="PPP2R5B" />);
    await userEvent.click(screen.getByRole("button", { name: /png/i }));
    await waitFor(() => expect(spy).toHaveBeenCalled());
    const lines = (spy.mock.calls.at(-1)![2] as { statLines: string[] }).statLines;
    expect(lines.find((l) => l.startsWith("A vs B"))).toContain("smallest possible p at n = 5 vs 2: 0.095");
    expect(lines.find((l) => l.startsWith("C vs D"))).not.toContain("smallest possible");
  });
});

describe("Expression Compare export prints the FDR q it claims", () => {
  it("puts each pair's FDR q beside its raw p, as the on-page table does", async () => {
    // The subtitle reads "Mann–Whitney with FDR correction" but the image
    // printed only raw p, so a reader of the PNG saw p = 0.0357 for a pair
    // whose corrected q is 0.107.
    const r = result("gse66586", "GSE66586");
    r.pairwise_tests = [
      { group_a: "MYCN non-amplified", group_b: "MYCN-amplified", n_a: 3, n_b: 5, u_stat: 15, p_value: 2 / 56, q_value: 0.1071 },
    ] as CompareMultiDatasetResult["pairwise_tests"];
    render(<DatasetPanel result={r} gene="PPP2R5B" />);
    await userEvent.click(screen.getByRole("button", { name: /png/i }));
    await waitFor(() => expect(spy).toHaveBeenCalled());
    const lines = (spy.mock.calls.at(-1)![2] as { statLines: string[] }).statLines;
    expect(lines.find((l) => l.startsWith("MYCN non-amplified"))).toBe(
      "MYCN non-amplified vs MYCN-amplified: p = 0.0357, FDR q = 0.1071",
    );
  });
});
