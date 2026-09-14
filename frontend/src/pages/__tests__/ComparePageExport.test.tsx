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
