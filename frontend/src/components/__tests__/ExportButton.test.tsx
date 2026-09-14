/**
 * Exported figures must identify their subject.
 *
 * Expression Compare wrote every PNG to "<dataset>-by-group.png" with the
 * dataset name as the only title, so exporting two genes silently
 * overwrote the first file, and the gene appeared nowhere inside the
 * image either -- a saved figure could not be identified afterwards.
 * See METHODS.md 5.7.
 *
 * This drives the real ExportButton and asserts on what it actually
 * hands the exporter, so it fails if the call site stops passing the
 * gene. Asserting against a re-implementation of the naming rule would
 * pass forever regardless of what the component does.
 */
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { useRef } from "react";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { ExportButton } from "../ExportButton";
import * as exportFigure from "../../lib/exportFigure";

/** Mirrors how ComparePage wires the button, including its own ref. */
function Harness({ gene, datasetId, displayName }: { gene: string; datasetId: string; displayName: string }) {
  const chartRef = useRef<SVGSVGElement>(null);
  return (
    <>
      <svg ref={chartRef} />
      <ExportButton
        svgRef={chartRef}
        filename={`${gene}-${datasetId}-by-group`}
        title={`${gene} — ${displayName}`}
        subtitle={`n = 10 · Mann–Whitney with FDR correction`}
        statLines={[]}
      />
    </>
  );
}

let spy: ReturnType<typeof vi.spyOn>;

beforeEach(() => {
  vi.restoreAllMocks();
  spy = vi.spyOn(exportFigure, "exportSvgAsPng").mockResolvedValue(undefined);
});

async function exportOnce(props: { gene: string; datasetId: string; displayName: string }) {
  const { unmount } = render(<Harness {...props} />);
  await userEvent.click(screen.getByRole("button", { name: /png/i }));
  await waitFor(() => expect(spy).toHaveBeenCalled());
  const [, filename, options] = spy.mock.calls.at(-1)!;
  unmount();
  return { filename: filename as string, options: options as { title: string } };
}

describe("ExportButton figure identity", () => {
  it("puts the gene in the filename so two genes do not collide", async () => {
    const a = await exportOnce({ gene: "MEF2C", datasetId: "depmap", displayName: "DepMap" });
    spy.mockClear();
    const b = await exportOnce({ gene: "DTX1", datasetId: "depmap", displayName: "DepMap" });

    expect(a.filename).toBe("MEF2C-depmap-by-group");
    expect(b.filename).toBe("DTX1-depmap-by-group");
    expect(a.filename).not.toBe(b.filename);
  });

  it("puts the gene in the in-image title, not just the filename", async () => {
    // A figure saved to disk has only its title to identify it by.
    const { options } = await exportOnce({
      gene: "MEF2C",
      datasetId: "depmap",
      displayName: "DepMap cell line RNA-seq",
    });
    expect(options.title).toContain("MEF2C");
    expect(options.title).toContain("DepMap cell line RNA-seq");
  });

  it("distinguishes one gene across two datasets", async () => {
    const a = await exportOnce({ gene: "MEF2C", datasetId: "depmap", displayName: "DepMap" });
    spy.mockClear();
    const b = await exportOnce({ gene: "MEF2C", datasetId: "gds4299", displayName: "GDS4299" });
    expect(a.filename).not.toBe(b.filename);
  });

  it("does not attempt an export when there is no chart to export", async () => {
    function Empty() {
      const ref = useRef<SVGSVGElement>(null);
      return <ExportButton svgRef={ref} filename="x" title="x" statLines={[]} />;
    }
    render(<Empty />);
    await userEvent.click(screen.getByRole("button", { name: /png/i }));
    expect(spy).not.toHaveBeenCalled();
  });
});
