/**
 * The context bar's dataset must actually reach a pane. Under StrictMode
 * (how the app runs in development) the dataset hook's "pick the first
 * dataset" fallback ran a second time with a stale closure and replaced
 * the adopted context dataset: Correlation sat on DepMap while the
 * context said GSE66586, with only a "use GSE66586" link to show for it.
 */
import { act, render, screen } from "@testing-library/react";
import { StrictMode } from "react";
import { afterEach, describe, expect, it } from "vitest";

import { useLinkedDataset, useLinkedGene, useWorkbench, WorkbenchProvider } from "../workbench";

const IDS = ["depmap", "target_all_p2", "gse66586"];

function DatasetProbe() {
  const [datasetId] = useLinkedDataset(IDS);
  return <p data-testid="dataset">{datasetId}</p>;
}

let setContext: ReturnType<typeof useWorkbench> extends infer T ? (T extends { setContext: infer S } ? S : never) : never;
function GeneProbe() {
  const ctx = useWorkbench();
  if (ctx) setContext = ctx.setContext;
  const [gene, setGene, link] = useLinkedGene("MYCN");
  return (
    <>
      <p data-testid="gene">{gene}</p>
      <p data-testid="linked">{String(link.linked)}</p>
      <button onClick={() => setGene("MYC")}>override</button>
    </>
  );
}

afterEach(() => localStorage.clear());

describe("workbench context", () => {
  it("a pane adopts the context dataset under StrictMode", () => {
    localStorage.setItem("expression-explorer:workbench-context", JSON.stringify({ genes: [], datasetId: "gse66586" }));
    render(
      <StrictMode>
        <WorkbenchProvider>
          <DatasetProbe />
        </WorkbenchProvider>
      </StrictMode>,
    );
    expect(screen.getByTestId("dataset")).toHaveTextContent("gse66586");
  });

  it("falls back to the pane's own default when the context dataset is not offered", () => {
    localStorage.setItem("expression-explorer:workbench-context", JSON.stringify({ genes: [], datasetId: "gds4299" }));
    render(
      <StrictMode>
        <WorkbenchProvider>
          <DatasetProbe />
        </WorkbenchProvider>
      </StrictMode>,
    );
    expect(screen.getByTestId("dataset")).toHaveTextContent("depmap");
  });

  it("keeps a pane's override until the context gene changes again", () => {
    render(
      <StrictMode>
        <WorkbenchProvider>
          <GeneProbe />
        </WorkbenchProvider>
      </StrictMode>,
    );
    expect(screen.getByTestId("gene")).toHaveTextContent("MYCN");
    act(() => setContext({ genes: ["LMO2"] }));
    expect(screen.getByTestId("gene")).toHaveTextContent("LMO2");
    expect(screen.getByTestId("linked")).toHaveTextContent("true");
    act(() => screen.getByRole("button", { name: "override" }).click());
    expect(screen.getByTestId("gene")).toHaveTextContent("MYC");
    expect(screen.getByTestId("linked")).toHaveTextContent("false");
    act(() => setContext({ genes: ["HHEX"] }));
    expect(screen.getByTestId("gene")).toHaveTextContent("HHEX");
  });
});
