/**
 * The blocker this suite exists for.
 *
 * mygene.info returns `alias` as a bare string when a gene has exactly
 * one (LYL1 -> "bHLHa18"). The API declares `aliases: list[str]`, and this
 * component did `aliases.length > 0 && aliases.join(", ")`. A string
 * passes `.length` (8, truthy) and then throws on `.join` -- and because
 * GeneAnnotation renders inside every results page with no error boundary
 * above it, the whole workspace unmounted to a blank page.
 *
 * LYL1 is in this app's own ETP-TF5 preset, so it was squarely in the
 * primary use case. It was found by an agent clicking the UI; no test
 * existed that could have caught it. These are those tests.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { GeneAnnotation } from "../GeneAnnotation";
import { api } from "../../lib/api";

function renderWithQuery(ui: React.ReactElement) {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0 } },
  });
  return render(<QueryClientProvider client={client}>{ui}</QueryClientProvider>);
}

const BASE = {
  symbol: "LYL1",
  name: "LYL1 basic helix-loop-helix family member",
  summary: null,
  ensembl_gene_id: "ENSG00000104903",
  entrez_gene_id: "4066",
};

beforeEach(() => {
  vi.restoreAllMocks();
});

describe("GeneAnnotation aliases handling", () => {
  it("renders a single alias delivered as a bare string without crashing", async () => {
    // The exact upstream shape that white-screened the app.
    vi.spyOn(api, "geneInfo").mockResolvedValue({
      ...BASE,
      aliases: "bHLHa18" as unknown as string[],
    });

    renderWithQuery(<GeneAnnotation gene="LYL1" />);

    expect(await screen.findByText(/bHLHa18/)).toBeInTheDocument();
    // Still renders the gene itself -- i.e. the component survived.
    expect(screen.getByText("LYL1")).toBeInTheDocument();
  });

  it("renders a normal list of aliases joined for display", async () => {
    vi.spyOn(api, "geneInfo").mockResolvedValue({
      ...BASE,
      symbol: "MEF2C",
      aliases: ["C5DELq14.3", "DEL5q14.3", "NEDHSIL"],
    });

    renderWithQuery(<GeneAnnotation gene="MEF2C" />);

    expect(await screen.findByText(/C5DELq14\.3, DEL5q14\.3, NEDHSIL/)).toBeInTheDocument();
  });

  it("omits the aliases row entirely when there are none", async () => {
    vi.spyOn(api, "geneInfo").mockResolvedValue({ ...BASE, symbol: "CD34", aliases: [] });

    renderWithQuery(<GeneAnnotation gene="CD34" />);

    await screen.findByText("CD34");
    expect(screen.queryByText(/Aliases/i)).not.toBeInTheDocument();
  });

  it.each([
    ["null", null],
    ["undefined", undefined],
    ["a number", 42],
    ["an object", { nope: true }],
  ])("survives a malformed aliases field: %s", async (_label, aliases) => {
    // The component must not be the thing that takes the page down, no
    // matter what the upstream sends. Asserting the render simply
    // completes is the point.
    vi.spyOn(api, "geneInfo").mockResolvedValue({
      ...BASE,
      aliases: aliases as unknown as string[],
    });

    renderWithQuery(<GeneAnnotation gene="LYL1" />);

    expect(await screen.findByText("LYL1")).toBeInTheDocument();
  });

  it("fails quiet when annotation is unavailable, rather than showing an error", async () => {
    // A gene with no external annotation isn't an error for the page --
    // the comparison itself doesn't depend on this panel.
    vi.spyOn(api, "geneInfo").mockRejectedValue(new Error("404"));

    const { container } = renderWithQuery(<GeneAnnotation gene="NOSUCHGENE" />);

    await waitFor(() => expect(container).toBeEmptyDOMElement());
  });
});
