/**
 * C1orf112 is now FIRRM. On GSE324984 (current symbols) a search for
 * C1orf112 404'd with a bare "not found" although the gene is there as
 * FIRRM. The annotation panel already knew the current symbol; it now
 * says so and offers the switch -- it never changes the search silently.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";

import { GeneAnnotation } from "../GeneAnnotation";
import { api } from "../../lib/api";

function renderWith(gene: string, symbol: string, onUseSymbol?: (s: string) => void) {
  vi.spyOn(api, "geneInfo").mockResolvedValue({
    symbol, name: "x", summary: null, aliases: [], ensembl_gene_id: null, entrez_gene_id: null,
  } as unknown as Awaited<ReturnType<typeof api.geneInfo>>);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <GeneAnnotation gene={gene} onUseSymbol={onUseSymbol} />
    </QueryClientProvider>,
  );
}

describe("GeneAnnotation renamed-symbol notice", () => {
  it("names the current symbol and offers it, without switching on its own", async () => {
    const onUse = vi.fn();
    renderWith("C1ORF112", "FIRRM", onUse);
    const note = await screen.findByTestId("renamed-symbol");
    expect(note.textContent).toContain("current symbol is FIRRM");
    expect(onUse).not.toHaveBeenCalled();
    await userEvent.click(screen.getByRole("button", { name: /search firrm instead/i }));
    expect(onUse).toHaveBeenCalledWith("FIRRM");
  });

  it("stays quiet when the searched name is already current", async () => {
    renderWith("mycn", "MYCN");
    await screen.findByText("MYCN");
    expect(screen.queryByTestId("renamed-symbol")).toBeNull();
  });
});
