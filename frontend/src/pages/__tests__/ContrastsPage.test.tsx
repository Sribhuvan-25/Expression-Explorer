/**
 * Guards the Perturbation Contrasts pane against the defects that matter
 * for the MYCL / PP2A question on GSE324984:
 *  - an untested gene (all-null stats; MYCN and PPP2R2B in C4-2B) must read
 *    "not tested" with its reason, never blank, 0, or a crash;
 *  - the "Across both experiments" call must follow sign and q exactly;
 *  - the preset must send current symbols (CIP2A, PABIR1) -- legacy names
 *    are aliases in some datasets and 404 in others -- and keep the
 *    researcher's categories;
 *  - a gene whose TPM means move opposite to the count-based log2FC must
 *    be flagged, or the table shows a contradictory pair with no warning.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";

import { ContrastsPage, classifyAcross, ACROSS_LABELS } from "../ContrastsPage";
import { api } from "../../lib/api";
import type { ContrastPanelCell, ContrastPanelResult, ContrastPanelRow, ContrastSummary, DatasetSummary } from "../../lib/api";

const DATASET: DatasetSummary = {
  dataset_id: "gse324984",
  display_name: "GSE324984",
  group_columns: ["arm"],
  supports_survival: false,
  n_samples: 8,
  assay_type: "rna_seq",
  expression_unit: "tpm",
  n_contrasts: 2,
  notes: "TPM file is truncated.",
};

function contrast(id: string, stratum: string, note: string): ContrastSummary {
  return {
    contrast_id: id,
    label: `${stratum} contrast`,
    group_column: "arm",
    reference: `${stratum} ctrl`,
    treatment: `${stratum} trt`,
    stratum_column: "cell_line",
    stratum_value: stratum,
    method: "limma-trend moderated t",
    direction_note: note,
    n_reference: 2,
    n_treatment: 2,
    reference_samples: ["R1", "R2"],
    treatment_samples: ["T1", "T2"],
    n_genes_input: 100,
    n_genes_tested: 80,
    n_dataset_total: 8,
    n_excluded: 4,
    exclusion_reason: `By design: only the ${stratum} samples enter this contrast.`,
  };
}
const CONTRASTS = [
  contrast("oe", "C4-2B", "log2FC > 0: higher with MYCL overexpressed"),
  contrast("kd", "NCI-H660", "log2FC > 0: higher with MYCL knocked down"),
];

function cell(id: string, over: Partial<ContrastPanelCell> = {}): ContrastPanelCell {
  return {
    contrast_id: id,
    in_contrast_table: true,
    tested: true,
    not_tested_reason: null,
    log2_fc: 1,
    ave_log_cpm: 5,
    t: 3,
    p_value: 0.001,
    q_value: 0.01,
    display_direction_agrees: true,
    reference: { label: "ref", mean_expression: 10, values: [{ sample_id: "R1", value: 9 }, { sample_id: "R2", value: 11 }] },
    treatment: { label: "trt", mean_expression: 20, values: [{ sample_id: "T1", value: 19 }, { sample_id: "T2", value: 21 }] },
    ...over,
  };
}
function row(gene: string, cells: ContrastPanelCell[]): ContrastPanelRow {
  return { gene, symbol: gene, feature_id: `ENSG_${gene}`, aliases: [], alias_note: null, contrasts: cells };
}

const UNTESTED_REASON = "Not tested: too few reads in this experiment's samples (edgeR filterByExpr).";

function mockApi(rows: (genes: string[]) => ContrastPanelRow[]) {
  vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [DATASET] });
  vi.spyOn(api, "contrasts").mockResolvedValue({ dataset_id: "gse324984", expression_unit: "tpm", contrasts: CONTRASTS });
  vi.spyOn(api, "contrastTop").mockResolvedValue({
    expression_unit: "tpm",
    contrast: CONTRASTS[0],
    direction: "both",
    n_significant_q05: 0,
    genes: [],
  });
  return vi.spyOn(api, "contrastPanel").mockImplementation(
    async (_id: string, genes: string[]): Promise<ContrastPanelResult> => ({
      expression_unit: "tpm",
      contrasts: CONTRASTS,
      rows: rows(genes),
      missing: [],
    }),
  );
}

function renderPage() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false, gcTime: 0 } } });
  return render(
    <QueryClientProvider client={client}>
      <ContrastsPage />
    </QueryClientProvider>,
  );
}

async function runWith(text: string | null) {
  const box = await screen.findByLabelText(/gene symbols/i);
  if (text !== null) {
    await userEvent.clear(box);
    await userEvent.type(box, text);
  } else {
    await userEvent.click(screen.getByRole("button", { name: /pp2a signalling panel/i }));
  }
  await waitFor(() => expect(screen.getByRole("button", { name: /run panel/i })).toBeEnabled());
  await userEvent.click(screen.getByRole("button", { name: /run panel/i }));
}

describe("Perturbation Contrasts gene panel", () => {
  it("renders an untested gene (all-null stats) as 'not tested' with its reason", async () => {
    mockApi(() => [
      row("MYCN", [
        cell("oe", {
          tested: false,
          not_tested_reason: UNTESTED_REASON,
          log2_fc: null,
          ave_log_cpm: null,
          t: null,
          p_value: null,
          q_value: null,
          display_direction_agrees: null,
          reference: { label: "ref", mean_expression: 0, values: [{ sample_id: "R1", value: 0 }, { sample_id: "R2", value: null }] },
        }),
        cell("kd", { log2_fc: -0.08, q_value: 0.9 }),
      ]),
    ]);
    renderPage();
    await runWith("MYCN");

    const geneRow = await screen.findByTestId("gene-row");
    const notTested = within(geneRow).getByText(/not tested/i);
    expect(notTested.closest("td")).toHaveAttribute("title", UNTESTED_REASON);
    // Reason is visible on the page, not only in a tooltip.
    expect(screen.getAllByText(new RegExp(UNTESTED_REASON.slice(0, 30))).length).toBeGreaterThan(0);
    // The other experiment still shows its real number.
    expect(within(geneRow).getByText(/0\.08/)).toBeInTheDocument();
  });

  it("preset sends current symbols and groups rows under the researcher's categories", async () => {
    const spy = mockApi((genes) => genes.map((g) => row(g, [cell("oe"), cell("kd")])));
    renderPage();
    await runWith(null);

    await waitFor(() => expect(spy).toHaveBeenCalled());
    const sent = spy.mock.calls.at(-1)![1] as string[];
    expect(sent).toContain("CIP2A");
    expect(sent).toContain("PABIR1");
    expect(sent).not.toContain("KIAA1524");
    expect(sent).not.toContain("FAM122A");
    expect(sent).toHaveLength(23);

    const headings = (await screen.findAllByTestId("category-heading")).map((h) => h.textContent);
    expect(headings).toEqual(["MYC family", "PP2A inhibitors", "PP2A regulators", "PP2A subunits and isoforms"]);
    // Rows follow their heading: first row after "PP2A inhibitors" is CIP2A.
    const tbodyRows = Array.from(document.querySelectorAll("tbody tr"));
    const inhib = tbodyRows.findIndex((r) => r.textContent === "PP2A inhibitors");
    expect(tbodyRows[inhib + 1].textContent).toMatch(/^CIP2A/);
    expect(tbodyRows[1].textContent).toMatch(/^MYC(?!L|N)/);
  });

  it("flags a gene whose TPM direction disagrees with the count-based log2FC, and only that one", async () => {
    // HSPA12A on GSE324984: log2FC +1.04 while TPM went 4.8 -> 3.6.
    const cpmArm = (m: number) => ({
      label: "x",
      mean_expression: 1,
      values: [{ sample_id: "S", value: 1 }],
      mean_norm_cpm: m,
      norm_cpm: [{ sample_id: "S", value: m }],
    });
    mockApi(() => [
      row("HSPA12A", [
        cell("oe", { log2_fc: 1.04, display_direction_agrees: false, reference: cpmArm(5.49), treatment: cpmArm(11.25) }),
        cell("kd"),
      ]),
      row("SET", [cell("oe"), cell("kd")]),
      // null = cannot judge (untested): must not warn.
      row("MYCN", [cell("oe", { tested: false, log2_fc: null, q_value: null, display_direction_agrees: null }), cell("kd")]),
    ]);
    renderPage();
    await runWith("HSPA12A SET MYCN");
    const rows = await screen.findAllByTestId("gene-row");
    const flag = within(rows[0]).getByTestId("direction-disagrees");
    expect(flag.textContent).toMatch(/TPM change points the other way from the count-based fold-change/);
    // The normalised counts the test used are shown beside the TPM.
    expect(within(rows[0]).getByTestId("norm-cpm").textContent).toMatch(/5\.49\s*→\s*11\.25/);
    expect(within(rows[1]).queryByTestId("direction-disagrees")).not.toBeInTheDocument();
    expect(within(rows[2]).queryByTestId("direction-disagrees")).not.toBeInTheDocument();
    expect(screen.getByTestId("direction-disagrees-legend").textContent).toMatch(/normalised counts/);
  });

  it("flags the same disagreement in the Top genes table", async () => {
    mockApi(() => []);
    vi.spyOn(api, "contrastTop").mockResolvedValue({
      expression_unit: "tpm",
      contrast: CONTRASTS[0],
      direction: "both",
      n_significant_q05: 2,
      genes: [
        { feature_id: "E1", symbol: "HSPA12A", log2_fc: 1.04, ave_log_cpm: 3, p_value: 1e-4, q_value: 0.01, in_expression_matrix: true, mean_reference: 4.8, mean_treatment: 3.6, mean_norm_cpm_reference: 5.5, mean_norm_cpm_treatment: 11.3, display_direction_agrees: false },
        { feature_id: "E2", symbol: "SET", log2_fc: -0.45, ave_log_cpm: 9, p_value: 2e-4, q_value: 0.02, in_expression_matrix: true, mean_reference: 583, mean_treatment: 372, mean_norm_cpm_reference: 865, mean_norm_cpm_treatment: 634, display_direction_agrees: null },
      ],
    });
    renderPage();
    await userEvent.click(await screen.findByRole("button", { name: /top genes/i }));
    const flags = await screen.findAllByTestId("direction-disagrees");
    expect(flags).toHaveLength(1);
    expect(flags[0].closest("tr")!.textContent).toContain("HSPA12A");
    expect(screen.getByTestId("top-direction-legend")).toBeInTheDocument();
  });
});

describe("Across both experiments classification", () => {
  it("covers reciprocal, opposite-not-significant, same, and untested", () => {
    expect(classifyAcross([cell("oe", { log2_fc: 2, q_value: 0.01 }), cell("kd", { log2_fc: -1, q_value: 0.02 })])).toBe(
      "reciprocal",
    );
    // CIP2A: OE -0.94 q=0.0002, KD +0.008 q=0.97 -> opposite but not both significant.
    expect(
      classifyAcross([cell("oe", { log2_fc: -0.94, q_value: 0.0002 }), cell("kd", { log2_fc: 0.008, q_value: 0.97 })]),
    ).toBe("opposite_ns");
    // q exactly 0.05 is not < 0.05.
    expect(classifyAcross([cell("oe", { log2_fc: 2, q_value: 0.05 }), cell("kd", { log2_fc: -1, q_value: 0.01 })])).toBe(
      "opposite_ns",
    );
    expect(classifyAcross([cell("oe", { log2_fc: -0.18 }), cell("kd", { log2_fc: -0.12 })])).toBe("same");
    expect(
      classifyAcross([cell("oe", { tested: false, log2_fc: null, q_value: null }), cell("kd", { log2_fc: -1 })]),
    ).toBe("none");
    expect(ACROSS_LABELS.reciprocal).toBe("Reciprocal (both q<0.05)");
    expect(ACROSS_LABELS.none).toBe("—");
  });
});

/**
 * Identity and accounting fields added after review: count-only genes
 * (CAST has no TPM), names shared by several loci (GTPBP6 / ASMTL, PAR
 * X/Y copies) and the per-contrast sample accounting CLAUDE.md requires.
 * Each was correctly computed by the API and invisible without a reader.
 */
describe("identity notes, shared names and sample accounting", () => {
  it("shows identity notes on a count-only gene's row and its sample accounting in the footnote", async () => {
    const countOnly = (id: string) =>
      cell(id, {
        reference: { label: "ref", mean_expression: null, values: [], mean_norm_cpm: 159.5, norm_cpm: [] },
        treatment: { label: "trt", mean_expression: null, values: [], mean_norm_cpm: 184.4, norm_cpm: [] },
        display_direction_agrees: null,
      });
    mockApi(() => [
      {
        ...row("CAST", [countOnly("oe"), countOnly("kd")]),
        gene_name: "CAST",
        in_expression_matrix: false,
        n_loci_with_name: 2,
        identity_notes: ["No tpm value: this gene is absent from the expression matrix.", "The name CAST is shared by 2 loci."],
      },
      row("SET", [cell("oe"), cell("kd")]),
    ]);
    renderPage();
    await runWith("CAST SET");
    const rows = await screen.findAllByTestId("gene-row");
    const notes = within(rows[0]).getAllByTestId("identity-note").map((n) => n.textContent);
    expect(notes).toEqual(["No tpm value: this gene is absent from the expression matrix.", "The name CAST is shared by 2 loci."]);
    expect(within(rows[1]).queryByTestId("identity-note")).not.toBeInTheDocument();
    // Count-only: no TPM shown as such, never as 0 -- and the stats still render.
    expect(within(rows[0]).getAllByText(/no TPM/i).length).toBeGreaterThan(0);

    const acct = screen.getAllByTestId("sample-accounting").map((a) => a.textContent);
    expect(acct[0]).toMatch(/4 of 8/);
    expect(acct[0]).toMatch(/By design: only the C4-2B samples/);
  });

  it("tags shared-name loci in Top genes and shows the name beside an Ensembl id", async () => {
    mockApi(() => []);
    const base = { ave_log_cpm: 1, p_value: 1e-6, q_value: 1e-4, in_expression_matrix: true, mean_reference: 1, mean_treatment: 2, display_direction_agrees: true };
    vi.spyOn(api, "contrastTop").mockResolvedValue({
      expression_unit: "tpm",
      contrast: CONTRASTS[0],
      direction: "both",
      n_significant_q05: 3,
      genes: [
        { ...base, feature_id: "ENSG00000277053", symbol: "ENSG00000277053", gene_name: "GTF2IP1", n_loci_with_name: 2, log2_fc: -5.1 },
        { ...base, feature_id: "ENSG00000116990", symbol: "MYCL", gene_name: "MYCL", n_loci_with_name: 1, log2_fc: 8.2 },
      ],
    });
    renderPage();
    await userEvent.click(await screen.findByRole("button", { name: /top genes/i }));
    const tags = await screen.findAllByTestId("shared-name");
    expect(tags).toHaveLength(1);
    const tr = tags[0].closest("tr")!;
    expect(tr.textContent).toContain("ENSG00000277053");
    expect(tr.textContent).toContain("GTF2IP1");
    // And the footnote carries the contrast's sample accounting here too.
    expect(screen.getAllByTestId("sample-accounting")[0].textContent).toMatch(/4 of 8/);
  });
});
