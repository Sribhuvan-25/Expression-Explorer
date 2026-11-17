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

import { ContrastsPage, describeAcross } from "../ContrastsPage";
import { api } from "../../lib/api";
import * as exportFigure from "../../lib/exportFigure";
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
  const T = ["C4-2B", "NCI-H660"];
  const pair = (fa: number, qa: number | null, fb: number, qb: number | null) =>
    describeAcross([cell("oe", { log2_fc: fa, q_value: qa }), cell("kd", { log2_fc: fb, q_value: qb })], T);

  it("names which experiment reaches q<0.05 instead of an ambiguous 'not significant in both'", () => {
    // Values from the served PP2A panel.
    expect(pair(1.26, 3.1e-4, -0.89, 0.0011)).toBe("Reciprocal (both q<0.05)"); // PPP2R5B
    expect(pair(-0.94, 2.1e-4, 0.01, 0.97)).toBe("Opposite direction; q<0.05 only in C4-2B"); // CIP2A
    expect(pair(-0.08, 0.44, 0.18, 0.19)).toBe("Opposite direction; q≥0.05 in both"); // PPME1
    expect(pair(1.12, 0.016, 0.4, 0.015)).toBe("Same direction (both q<0.05)"); // PPP2R2C
    expect(pair(-0.14, 0.19, -0.46, 0.0057)).toBe("Same direction; q<0.05 only in NCI-H660"); // PPP2R5D
    expect(pair(-0.18, 0.061, -0.12, 0.47)).toBe("Same direction; q≥0.05 in both"); // ARPP19
    // q exactly 0.05 is not < 0.05.
    expect(pair(2, 0.05, -1, 0.01)).toBe("Opposite direction; q<0.05 only in NCI-H660");
  });

  it("says which experiment a half-tested gene was tested in", () => {
    // MYCN: too few reads in C4-2B.
    expect(
      describeAcross([cell("oe", { tested: false, log2_fc: null, q_value: null }), cell("kd", { log2_fc: -0.08, q_value: 0.9 })], T),
    ).toBe("Tested only in NCI-H660");
    expect(
      describeAcross([cell("oe", { tested: false, log2_fc: null }), cell("kd", { tested: false, log2_fc: null })], T),
    ).toBe("—");
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

describe("dataset notes", () => {
  it("show the first sentence and expand on request", async () => {
    const { DatasetNotes } = await import("../ContrastsPage");
    const { render, screen, fireEvent } = await import("@testing-library/react");
    const design = "Two independent experiments, n=2 per arm: C4-2B MYCL overexpression vs empty vector, and NCI-H660 knockdown.";
    render(<DatasetNotes notes={`PMID 42001796. ${design} Long detail follows. And more.`} />);
    const notes = screen.getByTestId("dataset-notes");
    expect(notes).toHaveTextContent(`PMID 42001796. ${design}`);
    expect(notes).not.toHaveTextContent("Long detail follows");
    fireEvent.click(screen.getByRole("button", { name: "More about this dataset" }));
    expect(notes).toHaveTextContent("Long detail follows. And more.");
    expect(screen.getByRole("button", { name: "Less" })).toHaveAttribute("aria-expanded", "true");
  });
});

/**
 * The fold-change chart, from the researcher-committee review: drawn at a
 * fixed 640 px it used ~40% of a 1440 px pane. It now fills the width its
 * pane grants (the test setup stubs clientWidth to 800) and keeps a
 * readable minimum, scrolling sideways below it.
 */
describe("fold-change chart", () => {
  const chartSvg = () => document.querySelector("svg[aria-label='log2 fold-change per contrast']") as SVGSVGElement;

  it("fills the width the pane grants", async () => {
    mockApi(() => [row("SET", [cell("oe"), cell("kd")])]);
    renderPage();
    await runWith("SET");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    await waitFor(() => expect(chartSvg().getAttribute("width")).toBe("800"));
    expect(chartSvg().getAttribute("viewBox")).toMatch(/^0 0 800 /);
  });

  it("does not shrink below its minimum in a narrow pane, and scrolls sideways instead", async () => {
    const { FC_CHART_MIN_W } = await import("../ContrastsPage");
    vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockReturnValue(300);
    mockApi(() => [row("SET", [cell("oe"), cell("kd")])]);
    renderPage();
    await runWith("SET");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    expect(chartSvg().getAttribute("width")).toBe(String(FC_CHART_MIN_W));
    expect(chartSvg().parentElement!.className).toMatch(/overflow-x-auto/);
  });

  /**
   * Shipped defect: the axis was capped at ±4, so MYCL's +8.24 in C4-2B
   * (served value) was drawn as a bar of length 4 with an arrow -- in an
   * exported PNG, a bar half its value. Every bar must be drawn to its
   * value: bar lengths proportional to |log2FC|, inside the figure.
   */
  it("draws every bar to its true length -- no bar is clipped at the axis edge", async () => {
    mockApi(() => [
      row("MYCL", [cell("oe", { log2_fc: 8.238672784981318, q_value: 3.1e-10 }), cell("kd", { log2_fc: -3.0389494089903764, q_value: 2.4e-9 })]),
      row("PPP2R5B", [cell("oe", { log2_fc: 1.2612231780984975, q_value: 3.1e-4 }), cell("kd", { log2_fc: -0.8858491596285312, q_value: 0.0011 })]),
      row("PPP2R5E", [cell("oe", { log2_fc: -0.00670888680245163, q_value: 0.95 }), cell("kd", { log2_fc: -0.11047032000759316, q_value: 0.48 })]),
    ]);
    renderPage();
    await runWith("MYCL PPP2R5B PPP2R5E");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    const bars = Array.from(chartSvg().querySelectorAll("[data-testid=fc-bar]"));
    expect(bars).toHaveLength(6);
    const W = Number(chartSvg().getAttribute("width"));
    const bar = (gene: string, contrast: string) =>
      bars.find((b) => b.getAttribute("data-gene") === gene && b.getAttribute("data-contrast") === contrast)!;
    const len = (b: Element) => Number(b.getAttribute("width"));
    // Pixels per log2 unit, from a bar well inside any plausible cap.
    const ref = bar("PPP2R5B", "oe");
    const pxPerUnit = len(ref) / 1.2612231780984975;
    for (const b of bars) {
      const fc = Number(b.getAttribute("data-log2fc"));
      const x = Number(b.getAttribute("x"));
      // Bars under 1 px are drawn at 1 px so they stay visible.
      expect(len(b)).toBeCloseTo(Math.max(1, Math.abs(fc) * pxPerUnit), 3);
      expect(x).toBeGreaterThanOrEqual(0);
      expect(x + len(b)).toBeLessThanOrEqual(W);
    }
    // The largest bar is 8.24 / 1.26 = 6.5x the reference, not capped near 4/1.26.
    expect(len(bar("MYCL", "oe")) / len(ref)).toBeCloseTo(8.238672784981318 / 1.2612231780984975, 3);
    // Its value is printed, with no truncation arrow.
    expect(chartSvg().textContent).toContain("+8.24");
    expect(chartSvg().textContent).not.toContain("›");
  });

  it("titles the axis in sentence case in the body font, and keeps the researcher's headings as typed", async () => {
    mockApi((genes) => genes.map((g) => row(g, [cell("oe"), cell("kd")])));
    renderPage();
    await runWith("MYC family: MYC, MYCL");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    const title = within(chartSvg() as unknown as HTMLElement).getByTestId("fc-axis-title");
    expect(title.textContent).toMatch(/^log2 fold change \(treatment − reference\)/);
    expect(title.getAttribute("font-family")).toBe("var(--f-body)");
    expect(chartSvg().textContent).toContain("MYC family");
    expect(chartSvg().textContent).not.toMatch(/LOG2|MYC FAMILY/);
  });

  /**
   * Shipped without a key: colour, fade and the printed numbers were
   * explained only in an uppercase axis footer, and the TPM-vs-counts
   * direction flag (METHODS.md 11.3) not at all. The key is inside the
   * SVG so the exported PNG carries it.
   */
  it("has a key inside the figure for colour, significance, values and experiments", async () => {
    mockApi(() => [row("SET", [cell("oe", { log2_fc: 0.5, q_value: 0.01 }), cell("kd", { log2_fc: -0.2, q_value: 0.4 })])]);
    renderPage();
    await runWith("SET");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    const legend = chartSvg().querySelector("[data-testid=fc-legend]")!;
    expect(legend).not.toBeNull();
    const marks = Array.from(legend.querySelectorAll("[data-legend-mark]")).map((m) => m.getAttribute("data-legend-mark"));
    expect(marks).toEqual(["hot", "cool", "sig", "nonsig", "value", "tags"]);
    const text = legend.textContent!;
    expect(text).toMatch(/higher in treatment/);
    expect(text).toMatch(/q < 0\.05/);
    expect(text).toMatch(/q ≥ 0\.05/);
    expect(text).toMatch(/no interval drawn/);
    expect(text).toContain("C4-2B / NCI-H660");
    // The key's swatches are drawn as the bars are: solid for q < 0.05, faded otherwise.
    const opacity = (sel: string) => legend.querySelector(`[data-legend-mark=${sel}] rect`)!.getAttribute("fill-opacity");
    const barOpacity = (contrast: string) => chartSvg().querySelector(`[data-testid=fc-bar][data-contrast=${contrast}]`)!.getAttribute("fill-opacity");
    expect(opacity("sig")).toBe(barOpacity("oe"));
    expect(opacity("nonsig")).toBe(barOpacity("kd"));
    // Nothing untested, no direction flag: those entries are left out.
    expect(text).not.toMatch(/no test result|moves the other way/);
  });

  it("keys the direction flag and untested cells only when the figure shows them", async () => {
    mockApi(() => [
      row("HSPA12A", [cell("oe", { log2_fc: 1.04, display_direction_agrees: false }), cell("kd")]),
      row("MYCN", [cell("oe", { tested: false, log2_fc: null, q_value: null, display_direction_agrees: null }), cell("kd")]),
    ]);
    renderPage();
    await runWith("HSPA12A MYCN");
    await waitFor(() => expect(chartSvg()).toBeInTheDocument());
    const legend = chartSvg().querySelector("[data-testid=fc-legend]")!;
    expect(legend.textContent).toMatch(/mean TPM moves the other way from the count-based log2FC/);
    expect(legend.textContent).toMatch(/no test result in that experiment/);
    // The flagged bar's printed value carries the mark the key explains; the unflagged one does not.
    const value = (gene: string, contrast: string) =>
      chartSvg().querySelector(`[data-testid=fc-value][data-gene=${gene}][data-contrast=${contrast}]`)!.textContent;
    expect(value("HSPA12A", "oe")).toBe("+1.04 †");
    expect(value("HSPA12A", "kd")).toBe("+1.00");
  });

  /**
   * Shipped without a caption: the PNG named neither n per arm nor the
   * method's q scope in words, nor what solid/faded meant. The caption on
   * screen and the one printed into the PNG must be the same string
   * (METHODS.md 5.11), so this asserts on what reaches the exporter.
   */
  it("captions the figure with n, method, q and the marks, and exports exactly that caption", async () => {
    const exportSpy = vi.spyOn(exportFigure, "exportSvgAsPng").mockResolvedValue(undefined);
    mockApi(() => [
      row("MYCL", [cell("oe", { log2_fc: 8.24 }), cell("kd", { log2_fc: -3.04 })]),
      row("MYCN", [cell("oe", { tested: false, log2_fc: null, q_value: null, display_direction_agrees: null }), cell("kd")]),
    ]);
    renderPage();
    await runWith("MYCL MYCN");
    const caption = (await screen.findByText(/^log2 fold change \(treatment − reference\) for each gene/)).textContent!;
    expect(caption).toContain("C4-2B = C4-2B contrast (n = 2 vs 2)");
    expect(caption).toContain("NCI-H660 = NCI-H660 contrast (n = 2 vs 2)");
    expect(caption).toContain("fitted on its own samples only: limma-trend moderated t");
    expect(caption).toContain("q is Benjamini–Hochberg across all genes tested in that experiment (80 in C4-2B, 80 in NCI-H660)");
    expect(caption).toMatch(/solid bars have q < 0\.05, faded bars q ≥ 0\.05/);
    expect(caption).toMatch(/no interval is drawn/);
    expect(caption).toMatch(/“not tested”/);
    expect(caption).toContain("MYCL is the perturbed gene itself");
    // No direction flag in this figure, so no dagger sentence.
    expect(caption).not.toContain("†");

    await userEvent.click(screen.getAllByRole("button", { name: /png/i })[0]);
    await waitFor(() => expect(exportSpy).toHaveBeenCalled());
    const options = exportSpy.mock.calls.at(-1)![2] as { caption?: string };
    expect(options.caption).toBe(caption);
  });

  it("spans the axis to the data, with round ticks covering every value", async () => {
    const { foldChangeAxis } = await import("../ContrastsPage");
    const a = foldChangeAxis([8.24, -3.04, 1.26, -0.89], 906);
    expect(a.lo).toBeLessThanOrEqual(-3.04);
    expect(a.hi).toBeGreaterThanOrEqual(8.24);
    expect(a.ticks[0]).toBe(a.lo);
    expect(a.ticks.at(-1)).toBe(a.hi);
    expect(a.ticks).toContain(0);
    // Small panels keep at least -1..+1.
    const s = foldChangeAxis([0.2, -0.1], 600);
    expect([s.lo, s.hi]).toEqual([-1, 1]);
    // Narrow plots get coarser steps, never ticks closer than ~44 px.
    const n = foldChangeAxis([8.24, -3.04], 324);
    expect((n.step * 324) / (n.hi - n.lo)).toBeGreaterThanOrEqual(44);
  });
});
