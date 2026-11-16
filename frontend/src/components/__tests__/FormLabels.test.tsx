/**
 * Every form control needs an accessible name. The visible labels were
 * plain <label> elements not tied to their control, so a screen reader
 * announced "combo box" with no idea which (axe: select-name).
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { ComparePage } from "../../pages/ComparePage";
import { RankPage } from "../../pages/RankPage";
import { SurvivalPage } from "../../pages/SurvivalPage";
import { SignaturePage } from "../../pages/SignaturePage";
import { api } from "../../lib/api";

const DATASETS = {
  datasets: [
    {
      dataset_id: "target_all_p2",
      display_name: "TARGET",
      group_columns: ["etp_status"],
      supports_survival: true,
      n_samples: 469,
      assay_type: "rna_seq",
      expression_unit: "tpm",
    },
  ],
} as never;

function wrap(node: React.ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

describe("form controls have accessible names", () => {
  for (const [name, Page] of [
    ["Expression Compare", ComparePage],
    ["Sample Ranking", RankPage],
    ["Survival", SurvivalPage],
    ["Signature Score", SignaturePage],
  ] as const) {
    it(name, async () => {
      vi.spyOn(api, "listDatasets").mockResolvedValue(DATASETS);
      vi.spyOn(api, "groupValues").mockResolvedValue({ group_column: "etp_status", values: [], unavailable: false });
      wrap(<Page />);
      await screen.findAllByRole("button");
      const controls = [
        ...screen.queryAllByRole("combobox"),
        ...screen.queryAllByRole("textbox"),
        ...screen.queryAllByRole("spinbutton"),
      ];
      expect(controls.length).toBeGreaterThan(0);
      for (const c of controls) expect(c).toHaveAccessibleName();
    });
  }
});
