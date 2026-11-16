/**
 * Segmented toggles showed the selected option only as a background tint.
 * Without aria-pressed a screen reader announced "median, button",
 * "quartile, button" with no way to tell which cutoff was in force.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { BoxPlot } from "../BoxPlot";
import { PresetButton } from "../ui";
import { SurvivalPage } from "../../pages/SurvivalPage";
import { api } from "../../lib/api";

describe("toggle buttons expose which option is selected", () => {
  it("survival cutoff and tie-policy toggles", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({
      datasets: [
        {
          dataset_id: "target_all_p2",
          display_name: "TARGET ALL-P2",
          group_columns: ["vital_status"],
          supports_survival: true,
          n_samples: 469,
          assay_type: "rna_seq",
        },
      ],
    } as never);
    render(
      <QueryClientProvider client={new QueryClient()}>
        <SurvivalPage />
      </QueryClientProvider>,
    );
    const median = await screen.findByRole("button", { name: "median" });
    const quartile = screen.getByRole("button", { name: "quartile" });
    expect(median).toHaveAttribute("aria-pressed", "true");
    expect(quartile).toHaveAttribute("aria-pressed", "false");
    fireEvent.click(quartile);
    expect(quartile).toHaveAttribute("aria-pressed", "true");
    expect(median).toHaveAttribute("aria-pressed", "false");
    expect(screen.getByRole("button", { name: "Trim to size" })).toHaveAttribute("aria-pressed", "true");
  });

  it("box plot linear / log scale", () => {
    const points = [1, 2, 3, 50, 400].flatMap((v, i) => [
      { sample_id: `a${i}`, group: "A", value: v },
      { sample_id: `b${i}`, group: "B", value: v * 2 },
    ]);
    render(<BoxPlot points={points} valueLabel="expression" unit="tpm" />);
    const linear = screen.getByRole("button", { name: "Linear" });
    const log = screen.getByRole("button", { name: "Log scale" });
    fireEvent.click(linear);
    expect(linear).toHaveAttribute("aria-pressed", "true");
    expect(log).toHaveAttribute("aria-pressed", "false");
    fireEvent.click(log);
    expect(log).toHaveAttribute("aria-pressed", "true");
  });

  it("preset buttons", () => {
    render(<PresetButton label="ETP-TF5" active onClick={() => {}} />);
    expect(screen.getByRole("button", { name: "ETP-TF5" })).toHaveAttribute("aria-pressed", "true");
  });
});
