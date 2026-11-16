/**
 * Guards the degenerate-Cox-model rendering gap on the Survival page.
 *
 * A tied-scores policy (e.g. "exclude") on a quartile cutoff can empty
 * one arm entirely when many samples share the exact cutoff score --
 * common for a single gene, whose score is exactly 0 below the rank
 * threshold. The backend still returns a fitted Cox coefficient in that
 * case, but it's comparing a group against nothing: hazard ratio and CI
 * span many orders of magnitude (observed live: HR=18341.107,
 * CI 0.000-688692843566556.625, p=0.4294 on target_all_p2/MEF2C with
 * quartile+exclude). Rendered as an ordinary-looking statistics table,
 * that reads as a real result rather than a degenerate one. Found during
 * QA regression testing, 2026-09-23.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { SurvivalPage } from "../SurvivalPage";
import { api } from "../../lib/api";
import type { DatasetSummary, SurvivalResult } from "../../lib/api";

const TARGET_ALL_P2: DatasetSummary = {
  dataset_id: "target_all_p2",
  display_name: "TARGET ALL-P2 (paediatric T-ALL)",
  group_columns: ["vital_status"],
  supports_survival: true,
  n_samples: 469,
  assay_type: "rna_seq",
};

const DEGENERATE_RESULT: SurvivalResult = {
  genes: ["MEF2C"],
  n: 117,
  n_dataset_total: 469,
  n_excluded: 352,
  exclusion_reason: "no follow-up time recorded, or score fell in the excluded middle band",
  curves: { HIGH: { n: 117, points: [{ days: 0, survival_probability: 1 }] } },
  logrank_p_value: null as unknown as number,
  tie_note: {
    n_tied_at_cutoff: 218,
    tie_policy: "exclude",
    message: "218 samples share the exact score at a cutoff boundary. Applied tie policy: 'exclude'.",
  },
  cox: {
    n: 117,
    coefficients: {
      signature_score: {
        coef: 9.8169,
        "exp(coef)": 18341.107,
        "exp(coef) lower 95%": 4.8846e-7,
        "exp(coef) upper 95%": 688692843566556.6,
        p: 0.4294,
      },
    },
    log_likelihood_p_value: 0.4274,
  },
} as unknown as SurvivalResult;

const NORMAL_RESULT: SurvivalResult = {
  genes: ["MEF2C"],
  n: 466,
  n_dataset_total: 469,
  n_excluded: 3,
  exclusion_reason: "no follow-up time recorded",
  curves: {
    HIGH: { n: 233, points: [{ days: 0, survival_probability: 1 }] },
    LOW: { n: 233, points: [{ days: 0, survival_probability: 1 }] },
  },
  logrank_p_value: 1.93e-7,
  tie_note: null,
  cox: {
    n: 466,
    coefficients: {
      signature_score: {
        coef: 2.777,
        "exp(coef)": 16.068,
        "exp(coef) lower 95%": 4.244,
        "exp(coef) upper 95%": 60.838,
        p: 4.35e-5,
      },
    },
    log_likelihood_p_value: 2.77e-6,
  },
} as unknown as SurvivalResult;

function renderPage() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false, gcTime: 0 } } });
  return render(
    <QueryClientProvider client={client}>
      <SurvivalPage />
    </QueryClientProvider>,
  );
}

describe("SurvivalPage degenerate Cox model", () => {
  it("hides the Cox table and explains why when the tie policy leaves only one arm", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [TARGET_ALL_P2] });
    vi.spyOn(api, "survival").mockResolvedValue(DEGENERATE_RESULT);

    renderPage();

    fireEvent.click(await screen.findByRole("button", { name: /run survival analysis/i }));

    await waitFor(() => {
      expect(screen.getByText(/cox model not shown/i)).toBeInTheDocument();
    });

    // The absurd raw CI and coefficient table must never reach the page.
    expect(screen.queryByText(/688692843566556/)).not.toBeInTheDocument();
    expect(screen.queryByText(/18341\.107/)).not.toBeInTheDocument();
    expect(screen.queryByRole("columnheader", { name: /hazard ratio/i })).not.toBeInTheDocument();
    expect(screen.queryByText(/^cox model p$/i)).not.toBeInTheDocument();
  });

  it("shows the Cox table normally when both arms are populated", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [TARGET_ALL_P2] });
    vi.spyOn(api, "survival").mockResolvedValue(NORMAL_RESULT);

    renderPage();

    fireEvent.click(await screen.findByRole("button", { name: /run survival analysis/i }));

    await waitFor(() => {
      expect(screen.getByRole("columnheader", { name: /hazard ratio/i })).toBeInTheDocument();
    });
    expect(screen.getByText((_, node) => node?.textContent === "4.244 – 60.838")).toBeInTheDocument();
    expect(screen.queryByText(/cox model not shown/i)).not.toBeInTheDocument();
    // The HR is per 1.0 of a 0-1 score; without the unit it reads as a
    // huge per-unit effect.
    expect(screen.getByTestId("cox-unit-note")).toHaveTextContent("compares a score of 1 with a score of 0");
  });
});

describe("SurvivalPage sample accounting and captions follow the cutoff", () => {
  it("lists each exclusion reason on its own and names the quartile split", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [TARGET_ALL_P2] });
    vi.spyOn(api, "survival").mockResolvedValue({
      ...NORMAL_RESULT,
      n: 234,
      n_excluded: 235,
      exclusion_reason: "3 with no follow-up time recorded; 232 with a score between the low (25%) and high (25%) cutoffs",
      exclusion_breakdown: [
        { reason: "with no follow-up time recorded", n: 3 },
        { reason: "with a score between the low (25%) and high (25%) cutoffs", n: 232 },
      ],
    } as unknown as SurvivalResult);

    renderPage();
    fireEvent.click(await screen.findByRole("button", { name: "quartile" }));
    fireEvent.click(screen.getByRole("button", { name: /run survival analysis/i }));

    const list = await screen.findByTestId("survival-exclusions");
    expect(list.children).toHaveLength(2);
    expect(list).toHaveTextContent("3 with no follow-up time recorded");
    expect(list).toHaveTextContent("232 with a score between the low (25%) and high (25%) cutoffs");
    expect(screen.getByTestId("cox-unit-note")).toHaveTextContent("quartile");
    expect(screen.getByTestId("cox-unit-note")).not.toHaveTextContent("median-split");
  });
});

describe("SurvivalPage tie-break sensitivity", () => {
  it("shows how far the log-rank p moves with the choice of tied samples", async () => {
    vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [TARGET_ALL_P2] });
    vi.spyOn(api, "survival").mockResolvedValue({
      ...NORMAL_RESULT,
      logrank_p_value: 0.0041,
      tie_note: {
        n_tied_at_cutoff: 246,
        tie_policy: "trim",
        message: "246 samples share the exact score at a cutoff boundary.",
        sensitivity: { n_draws: 200, p_median: 0.044, p_low: 0.0034, p_high: 0.32 },
      },
    } as unknown as SurvivalResult);
    renderPage();
    fireEvent.click(await screen.findByRole("button", { name: "quartile" }));
    fireEvent.click(screen.getByRole("button", { name: /run survival analysis/i }));
    const note = await screen.findByTestId("tie-sensitivity");
    expect(note).toHaveTextContent("across 200 random choices, the log-rank p ranged 0.003 to 0.320");
    expect(note).toHaveTextContent("median 0.044");
    expect(note).toHaveTextContent("against 0.004 shown");
  });
});
