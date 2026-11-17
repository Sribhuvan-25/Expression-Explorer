/**
 * The Kaplan–Meier figure must name its endpoint, mark censored patients
 * and carry a numbers-at-risk table -- all inside the SVG, so the exported
 * PNG has them too. A clinician reviewer could not show or publish the
 * figure without them. Axis titles are sentence case in the body font
 * (they were uppercase tracked mono, "SURVIVAL PROBABILITY").
 */
import { render } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { SurvivalPlot, atRiskAt } from "../SurvivalCurve";
import type { SurvivalCurve } from "../../lib/api";

// Hand-worked arm (see backend/tests/test_survival_figure.py): deaths at
// 1, 3, 5, 10; censored at 3 and 8. At risk = follow-up >= t.
const AT_RISK = [
  { days: 0, n_at_risk: 6 },
  { days: 1, n_at_risk: 6 },
  { days: 3, n_at_risk: 5 },
  { days: 5, n_at_risk: 3 },
  { days: 8, n_at_risk: 2 },
  { days: 10, n_at_risk: 1 },
];

describe("atRiskAt", () => {
  it("counts patients with follow-up at least t, between and at observed times", () => {
    expect(atRiskAt(AT_RISK, 0)).toBe(6);
    expect(atRiskAt(AT_RISK, 2)).toBe(5); // between 1 and 3: the 5 with follow-up >= 3
    expect(atRiskAt(AT_RISK, 3)).toBe(5); // a death and a censoring at 3 both count at 3
    expect(atRiskAt(AT_RISK, 9)).toBe(1);
    expect(atRiskAt(AT_RISK, 10)).toBe(1);
    expect(atRiskAt(AT_RISK, 11)).toBe(0); // past the last follow-up
  });
});

function curve(n: number, extra: Partial<SurvivalCurve> = {}): SurvivalCurve {
  return {
    n,
    points: [
      { days: 0, survival_probability: 1 },
      { days: 4383, survival_probability: 0.6 },
    ],
    ...extra,
  };
}

const HIGH = curve(233, {
  censor_marks: [
    { days: 500, survival_probability: 0.9, n_censored: 1 },
    { days: 2500, survival_probability: 0.7, n_censored: 2 },
  ],
  // Steps placed just before, at and after ticks, so a tick that falls
  // between observed times must take the next step (999 -> 174 is not it).
  at_risk: [
    { days: 0, n_at_risk: 233 },
    { days: 999, n_at_risk: 174 },
    { days: 1001, n_at_risk: 173 },
    { days: 2000, n_at_risk: 127 },
    { days: 3001, n_at_risk: 57 },
    { days: 3990, n_at_risk: 9 },
    { days: 4383, n_at_risk: 1 },
  ],
});
const LOW = curve(233, {
  censor_marks: [{ days: 1200, survival_probability: 0.95, n_censored: 1 }],
  at_risk: [
    { days: 0, n_at_risk: 233 },
    { days: 1000, n_at_risk: 215 },
    { days: 2000, n_at_risk: 168 },
    { days: 3000, n_at_risk: 18 },
    { days: 3800, n_at_risk: 2 },
  ],
});

describe("Kaplan–Meier figure", () => {
  it("names the endpoint in sentence-case axis titles in the body font", () => {
    const { container } = render(
      <SurvivalPlot curves={{ HIGH, LOW }} endpointName="Overall survival" timeUnit="days" />,
    );
    const y = container.querySelector("[data-axis-title='y']")!;
    const x = container.querySelector("[data-axis-title='x']")!;
    expect(y.textContent).toBe("Overall survival probability");
    expect(x.textContent).toBe("Days");
    for (const t of [x, y]) {
      expect(t.getAttribute("font-family")).toBe("var(--f-body)");
      expect((t as SVGElement).style.letterSpacing).toBe("");
    }
  });

  it("draws one censoring tick per censoring time, in its curve's colour", () => {
    const { container } = render(<SurvivalPlot curves={{ HIGH, LOW }} endpointName="Overall survival" />);
    const high = container.querySelectorAll("svg [data-censor-mark='HIGH']");
    const low = container.querySelectorAll("svg [data-censor-mark='LOW']");
    expect(high).toHaveLength(2);
    expect(low).toHaveLength(1);
    expect(high[0].getAttribute("stroke")).toBe("var(--hot)");
    expect(low[0].getAttribute("stroke")).toBe("var(--cool)");
  });

  it("puts a numbers-at-risk row per curve inside the SVG, at the day ticks", () => {
    const { container } = render(<SurvivalPlot curves={{ HIGH, LOW }} endpointName="Overall survival" />);
    const table = container.querySelector("svg [data-testid='at-risk-table']")!;
    expect(table).not.toBeNull();
    const row = (label: string) =>
      [...table.querySelectorAll(`[data-at-risk-row='${label}'] [data-at-risk-day]`)].map((t) => [
        t.getAttribute("data-at-risk-day"),
        t.textContent,
      ]);
    expect(row("HIGH")).toEqual([
      ["0", "233"],
      ["1000", "173"],
      ["2000", "127"],
      ["3000", "57"],
      ["4000", "1"],
    ]);
    expect(row("LOW")).toEqual([
      ["0", "233"],
      ["1000", "215"],
      ["2000", "168"],
      ["3000", "18"],
      ["4000", "0"],
    ]);
    const labels = table.querySelectorAll("[data-at-risk-row] > text:first-child");
    expect([...labels].map((l) => [l.textContent, l.getAttribute("fill")])).toEqual([
      ["HIGH", "var(--hot)"],
      ["LOW", "var(--cool)"],
    ]);
  });

  it("still renders the curves from a response without the figure fields", () => {
    const { container } = render(<SurvivalPlot curves={{ HIGH: curve(233), LOW: curve(233) }} />);
    expect(container.querySelector("[data-testid='at-risk-table']")).toBeNull();
    expect(container.querySelectorAll("[data-censor-mark]")).toHaveLength(0);
    expect(container.querySelector("[data-axis-title='y']")!.textContent).toBe("Survival probability");
    expect(container.querySelectorAll("svg path")).toHaveLength(2);
  });
});
