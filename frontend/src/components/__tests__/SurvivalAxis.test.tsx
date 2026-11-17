import { render } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { SurvivalPlot } from "../SurvivalCurve";

describe("survival day axis", () => {
  it("puts day ticks on round values (TARGET follow-up ends at 4383 days)", () => {
    const curve = (n: number) => ({
      n,
      points: [
        { days: 0, survival_probability: 1 },
        { days: 4383, survival_probability: 0.6 },
      ],
    });
    const { container } = render(<SurvivalPlot curves={{ HIGH: curve(233), LOW: curve(233) }} />);
    const labels = [...container.querySelectorAll("svg text[data-x-tick]")]
      .map((t) => t.textContent ?? "")
      .filter((t) => /^\d+$/.test(t));
    expect(labels).toEqual(["0", "1000", "2000", "3000", "4000"]);
  });
});
