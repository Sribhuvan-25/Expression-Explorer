import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { ScoreDistribution } from "../SignaturePage";

describe("signature score histogram", () => {
  it("names the score and gives a count scale", () => {
    const scores = [0, 0, 0, 0, 0.1, 0.2, 0.25, 0.5];
    render(<ScoreDistribution scores={scores} shown={3} scoreLabel="AUCell score (0–1)" />);
    expect(screen.getByText("AUCell score (0–1)")).toBeInTheDocument();
    expect(screen.getByTestId("hist-peak")).toHaveTextContent("tallest bar 4");
    expect(screen.getByRole("img")).toHaveAccessibleName(/Histogram of 8 aucell score \(0–1\)s from 0.000 to 0.500/);
  });
});
