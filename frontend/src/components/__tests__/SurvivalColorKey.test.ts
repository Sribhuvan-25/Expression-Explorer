/**
 * Exported survival PNGs carried no legend: nothing in the image said
 * which curve was HIGH and which LOW. The text key must match the plot's
 * colours (HIGH = --hot/red, LOW = --cool/blue) and give each arm's n.
 */
import { describe, expect, it } from "vitest";

import { survivalColorKey } from "../SurvivalCurve";

describe("survivalColorKey", () => {
  it("names each curve's colour and size", () => {
    const curves = { HIGH: { n: 233, points: [] }, LOW: { n: 233, points: [] } };
    expect(survivalColorKey(curves)).toBe("red = HIGH (n=233), blue = LOW (n=233)");
  });
});
