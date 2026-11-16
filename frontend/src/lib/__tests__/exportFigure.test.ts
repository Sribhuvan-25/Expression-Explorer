/**
 * The PNG export rasterises a clone of the chart SVG outside the page,
 * where var(--token) no longer resolves. Paint set through an inline
 * style (style="fill: var(--hot)") overrode the literal fill attribute
 * written by inlineComputedPaint, so exported bars went black and a
 * white value label on a black bar disappeared (Perturbation Contrasts,
 * MYCL "+8.24 ›").
 */
import { describe, expect, it } from "vitest";

import { inlineComputedPaint } from "../exportFigure";

const NS = "http://www.w3.org/2000/svg";

describe("inlineComputedPaint", () => {
  it("overrides an unresolvable inline-style paint on the clone with the live resolved value", () => {
    const live = document.createElementNS(NS, "svg") as SVGSVGElement;
    const rect = document.createElementNS(NS, "rect");
    rect.setAttribute("style", "fill: rgb(200, 30, 30); stroke: rgb(10, 20, 30)");
    live.appendChild(rect);
    document.body.appendChild(live);

    const clone = live.cloneNode(true) as SVGSVGElement;
    // What the serialized clone effectively carries: a var() it cannot resolve.
    (clone.firstChild as SVGElement).setAttribute("style", "fill: var(--hot); stroke: var(--rule)");

    inlineComputedPaint(live, clone);
    const out = clone.firstChild as SVGElement;
    expect(out.style.fill).toBe("rgb(200, 30, 30)");
    expect(out.style.stroke).toBe("rgb(10, 20, 30)");
    expect(out.getAttribute("fill")).toBe("rgb(200, 30, 30)");
    live.remove();
  });
});

describe("wrapText", () => {
  const measure = (s: string) => s.length * 10; // 10px per character

  it("wraps a long title instead of letting it run off the image", async () => {
    const { wrapText } = await import("../exportFigure");
    const title = "PCA — MEF2C, LYL1, HHEX, LMO2 +1 more — GSE66586 (neuroblastoma cell lines, MYCN status)";
    const lines = wrapText(title, 300, measure);
    expect(lines.length).toBeGreaterThan(1);
    expect(lines.every((l) => measure(l) <= 300 || !l.includes(" "))).toBe(true);
    expect(lines.join(" ")).toBe(title); // nothing dropped
  });

  it("keeps a short line intact and never drops an over-long word", async () => {
    const { wrapText } = await import("../exportFigure");
    expect(wrapText("p = 0.0357", 300, measure)).toEqual(["p = 0.0357"]);
    expect(wrapText("a ENSG00000310517ENSG00000310517 b", 100, measure)).toEqual(["a", "ENSG00000310517ENSG00000310517", "b"]);
  });
});
