/**
 * The PNG export rasterises a clone of the chart SVG outside the page,
 * where var(--token) no longer resolves. Paint set through an inline
 * style (style="fill: var(--hot)") overrode the literal fill attribute
 * written by inlineComputedPaint, so exported bars went black and a
 * white value label on a black bar disappeared (Perturbation Contrasts,
 * MYCL "+8.24 ›").
 */
import { describe, expect, it } from "vitest";

import { inlineComputedPaint, localDate } from "../exportFigure";

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

describe("export date stamp", () => {
  it("uses the reader's local calendar date, not UTC", () => {
    // 21:30 local on 15 Nov: toISOString() would say 16 Nov in any zone
    // behind UTC. Built from local fields so the test holds in every zone.
    expect(localDate(new Date(2026, 10, 15, 21, 30))).toBe("2026-11-15");
    expect(localDate(new Date(2026, 0, 5, 0, 5))).toBe("2026-01-05");
  });
});

describe("exportSvgAsPng caption", () => {
  // A figure pasted on a slide is only self-explanatory if the legend
  // travels with it: the caption shown under the chart must be printed
  // into the PNG, not lost at export (METHODS.md 5.11). jsdom has no
  // canvas, so the 2D context, image decode and download are stubbed and
  // the text actually drawn is recorded.
  it("draws the full caption text into the exported image", async () => {
    const { exportSvgAsPng } = await import("../exportFigure");
    const drawn: string[] = [];
    const fakeCtx = new Proxy(
      {
        fillText: (t: string) => drawn.push(t),
        measureText: (t: string) => ({ width: t.length * 6 }),
      } as Record<string, unknown>,
      { get: (target, key: string) => (key in target ? target[key] : () => undefined), set: () => true },
    );
    const restore: (() => void)[] = [];
    const stub = <T extends object>(obj: T, key: keyof T, value: unknown) => {
      const orig = obj[key];
      (obj as Record<keyof T, unknown>)[key] = value;
      restore.push(() => ((obj as Record<keyof T, unknown>)[key] = orig));
    };
    stub(HTMLCanvasElement.prototype, "getContext", () => fakeCtx);
    stub(HTMLCanvasElement.prototype, "toBlob", (cb: (b: Blob) => void) => cb(new Blob(["png"])));
    stub(HTMLAnchorElement.prototype, "click", () => undefined);
    stub(URL, "createObjectURL", () => "blob:stub");
    stub(URL, "revokeObjectURL", () => undefined);
    class FakeImage {
      onload: (() => void) | null = null;
      width = 0;
      height = 0;
      set src(_v: string) {
        setTimeout(() => this.onload?.(), 0);
      }
    }
    stub(globalThis, "Image" as keyof typeof globalThis, FakeImage);

    const svg = document.createElementNS(NS, "svg") as SVGSVGElement;
    svg.setAttribute("viewBox", "0 0 400 200");
    document.body.appendChild(svg);
    const caption =
      "MYCN expression (TPM) by etp_status. Boxes show the median and interquartile range, whiskers the full range. " +
      "Kruskal–Wallis test across groups; pairwise two-sided Mann–Whitney U with Benjamini–Hochberg FDR.";
    try {
      await exportSvgAsPng(svg, "x.png", { title: "MYCN — TARGET", statLines: ["p = 2.89e-8"], caption });
    } finally {
      restore.reverse().forEach((r) => r());
      svg.remove();
    }
    expect(drawn).toContain("p = 2.89e-8");
    // Wrapped across lines, but every word drawn, in order.
    const start = drawn.findIndex((t) => t.startsWith("MYCN expression"));
    expect(start).toBeGreaterThan(-1);
    const captionDrawn = drawn.slice(start).filter((t) => !t.startsWith("Expression Explorer ·"));
    expect(captionDrawn.length).toBeGreaterThan(1);
    expect(captionDrawn.join(" ")).toBe(caption);
  });
});
