import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { BoxPlot, niceTicks, valueDomain } from "../BoxPlot";

// GSE66586 PPP2R5B by MYCN status, as served (log2 intensity).
const PPP2R5B = [
  ["non-malignant", 9.829],
  ["non-malignant", 10.037],
  ["MYCN-amplified", 8.307],
  ["MYCN-amplified", 7.534],
  ["MYCN-amplified", 8.032],
  ["MYCN-amplified", 8.309],
  ["MYCN-amplified", 8.195],
  ["MYCN non-amplified", 8.95],
  ["MYCN non-amplified", 9.499],
  ["MYCN non-amplified", 8.72],
].map(([group, value], i) => ({ sample_id: `GSM${i}`, group: group as string, value: value as number }));

function axisTicks(container: HTMLElement): number[] {
  // Tick labels are the mono-font texts left of the plot.
  return [...container.querySelectorAll("svg text[text-anchor='end'][dominant-baseline='middle']")]
    .map((t) => Number(t.textContent))
    .filter((v) => Number.isFinite(v));
}

describe("box plot value axis", () => {
  it("does not anchor log-scale units at 0, so group differences stay visible", () => {
    const { container } = render(<BoxPlot points={PPP2R5B} valueLabel="expression" unit="log2_intensity" />);
    const ticks = axisTicks(container);
    expect(Math.min(...ticks)).toBeGreaterThan(6);
    expect(Math.min(...ticks)).toBeLessThanOrEqual(7.534);
    expect(Math.max(...ticks)).toBeGreaterThanOrEqual(10.037);
    // The data fills most of the axis instead of its top quarter.
    expect((10.037 - 7.534) / (Math.max(...ticks) - Math.min(...ticks))).toBeGreaterThan(0.6);
  });

  it("names the unit on the axis and offers no log-of-log toggle", () => {
    render(<BoxPlot points={PPP2R5B} valueLabel="expression" unit="log2_intensity" />);
    // Sentence case in the body font: the title was uppercase tracked mono
    // ("EXPRESSION (LOG2 INTENSITY)"), which read as a code identifier.
    const title = screen.getByText("Expression (log2 intensity)");
    expect(title).toHaveAttribute("font-family", "var(--f-body)");
    expect(title.getAttribute("style") ?? "").not.toContain("letter-spacing");
    expect(screen.queryByText("Log scale")).toBeNull();
  });

  it("keeps 0 on the axis for linear units such as TPM", () => {
    const tpm = PPP2R5B.map((p) => ({ ...p, value: 2 ** p.value }));
    const { container } = render(<BoxPlot points={tpm} valueLabel="expression" unit="tpm" />);
    expect(Math.min(...axisTicks(container))).toBe(0);
    expect(screen.getByText("Expression (TPM)")).toBeInTheDocument();
    expect(screen.getByText("Log scale")).toBeInTheDocument();
  });

  it("puts ticks on round values", () => {
    expect(niceTicks(0, 10.4)).toEqual([0, 2.5, 5, 7.5, 10, 12.5]);
    expect(niceTicks(7.33, 10.24)).toEqual([7, 8, 9, 10, 11]);
    expect(niceTicks(-0.3, 0.7)).toEqual([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]);
    expect(niceTicks(5, 5).length).toBeGreaterThan(1);
  });

  it("pads a log-unit domain around the data, and a linear one from 0", () => {
    const [lo, hi] = valueDomain([7.534, 10.037], true);
    expect(lo).toBeLessThan(7.534);
    expect(lo).toBeGreaterThan(7);
    expect(hi).toBeGreaterThan(10.037);
    expect(valueDomain([3, 40], false)[0]).toBe(0);
  });
});

describe("box plot group colours", () => {
  it("gives three groups three distinguishable colours (no two blues)", () => {
    const { container } = render(<BoxPlot points={PPP2R5B} valueLabel="expression" unit="log2_intensity" />);
    const fills = [...container.querySelectorAll("svg rect[fill-opacity]")].map((r) => r.getAttribute("fill"));
    expect(fills).toHaveLength(3);
    expect(new Set(fills).size).toBe(3);
    // --accent and --cool render as near-identical blues.
    expect(fills.includes("var(--accent)") && fills.includes("var(--cool)")).toBe(false);
  });
});

describe("box plot tick labels", () => {
  it("prints ticks of one step at one precision", () => {
    const { container } = render(<BoxPlot points={PPP2R5B} valueLabel="expression" unit="log2_intensity" />);
    const labels = [...container.querySelectorAll("svg text[text-anchor='end'][dominant-baseline='middle']")].map(
      (t) => t.textContent,
    );
    expect(labels).toEqual(["7", "8", "9", "10", "11"]);
  });
});

describe("box plot axis for non-negative log units", () => {
  it("does not run below 0 when no value can", () => {
    const [lo] = valueDomain([0, 0, 0.2, 6.4], true);
    expect(lo).toBe(0);
    const pts = [0, 0, 0.1, 0.3, 2.2, 6.4].map((v, i) => ({ sample_id: `S${i}`, group: "Lymphoid", value: v }));
    const { container } = render(<BoxPlot points={pts} valueLabel="expression" unit="log2_tpm" />);
    expect(Math.min(...axisTicks(container))).toBe(0);
  });
});

describe("box plot log-scale ticks", () => {
  it("lands on round values", async () => {
    const { logTicks } = await import("../BoxPlot");
    expect(logTicks(148)).toEqual([0, 1, 3, 10, 30, 100]);
    expect(logTicks(6.4)).toEqual([0, 0.5, 1, 2, 5]);
    expect(logTicks(5000)).toEqual([0, 10, 100, 1000]);
    for (const max of [0.8, 3, 42, 900, 12000]) {
      const t = logTicks(max);
      expect(t.length).toBeGreaterThanOrEqual(2);
      expect(t.length).toBeLessThanOrEqual(7);
      expect(Math.max(...t)).toBeLessThanOrEqual(max);
    }
  });
});

describe("box plot with many groups on a narrow screen", () => {
  it("puts each group name above its box instead of in a side column", async () => {
    const { vi } = await import("vitest");
    const spy = vi.spyOn(HTMLElement.prototype, "clientWidth", "get").mockReturnValue(214);
    const pts = Array.from({ length: 12 }, (_, g) =>
      [1, 2, 3].map((k) => ({ sample_id: `S${g}-${k}`, group: `Subtype number ${g}`, value: g + k })),
    ).flat();
    const { container } = render(<BoxPlot points={pts} valueLabel="expression" unit="log2_tpm" />);
    const nameTexts = [...container.querySelectorAll("svg text")].filter((t) => /Subtype number/.test(t.textContent ?? ""));
    expect(nameTexts.length).toBe(12);
    expect(nameTexts.every((t) => t.getAttribute("text-anchor") === "start")).toBe(true);
    spy.mockRestore();
    const wide = render(<BoxPlot points={pts} valueLabel="expression" unit="log2_tpm" />);
    const wideNames = [...wide.container.querySelectorAll("svg text")].filter((t) => /Subtype number/.test(t.textContent ?? ""));
    expect(wideNames.every((t) => t.getAttribute("text-anchor") === "end")).toBe(true);
  });
});
