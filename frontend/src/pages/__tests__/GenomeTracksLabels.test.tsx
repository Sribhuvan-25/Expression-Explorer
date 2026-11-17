/**
 * igv.js leaves its chromosome menu and zoom slider unlabelled; the page
 * names them. Labelling used to start only once igv.createBrowser()
 * resolved -- which waits for every track to load -- while igv.js builds
 * the widgets earlier, so on a slow load (WebKit with several panes open)
 * axe found them unlabelled. It also stopped after its first success, so
 * a widget igv.js re-created later stayed unlabelled.
 */
import { act, render } from "@testing-library/react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";

vi.mock("igv", () => ({
  default: {
    // Never resolves: the tracks are "still loading".
    createBrowser: vi.fn((container: HTMLElement) => {
      const host = document.createElement("div");
      host.className = "igv-container";
      container.appendChild(host);
      const root = host.attachShadow({ mode: "open" });
      setTimeout(() => {
        root.innerHTML = '<select name="chromosome-select-widget"></select><input type="range" />';
      }, 100);
      return new Promise(() => {});
    }),
    removeBrowser: vi.fn(),
  },
}));

import { GenomeTracksPage, watchIgvWidgets } from "../GenomeTracksPage";

beforeEach(() => {
  vi.useFakeTimers();
  vi.spyOn(HTMLElement.prototype, "getBoundingClientRect").mockReturnValue({ width: 1200, height: 600 } as DOMRect);
});
afterEach(() => vi.useRealTimers());

const shadow = () => [...document.querySelectorAll(".igv-container")].find((el) => el.shadowRoot)!.shadowRoot!;

describe("igv.js widget labels", () => {
  it("labels the widgets while the tracks are still loading", async () => {
    render(<GenomeTracksPage />);
    await act(async () => {
      await vi.advanceTimersByTimeAsync(1500);
    });
    expect(shadow().querySelector("select")).toHaveAttribute("aria-label", "Chromosome");
    expect(shadow().querySelector('input[type="range"]')).toHaveAttribute("aria-label", "Zoom level");
  });

  it("labels a widget igv.js re-creates after the first labelling", async () => {
    vi.useRealTimers();
    const host = document.createElement("div");
    host.className = "igv-container";
    document.body.appendChild(host);
    const root = host.attachShadow({ mode: "open" });
    root.innerHTML = '<select name="chromosome-select-widget"></select><input type="range" />';
    const stop = watchIgvWidgets();
    try {
      expect(root.querySelector("select")).toHaveAttribute("aria-label", "Chromosome");
      await new Promise((r) => setTimeout(r, 600)); // polling has stopped by now
      root.querySelector("select")!.replaceWith(Object.assign(document.createElement("select"), { name: "chromosome-select-widget" }));
      await new Promise((r) => setTimeout(r, 20)); // observer callback
      expect(root.querySelector("select")).toHaveAttribute("aria-label", "Chromosome");
    } finally {
      stop();
      host.remove();
    }
  });
});
