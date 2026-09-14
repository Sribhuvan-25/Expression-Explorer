import "@testing-library/jest-dom/vitest";
import { cleanup } from "@testing-library/react";
import { afterEach, vi } from "vitest";

// jsdom implements neither of these, and the chart components use both to
// size themselves against their container. Without stubs, rendering any
// chart throws before a test can assert anything.
if (!("ResizeObserver" in globalThis)) {
  globalThis.ResizeObserver = class {
    observe() {}
    unobserve() {}
    disconnect() {}
  } as unknown as typeof ResizeObserver;
}

// jsdom gives every element a zero layout box, so charts that scale to
// their container would render at width 0. A fixed non-zero width keeps
// them measurable without pretending to do real layout.
Object.defineProperty(HTMLElement.prototype, "clientWidth", {
  configurable: true,
  get() {
    return 800;
  },
});

afterEach(() => {
  cleanup();
  vi.restoreAllMocks();
});
