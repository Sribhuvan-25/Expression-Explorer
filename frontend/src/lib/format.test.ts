import { describe, expect, it } from "vitest";

import { formatP } from "./format";

describe("formatP", () => {
  it("uses fixed notation from 0.001 up", () => {
    expect(formatP(0.751)).toBe("0.7510");
    expect(formatP(0.0009)).toBe("9.00e-4");
    expect(formatP(null)).toBe("—");
  });
});
