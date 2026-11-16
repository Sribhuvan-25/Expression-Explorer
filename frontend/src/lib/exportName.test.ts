import { describe, expect, it } from "vitest";

import { geneListTag } from "./exportName";

describe("geneListTag", () => {
  it("is stable for one list and differs between lists", () => {
    expect(geneListTag(["MYCN", "MYC"])).toBe(geneListTag(["MYCN", "MYC"]));
    expect(geneListTag(["MYCN", "MYC"])).not.toBe(geneListTag(["MYCN", "PPP2R5B"]));
    expect(geneListTag(["A", "B", "C"])).toMatch(/^3genes-[0-9a-z]+$/);
  });
});
