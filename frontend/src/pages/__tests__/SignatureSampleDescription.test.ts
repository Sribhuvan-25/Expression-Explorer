/**
 * The Signature Score table listed bare GSM ids for the cell-line datasets
 * ("GSM1625938 0.1996"), so a reader could not tell which line scored
 * highest without cross-referencing GEO. Rows now carry a description
 * built from the sample's own metadata.
 */
import { describe, expect, it } from "vitest";

import { sampleDescription } from "../SignaturePage";

describe("sampleDescription", () => {
  it("leads with the cell line and adds the grouping values", () => {
    expect(
      sampleDescription({
        sample_id: "GSM1625938", value: 0.2, mycn_status: "MYCN-amplified", cell_line: "BE(2)c",
        cell_type: "neuroblastoma", cell_line_identity: "SK-N-BE(2)-C",
      }),
    ).toBe("BE(2)c · MYCN-amplified · neuroblastoma · SK-N-BE(2)-C");
  });

  it("drops bookkeeping fields, follow-up numbers and duplicates", () => {
    expect(sampleDescription({ sample_id: "X", value: 1, sample_title: "H660 shControl replicate 1" })).toBe("H660 shControl replicate 1");
    expect(
      sampleDescription({
        sample_id: "GSM9589990", value: 0.4, sample_title: "C4-2B MYCL overexpression replicate 1",
        library_name: "C4-2B_MYCL_1", cell_line: "C4-2B", arm: "C4-2B MYCL overexpression", treatment: "MYCL overexpression",
      }),
    ).toBe("C4-2B · C4-2B MYCL overexpression · MYCL overexpression"); // title only used when no cell line
    expect(sampleDescription({ sample_id: "TARGET-1", value: 1, etp_status: "ETP", days_to_death: 300, vital_status: null })).toBe("ETP");
    expect(sampleDescription(undefined)).toBe("");
  });
});
