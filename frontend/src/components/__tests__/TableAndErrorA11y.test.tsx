import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";

import { DataTable, ErrorState } from "../ui";

describe("table and error accessibility", () => {
  it("table scroll region is focusable and empty headers are named", () => {
    render(<DataTable columns={[{ key: "gene", label: "Gene" }, { key: "score", label: "" }]} rows={[{ gene: "MYCN", score: 1 }]} />);
    expect(screen.getByRole("region", { name: "Results table" })).toHaveAttribute("tabindex", "0");
    expect(screen.getAllByRole("columnheader")[1]).toHaveTextContent("score");
  });

  it("errors are announced", () => {
    render(<ErrorState message="Gene 'X' not found in dataset." />);
    expect(screen.getByRole("alert")).toHaveTextContent("Gene 'X' not found");
  });
});
