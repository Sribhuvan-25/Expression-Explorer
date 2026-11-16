/**
 * Typing "CIP2A SET ANP32A ARPP19 PABIR1" and pressing Enter produced ONE
 * tag, "CIP2A SET ANP32A ARPP19 PABIR1", which no dataset can resolve, and
 * a pasted list did the same (found scoring the PP2A inhibitors on
 * GSE66586). Text left in the box when clicking Score was dropped.
 */
import { fireEvent, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { useState } from "react";
import { describe, expect, it } from "vitest";

import { GeneTagInput } from "../ui";

function Harness({ initial = [] as string[] }) {
  const [genes, setGenes] = useState<string[]>(initial);
  return (
    <>
      <GeneTagInput genes={genes} onChange={setGenes} />
      <output data-testid="genes">{genes.join("|")}</output>
    </>
  );
}

describe("gene tag input", () => {
  it("splits a space-separated list typed before Enter into separate genes", async () => {
    render(<Harness />);
    await userEvent.type(screen.getByRole("textbox"), "CIP2A SET ANP32A ARPP19 PABIR1{Enter}");
    expect(screen.getByTestId("genes")).toHaveTextContent("CIP2A|SET|ANP32A|ARPP19|PABIR1");
  });

  it("splits a pasted list with headings and former names, without duplicates", () => {
    render(<Harness initial={["SET"]} />);
    const box = screen.getByRole("textbox");
    fireEvent.paste(box, {
      clipboardData: { getData: () => "PP2A inhibitors: CIP2A (KIAA1524), SET, ANP32A\nPABIR1 (FAM122A)" },
    });
    expect(screen.getByTestId("genes")).toHaveTextContent("SET|CIP2A|ANP32A|PABIR1");
  });

  it("keeps text left in the box when focus moves away", async () => {
    render(
      <>
        <Harness />
        <button type="button">Score</button>
      </>,
    );
    await userEvent.type(screen.getByRole("textbox"), "mycn");
    await userEvent.click(screen.getByRole("button", { name: "Score" }));
    expect(screen.getByTestId("genes")).toHaveTextContent("MYCN");
  });
});
