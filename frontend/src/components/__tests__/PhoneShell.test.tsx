/**
 * On a phone (below 640px) the context bar folds behind a "Context"
 * button and the nav is a drawer opened from "Analyses" (the review found
 * the bar took a third of a 390px screen and the icon-only nav had no
 * labels). jsdom applies no CSS, so these tests pin the state the CSS
 * keys off: the disclosure/drawer flags exposed as aria-expanded.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { act, fireEvent, render, screen, waitFor, within } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { AppShell } from "../AppShell";
import { api } from "../../lib/api";
import { PANE_LABELS } from "../../panes/registry";

function mount() {
  vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [] } as never);
  return render(
    <QueryClientProvider client={new QueryClient()}>
      <AppShell />
    </QueryClientProvider>,
  );
}
const key = (k: string, target: EventTarget = document.body) =>
  act(() => {
    target.dispatchEvent(new KeyboardEvent("keydown", { key: k, bubbles: true, composed: true, cancelable: true }));
  });

beforeEach(() => localStorage.clear());

describe("folded context bar", () => {
  it('"/" unfolds the context bar before focusing the gene field', async () => {
    mount();
    const gene = await screen.findByLabelText("Context genes");
    const toggle = screen.getByRole("button", { name: /^Context/ });
    expect(toggle).toHaveAttribute("aria-expanded", "false");
    key("/");
    // Unfolded, or the field is display:none on a phone and focus is lost.
    expect(toggle).toHaveAttribute("aria-expanded", "true");
    expect(document.activeElement).toBe(gene);
  });
});

describe("phone nav drawer", () => {
  it("lists every analysis by name and question, and closes onto the pane picked", async () => {
    mount();
    const menu = await screen.findByRole("button", { name: "Analyses" });
    expect(menu).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(menu);
    expect(menu).toHaveAttribute("aria-expanded", "true");

    const nav = screen.getByRole("navigation", { name: "Analyses" });
    const rows = within(nav).getAllByRole("button");
    expect(rows).toHaveLength(8);
    for (const type of Object.keys(PANE_LABELS) as (keyof typeof PANE_LABELS)[]) {
      expect(rows.some((r) => r.textContent?.includes(PANE_LABELS[type]))).toBe(true);
    }
    await waitFor(() => expect(nav).toContainElement(document.activeElement as HTMLElement));

    const survival = within(nav).getByRole("button", { name: /^Survival/ });
    fireEvent.click(survival);
    expect(menu).toHaveAttribute("aria-expanded", "false");
    await waitFor(() => expect(survival).toHaveAttribute("aria-current", "page"));
  });

  it("closes on Escape and returns focus to the Analyses button", async () => {
    mount();
    const menu = await screen.findByRole("button", { name: "Analyses" });
    fireEvent.click(menu);
    const nav = screen.getByRole("navigation", { name: "Analyses" });
    await waitFor(() => expect(nav).toContainElement(document.activeElement as HTMLElement));
    fireEvent.keyDown(document.activeElement!, { key: "Escape" });
    expect(menu).toHaveAttribute("aria-expanded", "false");
    expect(document.activeElement).toBe(menu);
  });
});
