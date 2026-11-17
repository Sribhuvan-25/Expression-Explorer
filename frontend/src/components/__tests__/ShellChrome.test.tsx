/**
 * The app shell's chrome, from the researcher-committee review:
 *  - the key hints moved out of the header into a "Keyboard shortcuts"
 *    list, which must stay reachable (button and "?") or the shortcuts
 *    become undiscoverable;
 *  - the context bar says it is optional.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { act, fireEvent, render, screen } from "@testing-library/react";
import { beforeEach, describe, expect, it, vi } from "vitest";

import { AppShell } from "../AppShell";
import { api } from "../../lib/api";

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

describe("keyboard shortcuts list", () => {
  it("opens from its button, lists the shortcuts, and Escape closes it back onto the button", async () => {
    mount();
    const button = await screen.findByRole("button", { name: "Keyboard shortcuts" });
    expect(button).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(button);
    const list = screen.getByRole("dialog", { name: "Keyboard shortcuts" });
    for (const k of ["/", "Delete", "?"]) expect(list).toHaveTextContent(k);
    expect(list).toHaveTextContent(/1–8/);
    key("Escape");
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    expect(document.activeElement).toBe(button);
  });

  it('toggles with "?" but not while typing in a field', async () => {
    mount();
    const gene = await screen.findByLabelText("Context genes");
    key("?");
    expect(screen.getByRole("dialog", { name: "Keyboard shortcuts" })).toBeInTheDocument();
    key("?");
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    gene.focus();
    key("?", gene);
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });
});

describe("context bar", () => {
  it("states that it is optional, in its label and its fields' description", async () => {
    mount();
    const gene = await screen.findByLabelText("Context genes");
    expect(screen.getByRole("search", { name: "Analysis context" })).toHaveTextContent("(optional)");
    expect(gene).toHaveAccessibleDescription(/not required/i);
    expect(screen.getByLabelText("Context dataset")).toHaveAccessibleDescription(/not required/i);
  });
});
