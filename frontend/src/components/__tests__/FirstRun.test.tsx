/**
 * The review's first-year student could not tell which analysis answers
 * which question: the nav named methods only and nothing on first load
 * pointed at it. Each nav row now leads with a research question, and a
 * one-line cue points at the nav until dismissed -- once, not on every
 * visit.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { fireEvent, render, screen, within } from "@testing-library/react";
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

beforeEach(() => localStorage.clear());

describe("question-led nav", () => {
  it("shows a research question under every analysis name", async () => {
    mount();
    const nav = await screen.findByRole("navigation", { name: "Analyses" });
    const rows = within(nav).getAllByRole("button").filter((b) => b.hasAttribute("data-nav-row"));
    expect(rows).toHaveLength(8);
    for (const row of rows) expect(row.textContent).toMatch(/\?/);
    expect(within(nav).getByRole("button", { name: /^Survival — Does it predict survival\?$/ })).toHaveAccessibleDescription(
      /Kaplan–Meier/,
    );
  });
});

describe("first-run cue", () => {
  it("shows on a first visit and stays dismissed after Got it", async () => {
    const { unmount } = mount();
    fireEvent.click(await screen.findByRole("button", { name: "Got it" }));
    expect(screen.queryByText("New here?")).toBeNull();
    unmount();
    mount();
    await screen.findByLabelText("Context genes");
    expect(screen.queryByText("New here?")).toBeNull();
  });

  it("goes away once an analysis is opened from the nav", async () => {
    mount();
    expect(await screen.findByText("New here?")).toBeInTheDocument();
    const nav = screen.getByRole("navigation", { name: "Analyses" });
    fireEvent.click(within(nav).getByRole("button", { name: /^Survival/ }));
    expect(screen.queryByText("New here?")).toBeNull();
  });
});
