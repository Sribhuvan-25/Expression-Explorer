/**
 * dockview's default tab nested a close <button> inside the focusable tab
 * (axe: nested-interactive). The custom tab exposes no inner control to
 * assistive tech and closes on Delete / Backspace instead.
 */
import { fireEvent, render, screen } from "@testing-library/react";
import { describe, expect, it, vi } from "vitest";

import { PaneTab } from "../PaneTab";

function mount() {
  const close = vi.fn();
  const api = { title: "Survival", close, onDidTitleChange: () => ({ dispose() {} }) };
  const host = document.createElement("div");
  host.className = "dv-tab";
  host.tabIndex = 0;
  document.body.appendChild(host);
  render(<PaneTab {...({ api } as never)} />, { container: host });
  return { host, close };
}

describe("pane tab", () => {
  it("has no focusable control inside the tab and names how to close it", () => {
    const { host } = mount();
    expect(host.querySelector("button")).toBeNull();
    expect(screen.getByTitle("Close Survival")).toHaveAttribute("aria-hidden", "true");
    expect(host).toHaveAttribute("aria-label", "Survival — press Delete to close");
  });

  it("closes with Delete, Backspace or a click on the ×", () => {
    const { host, close } = mount();
    fireEvent.keyDown(host, { key: "Delete" });
    fireEvent.keyDown(host, { key: "Backspace" });
    fireEvent.click(screen.getByTitle("Close Survival"));
    expect(close).toHaveBeenCalledTimes(3);
  });
});
