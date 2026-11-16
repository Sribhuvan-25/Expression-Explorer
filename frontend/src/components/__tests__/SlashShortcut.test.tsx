/**
 * "/" jumps to the context gene field. It listened on window in the bubble
 * phase, so a library that stops key events at the document (igv.js does)
 * could swallow it; and inside a shadow root e.target is only the host, so
 * typing "/" in igv.js's own locus box would have jumped away.
 */
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen } from "@testing-library/react";
import { afterEach, describe, expect, it, vi } from "vitest";

import { AppShell } from "../AppShell";
import { api } from "../../lib/api";
import { WorkbenchProvider } from "../../lib/workbench";

function mount() {
  vi.spyOn(api, "listDatasets").mockResolvedValue({ datasets: [] } as never);
  return render(
    <QueryClientProvider client={new QueryClient()}>
      <WorkbenchProvider>
        <AppShell />
      </WorkbenchProvider>
    </QueryClientProvider>,
  );
}
const press = (target: EventTarget) =>
  target.dispatchEvent(new KeyboardEvent("keydown", { key: "/", bubbles: true, composed: true, cancelable: true }));

let swallow: ((e: Event) => void) | null = null;
afterEach(() => {
  if (swallow) document.removeEventListener("keydown", swallow);
  swallow = null;
  document.querySelectorAll("[data-shadow-test]").forEach((n) => n.remove());
});

describe('"/" shortcut', () => {
  it("still reaches the gene field when a document listener stops the event", async () => {
    mount();
    const gene = await screen.findByLabelText("Context genes");
    swallow = (e) => e.stopPropagation();
    document.addEventListener("keydown", swallow);
    press(document.body);
    expect(document.activeElement).toBe(gene);
  });

  it("leaves focus alone while typing in an input inside a shadow root", async () => {
    mount();
    await screen.findByLabelText("Context genes");
    const host = document.createElement("div");
    host.setAttribute("data-shadow-test", "");
    document.body.appendChild(host);
    const input = document.createElement("input");
    host.attachShadow({ mode: "open" }).appendChild(input);
    input.focus();
    press(input);
    expect(host.shadowRoot!.activeElement).toBe(input);
  });
});
