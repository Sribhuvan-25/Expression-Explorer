import type { PaneType } from "../panes/registry";

/**
 * Lets a page ask the workspace to open (or focus) another pane without
 * holding a reference to the dockview API. Pages render inside dockview
 * panels and are also rendered standalone in tests, so they cannot assume
 * a workspace exists: with no listener this is simply a no-op.
 */
export const OPEN_PANE_EVENT = "expression-explorer:open-pane";

export function requestOpenPane(type: PaneType): void {
  window.dispatchEvent(new CustomEvent<{ type: PaneType }>(OPEN_PANE_EVENT, { detail: { type } }));
}
