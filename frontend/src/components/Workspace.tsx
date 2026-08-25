import { useCallback } from "react";
import { DockviewReact, type DockviewApi, type DockviewReadyEvent } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import { PANE_COMPONENTS, PANE_LABELS, type PaneType } from "../panes/registry";

const STORAGE_KEY = "expression-explorer:workspace-layout";
let paneCounter = 0;

// After restoring a saved layout, panel ids like "genome-tracks-1" already
// exist from a PREVIOUS session's counter, but this session's paneCounter
// starts back at 0 -- the next openPane() call for that same type collides
// with the restored id and dockview throws ("panel with id ... already
// exists"), aborting the new pane before its component ever mounts. Seed
// the counter from the highest suffix actually present in the restored
// layout so a fresh session's ids never overlap with what fromJSON just
// recreated.
function seedPaneCounterFromLayout(layout: unknown): void {
  const json = JSON.stringify(layout);
  let max = 0;
  for (const m of json.matchAll(/"id":"[a-z-]+-(\d+)"/g)) {
    max = Math.max(max, Number(m[1]));
  }
  paneCounter = Math.max(paneCounter, max);
}

export function openPane(api: DockviewApi, type: PaneType, options?: { splitRight?: boolean }) {
  // A plain click focuses an existing pane of this type rather than piling
  // up duplicates -- every click on "Signature Score" used to open a brand
  // new tab even with one already open, so three clicks left three
  // identical tabs. Splitting is still always a deliberate new pane: the
  // whole point of alt/cmd-click is placing a second copy side by side
  // with what's already open, so that path is untouched.
  if (!options?.splitRight) {
    const existing = api.panels.find((p) => p.id === type || p.id.startsWith(`${type}-`));
    if (existing) {
      existing.api.setActive();
      return;
    }
  }
  const id = `${type}-${++paneCounter}`;
  api.addPanel({
    id,
    title: PANE_LABELS[type],
    component: type,
    position:
      options?.splitRight && api.activePanel
        ? { referencePanel: api.activePanel.id, direction: "right" }
        : undefined,
  });
}

export function Workspace({ onApiReady }: { onApiReady?: (api: DockviewApi) => void }) {
  const onReady = useCallback(
    (event: DockviewReadyEvent) => {
      onApiReady?.(event.api);

      const saved = localStorage.getItem(STORAGE_KEY);
      if (saved) {
        try {
          const parsed = JSON.parse(saved);
          event.api.fromJSON(parsed);
          seedPaneCounterFromLayout(parsed);
        } catch {
          openPane(event.api, "compare");
        }
      } else {
        openPane(event.api, "compare");
      }

      event.api.onDidLayoutChange(() => {
        localStorage.setItem(STORAGE_KEY, JSON.stringify(event.api.toJSON()));
      });
    },
    [onApiReady],
  );

  return (
    <div className="h-full w-full">
      <DockviewReact className="dockview-theme-expression" components={PANE_COMPONENTS} onReady={onReady} />
    </div>
  );
}

export type { PaneType };
