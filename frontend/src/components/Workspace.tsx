import { useCallback } from "react";
import { DockviewReact, type DockviewApi, type DockviewReadyEvent } from "dockview-react";
import "dockview-react/dist/styles/dockview.css";
import { PANE_COMPONENTS, PANE_LABELS, type PaneType } from "../panes/registry";

const STORAGE_KEY = "expression-explorer:workspace-layout";
let paneCounter = 0;

export function openPane(api: DockviewApi, type: PaneType, options?: { splitRight?: boolean }) {
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
          event.api.fromJSON(JSON.parse(saved));
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
