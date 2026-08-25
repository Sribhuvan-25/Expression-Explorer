import { useCallback, useEffect, useRef } from "react";
import { DockviewReact, type DockviewApi, type DockviewReadyEvent, type IWatermarkPanelProps } from "dockview-react";
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

// Shown by dockview itself whenever no pane is open -- closing the last
// tab used to leave a genuinely blank dark rectangle with no way to tell
// whether the app had crashed or was just empty, and no visible path back
// in short of clicking the sidebar (which still works, but nothing on
// screen said so).
function EmptyWorkspace({ containerApi }: IWatermarkPanelProps) {
  return (
    <div className="flex h-full w-full flex-col items-center justify-center gap-3 text-center">
      <svg width="28" height="28" viewBox="0 0 24 24" fill="none" className="text-ink-mute">
        <path d="M4 20V10M12 20V4M20 20V14" stroke="currentColor" strokeWidth="1.7" strokeLinecap="round" />
      </svg>
      <div>
        <p className="text-[13.5px] font-medium text-ink-soft">No analysis open</p>
        <p className="mt-1 max-w-[36ch] text-[12.5px] text-ink-mute">
          Choose one from the sidebar to get started.
        </p>
      </div>
      <button
        type="button"
        onClick={() => openPane(containerApi, "compare")}
        className="mt-1 rounded-[3px] border border-rule-firm bg-surface px-3 py-1.5 text-[12.5px] text-ink-soft transition-colors hover:border-accent hover:text-accent"
      >
        Open {PANE_LABELS.compare}
      </button>
    </div>
  );
}

export function Workspace({ onApiReady }: { onApiReady?: (api: DockviewApi) => void }) {
  const rootRef = useRef<HTMLDivElement>(null);
  const apiRef = useRef<DockviewApi | null>(null);

  const onReady = useCallback(
    (event: DockviewReadyEvent) => {
      apiRef.current = event.api;
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

  // dockview-react calls api.layout() exactly once, on mount -- there is
  // no ResizeObserver anywhere in the library watching its own root
  // element for later size changes; that's left entirely to the host.
  // Without this, resizing the actual browser window (not a dockview
  // pane drag, which the library does handle internally) left every pane
  // latched at its mount-time width: charts stayed the old width inside a
  // now-narrower ancestor, clipped with no scrollbar to reach the excess
  // since the overflow was absorbed by intermediate `overflow-x: visible`
  // containers rather than surfacing at the page level.
  //
  // Both a ResizeObserver on the root element AND a window 'resize'
  // listener call layout(), not just one: they catch different cases.
  // ResizeObserver fires for any size change of this element regardless
  // of cause (e.g. the sidebar's own width changing at a breakpoint,
  // independent of the window); 'resize' is the reliable signal for an
  // actual browser-window resize specifically. Relying on only the
  // observer left a real gap -- some resize paths (a window resize
  // delivered as a single native reflow rather than the OS's usual
  // incremental resize-drag events) did not reliably produce a
  // ResizeObserver callback in testing, while 'resize' always fired.
  useEffect(() => {
    if (!rootRef.current) return;
    const el = rootRef.current;
    const relayout = () => {
      if (!apiRef.current || !rootRef.current) return;
      const { clientWidth, clientHeight } = rootRef.current;
      apiRef.current.layout(clientWidth, clientHeight);
    };
    const observer = new ResizeObserver(relayout);
    observer.observe(el);
    window.addEventListener("resize", relayout);
    return () => {
      observer.disconnect();
      window.removeEventListener("resize", relayout);
    };
  }, []);

  return (
    <div ref={rootRef} className="h-full w-full">
      <DockviewReact
        className="dockview-theme-expression"
        components={PANE_COMPONENTS}
        watermarkComponent={EmptyWorkspace}
        onReady={onReady}
      />
    </div>
  );
}

export type { PaneType };
