import type { DockviewApi } from "dockview-react";
import { useCallback, useState } from "react";
import { openPane, Workspace } from "./Workspace";
import { Sources } from "./Sources";
import { PANE_LABELS, type PaneType } from "../panes/registry";
import { DatasetStatus } from "./DatasetStatus";

const NAV_ITEMS: { type: PaneType; hint: string }[] = [
  { type: "compare", hint: "Gene across groups" },
  { type: "signature", hint: "Gene set, ranked" },
  { type: "rank", hint: "Pick candidates" },
  { type: "survival", hint: "Kaplan–Meier + Cox" },
  { type: "genome-tracks", hint: "ChIP-seq / ATAC-seq" },
];

function MarkGlyph() {
  return (
    <svg width="22" height="22" viewBox="0 0 22 22" fill="none" aria-hidden="true">
      <path
        d="M4 18C4 18 6 4 11 4C16 4 18 18 18 18"
        stroke="var(--rail-accent)"
        strokeWidth="1.6"
        strokeLinecap="round"
      />
      <path
        d="M4 11C4 11 7 15 11 11C15 7 18 11 18 11"
        stroke="var(--rail-accent)"
        strokeWidth="1.6"
        strokeLinecap="round"
        opacity="0.55"
      />
    </svg>
  );
}

function NavRow({ type, hint, onOpen }: { type: PaneType; hint: string; onOpen: (splitRight: boolean) => void }) {
  return (
    <button
      type="button"
      onClick={(e) => onOpen(e.altKey || e.metaKey)}
      title="Click to open · Alt/Cmd-click to split beside the active pane"
      className="group flex flex-col gap-0.5 rounded-[3px] px-3 py-2.5 text-left text-rail-ink-mute transition-colors hover:bg-rail-active-bg hover:text-rail-ink"
    >
      <span className="flex items-center gap-2 text-[13.5px] font-medium">
        <span className="h-1.5 w-1.5 rounded-full bg-rail-rule transition-colors group-hover:bg-rail-accent" />
        {PANE_LABELS[type]}
      </span>
      <span className="pl-3.5 font-mono text-[10.5px] uppercase tracking-wider text-rail-ink-mute">{hint}</span>
    </button>
  );
}

export function AppShell() {
  const [api, setApi] = useState<DockviewApi | null>(null);

  const handleApiReady = useCallback((readyApi: DockviewApi) => {
    setApi(readyApi);
  }, []);

  return (
    <div className="flex h-screen overflow-hidden bg-ground text-ink">
      <aside className="flex w-[248px] shrink-0 flex-col border-r border-rail-rule bg-rail-bg">
        <div className="flex items-center gap-2.5 border-b border-rail-rule px-4 py-4">
          <MarkGlyph />
          <span className="font-display text-[15px] font-semibold tracking-tight text-rail-ink">
            Expression Explorer
          </span>
        </div>

        <nav className="flex flex-col gap-1 px-2.5 py-3">
          <span className="px-3 pb-1 pt-1 font-mono text-[10px] uppercase tracking-wider text-rail-ink-mute">
            Analyses — click to open, alt-click to split
          </span>
          {NAV_ITEMS.map((item) => (
            <NavRow
              key={item.type}
              {...item}
              onOpen={(splitRight) => api && openPane(api, item.type, { splitRight })}
            />
          ))}
        </nav>

        <div className="flex-1" />

        <DatasetStatus />

        <Sources />
      </aside>

      <main className="min-w-0 flex-1 overflow-hidden">
        <Workspace onApiReady={handleApiReady} />
      </main>
    </div>
  );
}
