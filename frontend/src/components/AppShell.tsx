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
  { type: "correlation", hint: "Gene vs gene, network, PCA" },
  { type: "genome-tracks", hint: "ChIP-seq / ATAC-seq" },
];

// One small distinct glyph per analysis type so the collapsed, icon-only
// rail (below 1024px — see NavRow) still reads as different destinations
// rather than five identical dots.
const NAV_ICONS: Record<PaneType, React.ReactNode> = {
  compare: (
    <path d="M4 20V10M12 20V4M20 20V14" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
  ),
  signature: (
    <>
      <circle cx="11" cy="11" r="7.5" stroke="currentColor" strokeWidth="1.6" />
      <path d="M11 6.5V11L14 13" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
    </>
  ),
  rank: <path d="M5 6h6M5 11h9M5 16h13" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />,
  survival: <path d="M4 18C9 18 9 6 18 6" stroke="currentColor" strokeWidth="1.6" fill="none" strokeLinecap="round" />,
  correlation: (
    <>
      <circle cx="5" cy="16" r="1.7" fill="currentColor" />
      <circle cx="10" cy="10" r="1.7" fill="currentColor" />
      <circle cx="14" cy="13" r="1.7" fill="currentColor" />
      <circle cx="18" cy="6" r="1.7" fill="currentColor" />
    </>
  ),
  "genome-tracks": (
    <>
      <rect x="3.5" y="9" width="15" height="4" rx="1" stroke="currentColor" strokeWidth="1.4" />
      <path d="M6 9V6M10 9V4M14 9V7M6 13V16M10 13V18M14 13V15" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" />
    </>
  ),
};

function NavIcon({ type }: { type: PaneType }) {
  return (
    <svg width="16" height="16" viewBox="0 0 22 22" fill="none" className="shrink-0">
      {NAV_ICONS[type]}
    </svg>
  );
}

function MarkGlyph() {
  return (
    <svg width="22" height="22" viewBox="0 0 22 22" fill="none" aria-hidden="true" className="shrink-0">
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
      title={`${PANE_LABELS[type]} — click to open, alt/cmd-click to split`}
      className="group flex items-center gap-2.5 rounded-[3px] px-3 py-2.5 text-left text-rail-ink-mute transition-colors hover:bg-rail-active-bg hover:text-rail-ink max-lg:justify-center max-lg:px-0"
    >
      <NavIcon type={type} />
      <span className="flex min-w-0 flex-col gap-0.5 max-lg:hidden">
        <span className="text-[13.5px] font-medium">{PANE_LABELS[type]}</span>
        <span className="font-mono text-[10.5px] uppercase tracking-wider text-rail-ink-mute">{hint}</span>
      </span>
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
      <aside className="flex w-[248px] shrink-0 flex-col border-r border-rail-rule bg-rail-bg max-lg:w-[68px]">
        <div className="flex items-center gap-2.5 border-b border-rail-rule px-4 py-4 max-lg:justify-center max-lg:px-0">
          <MarkGlyph />
          <span className="font-display text-[15px] font-semibold tracking-tight text-rail-ink max-lg:hidden">
            Expression Explorer
          </span>
        </div>

        <nav className="flex flex-col gap-1 px-2.5 py-3">
          <span className="px-3 pb-1 pt-1 font-mono text-[10px] uppercase tracking-wider text-rail-ink-mute max-lg:hidden">
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

        <div className="max-lg:hidden">
          <DatasetStatus />
          <Sources />
        </div>
      </aside>

      <main className="min-w-0 flex-1 overflow-hidden">
        <Workspace onApiReady={handleApiReady} />
      </main>
    </div>
  );
}
