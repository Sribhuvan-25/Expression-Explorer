import type { DockviewApi } from "dockview-react";
import { useCallback, useEffect, useRef, useState } from "react";
import { openPane, Workspace } from "./Workspace";
import { Sources } from "./Sources";
import { PANE_LABELS, type PaneType } from "../panes/registry";
import { DatasetStatus } from "./DatasetStatus";
import { ContextBar } from "./ContextBar";
import { WorkbenchProvider } from "../lib/workbench";

const NAV_ITEMS: { type: PaneType; hint: string }[] = [
  { type: "compare", hint: "Gene across groups" },
  { type: "contrasts", hint: "Designed experiments, per gene" },
  { type: "signature", hint: "Gene set, ranked" },
  { type: "rank", hint: "Pick candidates" },
  { type: "survival", hint: "Kaplan–Meier + Cox" },
  { type: "correlation", hint: "Gene vs gene, network, PCA" },
  { type: "multi-omics", hint: "Mutations, CRISPR, drugs" },
  { type: "genome-tracks", hint: "ChIP-seq / ATAC-seq" },
];

// One small distinct glyph per analysis type so the collapsed, icon-only
// nav (below 1024px) still reads as different destinations.
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
  "multi-omics": (
    <>
      <circle cx="6" cy="6" r="2.2" stroke="currentColor" strokeWidth="1.5" />
      <circle cx="16" cy="6" r="2.2" stroke="currentColor" strokeWidth="1.5" />
      <circle cx="11" cy="16" r="2.2" stroke="currentColor" strokeWidth="1.5" />
      <path d="M7.6 7.7L9.6 14M14.4 7.7L12.4 14M8.2 6h5.6" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" />
    </>
  ),
  contrasts: (
    <path d="M11 3V19M6 7H11M11 11H17M8 15H11" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
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
    <svg width="16" height="16" viewBox="0 0 22 22" fill="none" aria-hidden="true" className="shrink-0">
      {NAV_ICONS[type]}
    </svg>
  );
}

function MarkGlyph() {
  return (
    <svg width="18" height="18" viewBox="0 0 22 22" fill="none" aria-hidden="true" className="shrink-0">
      <path d="M3 19H19" stroke="var(--ink-mute)" strokeWidth="1.4" strokeLinecap="round" />
      <path d="M4 16C7 16 8 5 11 5C14 5 15 16 18 16" stroke="var(--accent)" strokeWidth="1.7" strokeLinecap="round" />
    </svg>
  );
}

// macOS reports the Option key as Alt; show the glyph its users read.
const IS_MAC = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.platform || navigator.userAgent);
const ALT = IS_MAC ? "⌥" : "Alt+";

function isTypingTarget(el: EventTarget | null): boolean {
  if (!(el instanceof HTMLElement)) return false;
  return el.isContentEditable || ["INPUT", "TEXTAREA", "SELECT"].includes(el.tagName);
}

function NavRow({
  type,
  hint,
  index,
  active,
  onOpen,
}: {
  type: PaneType;
  hint: string;
  index: number;
  active: boolean;
  onOpen: (splitRight: boolean) => void;
}) {
  return (
    <button
      type="button"
      data-nav-row
      onClick={(e) => onOpen(e.altKey || e.metaKey)}
      aria-current={active ? "page" : undefined}
      aria-keyshortcuts={`Alt+${index + 1}`}
      aria-label={`${PANE_LABELS[type]} — ${hint}`}
      title={`${PANE_LABELS[type]} — click to open, alt/cmd-click to split (${ALT}${index + 1})`}
      className={[
        "relative flex w-full items-center gap-2.5 rounded-[3px] px-2.5 py-1.5 text-left transition-colors max-lg:justify-center max-lg:px-0 max-lg:py-2",
        active ? "bg-rail-active-bg text-rail-active-ink" : "text-ink-soft hover:bg-surface-hover hover:text-ink",
      ].join(" ")}
    >
      {active && <span aria-hidden="true" className="absolute inset-y-1 left-0 w-[2px] rounded-full bg-accent" />}
      <NavIcon type={type} />
      <span className="flex min-w-0 flex-1 flex-col max-lg:hidden">
        <span className={`truncate text-[13px] ${active ? "font-semibold" : "font-medium"}`}>{PANE_LABELS[type]}</span>
        <span className="truncate text-[11px] text-ink-mute">{hint}</span>
      </span>
      <span aria-hidden="true" className="wb-kbd max-lg:hidden">
        {ALT}
        {index + 1}
      </span>
    </button>
  );
}

function focusPaneHeading(api: DockviewApi) {
  const id = api.activePanel?.id;
  if (!id) return;
  // Two frames: dockview mounts the panel's React portal after addPanel
  // returns, and the heading only exists once that has rendered.
  requestAnimationFrame(() =>
    requestAnimationFrame(() => {
      document.querySelector<HTMLElement>(`[data-pane-id="${id}"] [data-pane-heading]`)?.focus();
    }),
  );
}

export function AppShell() {
  return (
    <WorkbenchProvider>
      <Shell />
    </WorkbenchProvider>
  );
}

function Shell() {
  const [api, setApi] = useState<DockviewApi | null>(null);
  const [activeType, setActiveType] = useState<PaneType | null>(null);
  const contextGeneRef = useRef<HTMLInputElement>(null);
  const navRef = useRef<HTMLElement>(null);

  const handleApiReady = useCallback((readyApi: DockviewApi) => {
    setApi(readyApi);
  }, []);

  // Which analysis is in front, for the nav's current marker.
  useEffect(() => {
    if (!api) return;
    const sync = () => {
      const id = api.activePanel?.id ?? "";
      setActiveType(NAV_ITEMS.find((n) => id === n.type || id.startsWith(`${n.type}-`))?.type ?? null);
    };
    sync();
    const sub = api.onDidActivePanelChange(sync);
    return () => sub.dispose();
  }, [api]);

  const open = useCallback(
    (type: PaneType, splitRight = false, focus = false) => {
      if (!api) return;
      openPane(api, type, { splitRight });
      if (focus) focusPaneHeading(api);
    },
    [api],
  );

  // Keyboard: Alt+1..8 opens (or brings forward) an analysis; "/" jumps
  // to the context gene field from anywhere you are not typing. e.code,
  // not e.key: on macOS Option+1 types "¡".
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      if (e.altKey && !e.ctrlKey && !e.metaKey && /^Digit[1-8]$/.test(e.code)) {
        const item = NAV_ITEMS[Number(e.code.slice(5)) - 1];
        if (item) {
          e.preventDefault();
          open(item.type, false, true);
        }
        return;
      }
      if (e.key === "/" && !e.altKey && !e.ctrlKey && !e.metaKey && !isTypingTarget(e.target)) {
        e.preventDefault();
        contextGeneRef.current?.focus();
        contextGeneRef.current?.select();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [open]);

  // Arrow keys move between analyses in the nav list.
  const onNavKey = (e: React.KeyboardEvent) => {
    if (!["ArrowDown", "ArrowUp", "Home", "End"].includes(e.key)) return;
    const rows = Array.from(navRef.current?.querySelectorAll<HTMLElement>("[data-nav-row]") ?? []);
    const i = rows.indexOf(document.activeElement as HTMLElement);
    if (i < 0) return;
    e.preventDefault();
    const next =
      e.key === "Home"
        ? 0
        : e.key === "End"
          ? rows.length - 1
          : (i + (e.key === "ArrowDown" ? 1 : -1) + rows.length) % rows.length;
    rows[next]?.focus();
  };

  return (
    <div className="flex h-screen flex-col overflow-hidden bg-ground text-ink">
      <header className="flex shrink-0 flex-wrap items-center gap-x-4 gap-y-2 border-b border-rule-firm bg-surface px-3 py-2">
        <div className="flex w-[208px] shrink-0 items-center gap-2 max-lg:w-auto">
          <MarkGlyph />
          <span className="text-[13.5px] font-semibold text-ink">Expression Explorer</span>
        </div>
        <ContextBar ref={contextGeneRef} />
        <p className="flex shrink-0 items-center gap-1.5 text-[11px] text-ink-mute max-xl:hidden">
          <span className="wb-kbd">/</span> context gene
          <span className="wb-kbd ml-2">{ALT}1–8</span> analyses
        </p>
      </header>

      <div className="flex min-h-0 flex-1">
        {/* overflow-y-auto: the dataset list below the nav runs past short
            screens, so the column scrolls rather than clipping it. */}
        <aside className="flex w-[232px] shrink-0 flex-col overflow-y-auto border-r border-rule bg-rail-bg max-lg:w-[48px]">
          <nav
            ref={navRef}
            aria-label="Analyses"
            onKeyDown={onNavKey}
            className="flex flex-col gap-px px-1.5 py-2 max-lg:px-1"
          >
            <p className="px-2.5 pb-1 text-[11px] font-semibold text-ink-mute max-lg:sr-only">
              Analyses <span className="font-normal">— alt-click to split</span>
            </p>
            {NAV_ITEMS.map((item, i) => (
              <NavRow
                key={item.type}
                {...item}
                index={i}
                active={activeType === item.type}
                onOpen={(splitRight) => open(item.type, splitRight)}
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
    </div>
  );
}
