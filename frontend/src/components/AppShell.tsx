import type { DockviewApi } from "dockview-react";
import { useCallback, useEffect, useRef, useState } from "react";
import { flushSync } from "react-dom";
import { openPane, Workspace } from "./Workspace";
import { Sources } from "./Sources";
import { PANE_LABELS, type PaneType } from "../panes/registry";
import { DatasetStatus } from "./DatasetStatus";
import { ContextBar } from "./ContextBar";
import { useWorkbench, WorkbenchProvider } from "../lib/workbench";

// Each analysis is introduced by the question a biologist would bring to
// it (the review's first-year student could not tell which analysis
// answered which question from method names alone); the name stays the
// row's title and the method sits in its tooltip and description.
// `short` labels the icon in the narrow rail (640–1023px).
const NAV_ITEMS: { type: PaneType; short: string; question: string; method: string }[] = [
  { type: "compare", short: "Compare", question: "Does a gene differ by group?", method: "Box plot, Mann–Whitney + FDR; top differential genes" },
  { type: "contrasts", short: "Contrasts", question: "What does a perturbation change?", method: "Precomputed limma contrasts, per gene" },
  { type: "signature", short: "Signature", question: "How active is a gene set?", method: "AUCell or log2-mean signature score per sample" },
  { type: "rank", short: "Ranking", question: "Which samples express it most?", method: "Samples sorted by a gene or signature" },
  { type: "survival", short: "Survival", question: "Does it predict survival?", method: "Kaplan–Meier, log-rank + Cox" },
  { type: "correlation", short: "Correlation", question: "Which genes move together?", method: "Gene vs gene, co-expression network, PCA" },
  { type: "multi-omics", short: "Multi-omics", question: "Tied to mutations, CRISPR, drugs?", method: "Mutations, CRISPR dependency, drug response" },
  { type: "genome-tracks", short: "Tracks", question: "Is a region bound or open?", method: "ChIP-seq / ATAC-seq signal tracks" },
];

// One small distinct glyph per analysis type; in the narrow rail
// (640–1023px) each sits above its short name.
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

const SHORTCUTS: { keys: string; does: string }[] = [
  { keys: "/", does: "Jump to the context gene field" },
  { keys: `${ALT}1–8`, does: "Open an analysis (in list order)" },
  { keys: `${ALT}-click`, does: "Open an analysis side by side" },
  { keys: "↑ ↓", does: "Move between analyses in the list" },
  { keys: "Delete", does: "Close the focused tab" },
  { keys: "Esc", does: "Undo an unapplied context edit" },
  { keys: "?", does: "Show or hide this list" },
];

/**
 * The keyboard shortcuts, behind one quiet "?" button. They used to sit
 * as key hints in the header and on every nav row, which first-time
 * users read as power-user chrome they had to understand first; the
 * shortcuts themselves are unchanged.
 */
function ShortcutsHelp({ open, setOpen }: { open: boolean; setOpen: (open: boolean) => void }) {
  const rootRef = useRef<HTMLDivElement>(null);
  const buttonRef = useRef<HTMLButtonElement>(null);

  useEffect(() => {
    if (!open) return;
    const onPointer = (e: PointerEvent) => {
      if (!rootRef.current?.contains(e.target as Node)) setOpen(false);
    };
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") {
        setOpen(false);
        buttonRef.current?.focus();
      }
    };
    document.addEventListener("pointerdown", onPointer);
    document.addEventListener("keydown", onKey);
    return () => {
      document.removeEventListener("pointerdown", onPointer);
      document.removeEventListener("keydown", onKey);
    };
  }, [open, setOpen]);

  return (
    <div ref={rootRef} className="relative ml-auto shrink-0 max-sm:hidden">
      <button
        ref={buttonRef}
        type="button"
        aria-expanded={open}
        aria-controls="shortcuts-help"
        aria-keyshortcuts="?"
        onClick={() => setOpen(!open)}
        title="Keyboard shortcuts (?)"
        className="flex h-[26px] items-center gap-1.5 rounded-[3px] px-1.5 text-[11.5px] text-ink-mute hover:bg-surface-hover hover:text-ink"
      >
        <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden="true">
          <rect x="1.5" y="4" width="13" height="8.5" rx="1.5" stroke="currentColor" strokeWidth="1.2" />
          <path d="M4 7h1M7.5 7h1M11 7h1M5 10h6" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" />
        </svg>
        <span className="max-lg:sr-only">Keyboard shortcuts</span>
      </button>
      {open && (
        <div
          id="shortcuts-help"
          role="dialog"
          aria-labelledby="shortcuts-help-title"
          className="absolute right-0 top-[calc(100%+6px)] z-50 w-[300px] rounded-[3px] border border-rule-firm bg-surface p-3 shadow-lg"
        >
          <p id="shortcuts-help-title" className="mb-2 text-[12px] font-semibold text-ink">
            Keyboard shortcuts
          </p>
          <dl className="grid grid-cols-[auto_1fr] items-center gap-x-3 gap-y-1.5 text-[12px]">
            {SHORTCUTS.map((s) => (
              <div key={s.keys} className="contents">
                <dt>
                  <span className="wb-kbd">{s.keys}</span>
                </dt>
                <dd className="m-0 text-ink-soft">{s.does}</dd>
              </div>
            ))}
          </dl>
        </div>
      )}
    </div>
  );
}

function NavRow({
  type,
  short,
  question,
  method,
  index,
  active,
  onOpen,
}: {
  type: PaneType;
  short: string;
  question: string;
  method: string;
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
      aria-label={`${PANE_LABELS[type]} — ${question}`}
      aria-describedby={`nav-method-${type}`}
      title={`${PANE_LABELS[type]}: ${method}. Click to open, alt/cmd-click to split (${ALT}${index + 1})`}
      className={[
        // sm:max-lg = the narrow rail (icon over a short name); below sm
        // the row is shown full width in the nav drawer.
        "group relative flex w-full items-center gap-2.5 rounded-[3px] px-2.5 py-1.5 text-left transition-colors sm:max-lg:flex-col sm:max-lg:gap-1 sm:max-lg:px-0.5 sm:max-lg:py-2 max-sm:py-2.5",
        active ? "bg-rail-active-bg text-rail-active-ink" : "text-ink-soft hover:bg-surface-hover hover:text-ink",
      ].join(" ")}
    >
      {active && <span aria-hidden="true" className="absolute inset-y-1 left-0 w-[2px] rounded-full bg-accent" />}
      <NavIcon type={type} />
      <span className="flex min-w-0 flex-1 flex-col sm:max-lg:hidden">
        <span className={`truncate text-[13px] ${active ? "font-semibold" : "font-medium"}`}>{PANE_LABELS[type]}</span>
        <span className="text-[11px] leading-snug text-ink-mute">{question}</span>
      </span>
      <span aria-hidden="true" className="hidden text-center text-[10.5px] leading-tight sm:max-lg:block">
        {short}
      </span>
      <span id={`nav-method-${type}`} hidden>
        {method}
      </span>
      <span
        aria-hidden="true"
        className="wb-kbd absolute right-2 top-2 opacity-0 transition-opacity group-hover:opacity-100 group-focus-visible:opacity-100 max-lg:hidden"
      >
        {ALT}
        {index + 1}
      </span>
    </button>
  );
}

const CUE_KEY = "expression-explorer:first-run-cue-dismissed";

function cueDismissed(): boolean {
  try {
    return localStorage.getItem(CUE_KEY) === "1";
  } catch {
    return false;
  }
}

/**
 * One line, shown until dismissed or until the researcher opens an
 * analysis from the nav: it points a first-time user at the question-led
 * Analyses list and says the context bar is optional. Deliberately not a
 * tour -- experienced users dismiss it once and never see it again.
 */
function FirstRunCue({ onDismiss }: { onDismiss: () => void }) {
  return (
    <div className="flex shrink-0 items-center gap-2.5 border-b border-rule bg-accent-soft px-3 py-1.5 text-[12.5px] text-accent-ink">
      <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden="true" className="shrink-0">
        <path d="M13 8H3M7 4L3 8L7 12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" />
      </svg>
      <p className="min-w-0 flex-1">
        <span className="font-semibold">New here?</span>{" "}
        <span className="max-sm:hidden">
          Each analysis on the left answers one question — pick the one closest to yours. Setting a context at the top is
          optional.
        </span>
        <span className="sm:hidden">Tap Analyses to see the question each one answers. Context is optional.</span>
      </p>
      <button
        type="button"
        onClick={onDismiss}
        className="h-[24px] shrink-0 rounded-[3px] border border-accent px-2 text-[12px] font-medium text-accent-ink hover:bg-surface"
      >
        Got it
      </button>
    </div>
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
  const [helpOpen, setHelpOpen] = useState(false);
  const [cueOpen, setCueOpen] = useState(() => !cueDismissed());
  const dismissCue = useCallback(() => {
    setCueOpen(false);
    try {
      localStorage.setItem(CUE_KEY, "1");
    } catch {
      /* blocked storage: the cue simply shows again next visit */
    }
  }, []);
  // Below 640px the context bar folds behind a "Context" button and the
  // nav becomes a drawer; both flags only change classes under max-sm.
  const [contextOpen, setContextOpen] = useState(false);
  const [navOpen, setNavOpen] = useState(false);
  const contextGeneRef = useRef<HTMLInputElement>(null);
  const navRef = useRef<HTMLElement>(null);
  const asideRef = useRef<HTMLElement>(null);
  const menuButtonRef = useRef<HTMLButtonElement>(null);
  const contextGenes = useWorkbench()?.genes ?? [];

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
      // They found the nav: the first-run cue has done its job.
      dismissCue();
    },
    [api, dismissCue],
  );

  // Keyboard: Alt+1..8 opens (or brings forward) an analysis; "/" jumps
  // to the context gene field and "?" lists the shortcuts, from anywhere
  // you are not typing. e.code,
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
      // composedPath()[0] is the element really focused: inside igv.js's
      // shadow DOM e.target is only its host, so typing "/" in igv's own
      // locus box would have jumped to the gene field.
      if (e.altKey || e.ctrlKey || e.metaKey || isTypingTarget(e.composedPath()[0] ?? e.target)) return;
      if (e.key === "/") {
        e.preventDefault();
        // On a phone the field may be folded away (display: none, which
        // cannot take focus): unfold it synchronously, then focus.
        flushSync(() => {
          setContextOpen(true);
          setNavOpen(false);
        });
        contextGeneRef.current?.focus();
        contextGeneRef.current?.select();
      } else if (e.key === "?") {
        e.preventDefault();
        setHelpOpen((o) => !o);
      }
    };
    // Capture phase: igv.js listens for keys on the document and can stop
    // them before a bubbling window listener sees them, which made "/"
    // intermittently do nothing once Genome Tracks had loaded.
    window.addEventListener("keydown", onKey, true);
    return () => window.removeEventListener("keydown", onKey, true);
  }, [open]);

  // Opening the drawer moves focus into it, onto the analysis in front.
  useEffect(() => {
    if (!navOpen) return;
    const rows = Array.from(navRef.current?.querySelectorAll<HTMLElement>("[data-nav-row]") ?? []);
    (rows.find((r) => r.getAttribute("aria-current") === "page") ?? rows[0])?.focus();
  }, [navOpen]);

  const closeNav = (refocus: boolean) => {
    setNavOpen(false);
    if (refocus) menuButtonRef.current?.focus();
  };

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
      <header className="flex shrink-0 flex-wrap items-center gap-x-4 gap-y-2 border-b border-rule-firm bg-surface px-3 py-2 max-sm:gap-x-2 max-sm:py-1.5">
        <button
          ref={menuButtonRef}
          type="button"
          aria-expanded={navOpen}
          aria-controls="app-nav"
          onClick={() => setNavOpen((o) => !o)}
          className="flex h-[30px] shrink-0 items-center gap-1.5 rounded-[3px] border border-field-edge px-2 text-[12.5px] font-medium text-ink-soft sm:hidden"
        >
          <svg width="14" height="14" viewBox="0 0 16 16" fill="none" aria-hidden="true">
            <path d="M2 4h12M2 8h12M2 12h12" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" />
          </svg>
          Analyses
        </button>
        <div className="flex w-[208px] shrink-0 items-center gap-2 max-lg:w-auto max-sm:min-w-0 max-sm:flex-1 max-sm:shrink">
          <MarkGlyph />
          <span className="truncate text-[13.5px] font-semibold text-ink">Expression Explorer</span>
        </div>
        <button
          type="button"
          aria-expanded={contextOpen}
          aria-controls="context-region"
          onClick={() => setContextOpen((o) => !o)}
          className="flex h-[30px] min-w-0 shrink-0 items-center gap-1.5 rounded-[3px] border border-field-edge px-2 text-[12.5px] font-medium text-ink-soft sm:hidden"
        >
          Context
          {contextGenes.length > 0 && (
            <span className="wb-gene max-w-[9ch] truncate text-[12px] font-normal text-accent-ink">
              {contextGenes[0]}
              {contextGenes.length > 1 ? ` +${contextGenes.length - 1}` : ""}
            </span>
          )}
          <svg width="10" height="10" viewBox="0 0 10 10" fill="none" aria-hidden="true" className={contextOpen ? "rotate-180" : ""}>
            <path d="M2 3.5L5 6.5L8 3.5" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
        </button>
        <div id="context-region" className={`flex min-w-0 flex-1 ${contextOpen ? "max-sm:basis-full" : "max-sm:hidden"}`}>
          <ContextBar ref={contextGeneRef} />
        </div>
        <ShortcutsHelp open={helpOpen} setOpen={setHelpOpen} />
      </header>

      <div className="flex min-h-0 flex-1">
        {navOpen && (
          <div aria-hidden="true" onClick={() => closeNav(false)} className="fixed inset-0 z-30 bg-black/30 sm:hidden" />
        )}
        {/* overflow-y-auto: the dataset list below the nav runs past short
            screens, so the column scrolls rather than clipping it. Below
            640px the column is a drawer opened from the header. */}
        <aside
          ref={asideRef}
          id="app-nav"
          onKeyDown={(e) => {
            if (e.key === "Escape" && navOpen) {
              e.stopPropagation();
              closeNav(true);
            }
          }}
          onBlur={(e) => {
            // Tabbing out of the open drawer closes it. A null
            // relatedTarget (a tap on iOS, which does not focus buttons)
            // must not, or the tap that follows would land on nothing.
            const to = e.relatedTarget as Node | null;
            if (navOpen && to && !asideRef.current?.contains(to) && to !== menuButtonRef.current) setNavOpen(false);
          }}
          className={[
            "flex w-[232px] shrink-0 flex-col overflow-y-auto border-r border-rule bg-rail-bg sm:max-lg:w-[76px]",
            "max-sm:fixed max-sm:inset-y-0 max-sm:left-0 max-sm:z-40 max-sm:w-[min(300px,86vw)] max-sm:shadow-xl",
            navOpen ? "" : "max-sm:hidden",
          ].join(" ")}
        >
          <div className="flex items-center justify-between border-b border-rule px-3 py-2 sm:hidden">
            <span className="text-[13px] font-semibold text-ink">Analyses</span>
            <button
              type="button"
              onClick={() => closeNav(true)}
              aria-label="Close analyses"
              className="flex h-[30px] w-[30px] items-center justify-center rounded-[3px] text-ink-mute hover:bg-surface-hover hover:text-ink"
            >
              <svg width="11" height="11" viewBox="0 0 10 10" aria-hidden="true">
                <path d="M1 1L9 9M9 1L1 9" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" />
              </svg>
            </button>
          </div>
          <nav
            ref={navRef}
            aria-label="Analyses"
            onKeyDown={onNavKey}
            className="flex flex-col gap-px px-1.5 py-2 sm:max-lg:px-1"
          >
            <p className="px-2.5 pb-1 text-[11px] font-semibold text-ink-mute sm:max-lg:sr-only max-sm:hidden">Analyses</p>
            {NAV_ITEMS.map((item, i) => (
              <NavRow
                key={item.type}
                {...item}
                index={i}
                active={activeType === item.type}
                onOpen={(splitRight) => {
                  if (navOpen) {
                    // From the phone drawer: close it and move focus to the
                    // pane that opened, not back to a now-hidden row.
                    setNavOpen(false);
                    open(item.type, false, true);
                  } else {
                    open(item.type, splitRight);
                  }
                }}
              />
            ))}
          </nav>

          <div className="flex-1" />

          <div className="sm:max-lg:hidden">
            <DatasetStatus />
            <Sources />
          </div>
        </aside>

        <main className="flex min-w-0 flex-1 flex-col overflow-hidden">
          {cueOpen && <FirstRunCue onDismiss={dismissCue} />}
          <div className="min-h-0 flex-1">
            <Workspace onApiReady={handleApiReady} />
          </div>
        </main>
      </div>
    </div>
  );
}
