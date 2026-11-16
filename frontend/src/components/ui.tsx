import { useEffect, useId, useRef, useState, type ReactNode } from "react";
import { parseGeneInput } from "../lib/geneInput";
import type { LinkState } from "../lib/workbench";

/**
 * Compact page header. The pane tab already names the analysis, so the
 * header spends one line on orientation ("Expression Compare / what it
 * does") and one on method, then gets out of the way of the result.
 * The h1 takes programmatic focus (tabIndex -1) so keyboard navigation
 * between analyses lands somewhere meaningful.
 */
export function PageHeader({
  eyebrow,
  title,
  description,
}: {
  eyebrow: string;
  title: string;
  description: string;
}) {
  return (
    <header className="border-b border-rule bg-surface px-5 py-3 max-sm:px-3">
      <div className="flex flex-wrap items-baseline gap-x-2 gap-y-0.5">
        <p className="text-[12px] font-medium text-ink-mute">
          {eyebrow}
          <span aria-hidden="true" className="px-1.5 text-rule-firm">/</span>
        </p>
        <h1 tabIndex={-1} data-pane-heading className="text-[15px] font-semibold leading-tight text-ink outline-none">
          {title}
        </h1>
      </div>
      <p className="mt-1 max-w-[100ch] text-[12.5px] leading-snug text-ink-mute">{description}</p>
    </header>
  );
}

export function Panel({
  title,
  action,
  children,
  className = "",
}: {
  title?: string;
  action?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  return (
    <section className={`rounded-[3px] border border-rule bg-surface ${className}`}>
      {title && (
        // flex-wrap: on a phone-width pane the actions (unit, n, PNG) were
        // squeezed into narrow columns beside the title; they now drop below.
        <div className="flex min-h-[34px] flex-wrap items-center justify-between gap-x-3 gap-y-1 border-b border-rule px-3 py-1.5">
          <h2 className="min-w-0 text-[12.5px] font-semibold text-ink">{title}</h2>
          {action && <div className="max-w-full [&>div]:flex-wrap">{action}</div>}
        </div>
      )}
      <div className="p-3">{children}</div>
    </section>
  );
}

/**
 * The query panel of an analysis. Once a run has produced a result it
 * collapses to a one-line summary of what ran, with an Edit control, so
 * the result gets the screen. The controls stay mounted while collapsed
 * (hidden with CSS) -- their state, and anything that reads it, is
 * untouched. It stays open while `dirty` (the live inputs no longer
 * match the run, e.g. the context gene changed), so the change is never
 * hidden behind a summary of an older query.
 */
export function QueryPanel({
  title = "Query",
  summary,
  done,
  runKey,
  dirty = false,
  children,
}: {
  title?: string;
  summary: ReactNode;
  done: boolean;
  runKey: unknown;
  dirty?: boolean;
  children: ReactNode;
}) {
  const bodyId = useId();
  const [editing, setEditing] = useState(false);
  const sectionRef = useRef<HTMLElement>(null);
  const editRef = useRef<HTMLButtonElement>(null);
  useEffect(() => setEditing(false), [runKey]);
  // Once a change forces the panel open, keep it open until the next run
  // -- otherwise typing the old value back would snap it shut mid-edit.
  useEffect(() => {
    if (dirty && done) setEditing(true);
  }, [dirty, done]);
  const collapsed = done && !editing && !dirty;

  // Collapsing hides whatever had focus (the Run button, or the input the
  // user pressed Enter in). Hand focus to Edit rather than dropping it on
  // <body>, where a keyboard user would have to start over.
  const wasCollapsed = useRef(collapsed);
  useEffect(() => {
    if (collapsed && !wasCollapsed.current) {
      const active = document.activeElement;
      if (!active || active === document.body || sectionRef.current?.contains(active)) editRef.current?.focus();
    }
    if (!collapsed && wasCollapsed.current && editing) {
      sectionRef.current?.querySelector<HTMLElement>(`[id="${bodyId}"] input, [id="${bodyId}"] select, [id="${bodyId}"] button`)?.focus();
    }
    wasCollapsed.current = collapsed;
  }, [collapsed, editing, bodyId]);

  return (
    <section ref={sectionRef} className="rounded-[3px] border border-rule bg-surface">
      {collapsed ? (
        <div className="flex min-h-[38px] flex-wrap items-center gap-x-3 gap-y-1 px-3 py-1.5">
          <h2 className="text-[12.5px] font-semibold text-ink">{title}</h2>
          <div data-testid="query-summary" className="min-w-0 flex-1 text-[12.5px] text-ink-soft">
            {summary}
          </div>
          <button
            ref={editRef}
            type="button"
            aria-expanded={false}
            aria-controls={bodyId}
            onClick={() => setEditing(true)}
            className="wb-btn-quiet"
          >
            Edit query
          </button>
        </div>
      ) : (
        <div className="flex min-h-[34px] flex-wrap items-center justify-between gap-x-3 gap-y-1 border-b border-rule px-3 py-1.5">
          <h2 className="text-[12.5px] font-semibold text-ink">{title}</h2>
          {done && !dirty && (
            <button
              type="button"
              aria-expanded={true}
              aria-controls={bodyId}
              onClick={() => setEditing(false)}
              className="text-[11.5px] font-medium text-ink-mute underline decoration-rule-firm underline-offset-2 hover:text-accent"
            >
              Collapse
            </button>
          )}
        </div>
      )}
      <div id={bodyId} className={collapsed ? "hidden" : "p-3"}>
        {children}
      </div>
    </section>
  );
}

/** Inline separator for query summaries: "MYCN · etp_status · 2 datasets". */
export function SummaryParts({ parts }: { parts: (ReactNode | null | false | undefined)[] }) {
  const shown = parts.filter((p) => p != null && p !== false && p !== "");
  return (
    <span className="flex flex-wrap items-baseline gap-x-2 gap-y-0.5">
      {shown.map((p, i) => (
        <span key={i} className="flex items-baseline gap-x-2">
          {i > 0 && <span aria-hidden="true" className="text-rule-firm">·</span>}
          <span>{p}</span>
        </span>
      ))}
    </span>
  );
}

/** A gene symbol, set the way gene symbols are set. */
export function Gene({ children }: { children: ReactNode }) {
  return <span className="font-mono italic text-ink">{children}</span>;
}

/**
 * How a pane input relates to the context bar: nothing when no context is
 * set, a quiet "from context" when it follows it, and a one-click way
 * back when the pane has overridden it.
 */
export function LinkHint({ link, display, noun = "gene" }: { link?: LinkState; display?: string; noun?: string }) {
  if (!link || link.contextValue == null) return null;
  const shown = display ?? link.contextValue;
  if (link.linked) {
    return (
      <span className="text-[11px] text-ink-mute" title={`Follows the context ${noun} set in the top bar`}>
        from context
      </span>
    );
  }
  return (
    <button
      type="button"
      onClick={link.relink}
      aria-label={`Use context ${noun} ${shown}`}
      className="max-w-[24ch] truncate text-[11px] font-medium text-accent underline decoration-accent/40 underline-offset-2 hover:decoration-accent"
    >
      use {shown}
    </button>
  );
}

/** Field label with an optional context-link hint on the right. */
export function FieldLabel({ children, link, display, noun }: { children: ReactNode; link?: LinkState; display?: string; noun?: string }) {
  if (!link || link.contextValue == null) return <label className="wb-label">{children}</label>;
  return (
    <div className="flex items-baseline justify-between gap-2">
      <label className="wb-label">{children}</label>
      <LinkHint link={link} display={display} noun={noun} />
    </div>
  );
}

export function GeneTagInput({
  genes,
  onChange,
  placeholder = "Type a gene symbol and press Enter",
}: {
  genes: string[];
  onChange: (genes: string[]) => void;
  placeholder?: string;
}) {
  // "CIP2A SET ANP32A" typed or pasted used to become ONE tag that no
  // dataset could resolve. Split with the gene panel's parser: spaces,
  // commas, newlines, "Heading:" prefixes and "(old name)" notes.
  const commit = (text: string) => {
    const add = parseGeneInput(text)
      .map((g) => g.gene.toUpperCase())
      .filter((g) => !genes.includes(g));
    if (add.length) onChange([...genes, ...add]);
  };
  return (
    <div>
      <div className="flex min-h-[30px] flex-wrap items-center gap-1 rounded-[3px] border border-field-edge bg-surface px-1.5 py-1 focus-within:border-accent">
        {genes.map((g) => (
          <span
            key={g}
            className="flex items-center gap-1 rounded-[2px] border border-rule bg-surface-sunk py-px pl-1.5 pr-0.5 font-mono text-[12px] italic text-ink"
          >
            {g}
            <button
              type="button"
              onClick={() => onChange(genes.filter((x) => x !== g))}
              className="flex h-4 w-4 items-center justify-center rounded-[2px] not-italic text-ink-mute hover:bg-surface-hover hover:text-ink"
              aria-label={`Remove ${g}`}
            >
              &times;
            </button>
          </span>
        ))}
        <input
          aria-label="Add a gene"
          type="text"
          placeholder={genes.length === 0 ? placeholder : "Add another"}
          className="min-w-[10ch] flex-1 bg-transparent px-1 py-0.5 font-mono text-[12.5px] italic text-ink outline-none placeholder:text-ink-mute placeholder:font-sans placeholder:not-italic"
          onKeyDown={(e) => {
            if ((e.key === "Enter" || e.key === ",") && e.currentTarget.value.trim()) {
              e.preventDefault();
              commit(e.currentTarget.value);
              e.currentTarget.value = "";
            } else if (e.key === "Backspace" && !e.currentTarget.value && genes.length) {
              onChange(genes.slice(0, -1));
            }
          }}
          // Text left in the box when the user clicks Score would
          // otherwise be silently dropped from the gene set.
          onBlur={(e) => {
            if (e.currentTarget.value.trim()) {
              commit(e.currentTarget.value);
              e.currentTarget.value = "";
            }
          }}
          onPaste={(e) => {
            const text = e.clipboardData.getData("text");
            if (!/[\s,;:()]/.test(text)) return; // a single symbol: let it paste normally
            e.preventDefault();
            commit(`${e.currentTarget.value} ${text}`);
            e.currentTarget.value = "";
          }}
        />
      </div>
    </div>
  );
}

export function PresetButton({
  label,
  active,
  onClick,
}: {
  label: string;
  active?: boolean;
  onClick: () => void;
}) {
  return (
    <button
      aria-pressed={active}
      type="button"
      onClick={onClick}
      className={[
        "h-[26px] whitespace-nowrap rounded-[3px] border px-2.5 font-mono text-[11.5px] transition-colors",
        active
          ? "border-accent bg-accent-soft font-medium text-accent-ink"
          : "border-rule-firm bg-surface text-ink-soft hover:border-ink-mute hover:text-ink",
      ].join(" ")}
    >
      {label}
    </button>
  );
}

export function PrimaryButton({
  children,
  onClick,
  disabled,
  loading,
}: {
  children: ReactNode;
  onClick: () => void;
  disabled?: boolean;
  loading?: boolean;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      disabled={disabled || loading}
      className="flex h-[30px] items-center justify-center gap-2 whitespace-nowrap rounded-[3px] bg-accent px-3.5 text-[13px] font-semibold text-white transition-colors hover:bg-accent-strong disabled:cursor-not-allowed disabled:bg-ink-mute disabled:opacity-60"
    >
      {loading && (
        <span className="h-3 w-3 animate-spin rounded-full border-2 border-white/40 border-t-white" />
      )}
      {children}
    </button>
  );
}

export function Stat({
  label,
  value,
  tone = "neutral",
}: {
  label: string;
  value: string;
  tone?: "neutral" | "hot" | "cool" | "accent";
}) {
  const toneClass = {
    neutral: "text-ink",
    hot: "text-hot",
    cool: "text-cool",
    accent: "text-accent",
  }[tone];
  return (
    <div className="flex flex-col gap-0.5">
      <span className="text-[11px] font-medium text-ink-mute">{label}</span>
      <span className={`font-mono text-[15px] font-semibold leading-tight tabular-nums ${toneClass}`}>{value}</span>
    </div>
  );
}

export function EmptyState({
  icon,
  title,
  description,
  hint,
}: {
  icon: ReactNode;
  title: string;
  description: string;
  hint?: ReactNode;
}) {
  return (
    <div className="flex flex-col items-center justify-center gap-1.5 rounded-[3px] border border-dashed border-rule-firm px-4 py-10 text-center">
      <div className="text-ink-mute">{icon}</div>
      <p className="text-[13.5px] font-medium text-ink-soft">{title}</p>
      <p className="max-w-[40ch] text-[12.5px] text-ink-mute">{description}</p>
      {hint && <p className="mt-2 max-w-[52ch] text-[12px] text-ink-mute">{hint}</p>}
    </div>
  );
}

export function ErrorState({ message }: { message: string }) {
  return (
    // role="alert" so a screen reader announces the failure.
    <div role="alert" className="rounded-[3px] border-l-[3px] border-hot bg-hot-soft px-4 py-3">
      <p className="text-[11.5px] font-semibold text-hot">Request failed</p>
      <p className="mt-1 text-[13px] text-ink">{message}</p>
    </div>
  );
}

export function DataTable({
  columns,
  rows,
}: {
  // srLabel names a column whose visible header is blank (a bar column).
  columns: { key: string; label: string; align?: "left" | "right"; srLabel?: string }[];
  rows: Record<string, ReactNode>[];
}) {
  return (
    // `sticky` pins relative to the nearest scrolling ancestor, not just
    // any ancestor with overflow set -- this div previously had
    // overflow-x-auto only, so on tall tables (TARGET's 469 rows) the
    // actual vertical scroll happened on a distant dockview pane several
    // levels up, and the header scrolled off-screen with the body instead
    // of staying pinned. max-h + overflow-y-auto makes this div itself the
    // scrolling ancestor sticky needs.
    // tabIndex + role: a scrolling region must be reachable by keyboard
    // (axe: scrollable-region-focusable), or a keyboard user cannot see
    // rows past the first screen.
    <div
      tabIndex={0}
      role="region"
      aria-label="Results table"
      className="max-h-[70vh] overflow-x-auto overflow-y-auto rounded-[3px] border border-rule"
    >
      <table className="w-full min-w-[28rem] border-collapse text-[12.5px]">
        <thead>
          <tr>
            {columns.map((c) => (
              <th
                key={c.key}
                className={`sticky top-0 z-10 border-b border-rule-firm bg-surface-sunk px-2.5 py-1.5 text-[11px] font-semibold text-ink-soft ${
                  c.align === "right" ? "text-right" : "text-left"
                }`}
              >
                {c.label || <span className="sr-only">{c.srLabel ?? c.key}</span>}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row, i) => (
            <tr key={i} className="border-b border-rule last:border-none hover:bg-ground">
              {columns.map((c) => (
                <td
                  key={c.key}
                  className={`px-2.5 py-1 tabular-nums ${c.align === "right" ? "text-right font-mono text-[12px]" : "text-left"}`}
                >
                  {row[c.key]}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
