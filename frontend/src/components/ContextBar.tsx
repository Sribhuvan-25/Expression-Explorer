import { useQuery } from "@tanstack/react-query";
import { forwardRef, useEffect, useState } from "react";
import { api } from "../lib/api";
import { parseGeneInput } from "../lib/geneInput";
import { useWorkbench } from "../lib/workbench";

/**
 * The persistent context bar: a dataset and gene(s) set once and picked
 * up by every analysis pane (see lib/workbench.tsx for the rules). It is
 * a convenience for filling in inputs -- each pane still shows, and runs
 * on, its own inputs, so what ran is always stated on the pane itself.
 */
export const ContextBar = forwardRef<HTMLInputElement>(function ContextBar(_props, geneRef) {
  const ctx = useWorkbench();
  const { data } = useQuery({ queryKey: ["datasets"], queryFn: api.listDatasets });
  const datasets = data?.datasets ?? [];

  const applied = ctx?.genes.join(", ") ?? "";
  const [draft, setDraft] = useState(applied);
  useEffect(() => setDraft(applied), [applied]);

  if (!ctx) return null;

  const parsedDraft = parseGeneInput(draft).map((g) => g.gene.toUpperCase());
  const dirty = parsedDraft.join(", ") !== applied;
  const apply = () => ctx.setContext({ genes: parsedDraft });
  const clear = () => {
    setDraft("");
    ctx.setContext({ genes: [] });
  };

  return (
    <form
      role="search"
      aria-label="Analysis context"
      onSubmit={(e) => {
        e.preventDefault();
        apply();
      }}
      className="flex min-w-0 flex-1 flex-wrap items-center gap-x-3 gap-y-1.5 max-sm:basis-full"
    >
      {/* "optional" is part of the label: first-time users read an
          unlabelled bar at the top of every page as a required step. */}
      <span className="text-[11.5px] font-semibold text-ink-soft max-sm:hidden">
        Context <span className="font-normal text-ink-mute">(optional)</span>
      </span>

      <label className="flex min-w-0 items-center gap-1.5 max-sm:w-full">
        <span className="text-[11.5px] text-ink-mute max-sm:w-[52px] max-sm:shrink-0">Dataset</span>
        <select
          aria-label="Context dataset"
          aria-describedby={dirty ? undefined : "context-help"}
          value={ctx.datasetId}
          onChange={(e) => ctx.setContext({ datasetId: e.target.value })}
          className="wb-field wb-field-sm !w-auto min-w-0 max-w-[30ch] max-sm:!w-full max-sm:max-w-none"
        >
          <option value="">Each analysis picks its own</option>
          {datasets.map((d) => (
            <option key={d.dataset_id} value={d.dataset_id}>
              {d.display_name}
            </option>
          ))}
        </select>
      </label>

      <label className="flex min-w-0 items-center gap-1.5 max-sm:w-full">
        <span className="text-[11.5px] text-ink-mute max-sm:w-[52px] max-sm:shrink-0">Gene(s)</span>
        <span className="relative flex min-w-0 flex-1 items-center">
          <input
            ref={geneRef}
            aria-label="Context genes"
            aria-describedby={dirty ? undefined : "context-help"}
            value={draft}
            onChange={(e) => setDraft(e.target.value)}
            onKeyDown={(e) => {
              if (e.key === "Escape") {
                setDraft(applied);
                e.currentTarget.blur();
              }
            }}
            placeholder="Optional, e.g. MYCN"
            spellCheck={false}
            autoComplete="off"
            className="wb-field wb-field-sm wb-gene w-[26ch] min-w-0 pr-6 max-sm:w-full"
          />
          {applied && (
            <button
              type="button"
              onClick={clear}
              aria-label="Clear context genes"
              className="absolute right-1 flex h-5 w-5 items-center justify-center rounded-[2px] text-ink-mute hover:bg-surface-hover hover:text-ink"
            >
              <svg width="9" height="9" viewBox="0 0 10 10" aria-hidden="true">
                <path d="M1 1L9 9M9 1L1 9" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" />
              </svg>
            </button>
          )}
        </span>
      </label>

      {dirty ? (
        <button
          type="submit"
          className="h-[26px] rounded-[3px] bg-accent px-2.5 text-[12px] font-semibold text-white hover:bg-accent-strong"
        >
          {parsedDraft.length ? "Apply to all analyses" : "Clear"}
        </button>
      ) : (
        <span id="context-help" className="text-[11.5px] text-ink-mute sm:max-md:hidden">
          {applied
            ? "Shared by every analysis — a pane can override it."
            : "Not required — set a gene once and every analysis fills it in."}
        </span>
      )}
    </form>
  );
});
