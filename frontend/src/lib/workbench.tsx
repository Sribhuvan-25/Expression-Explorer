/**
 * The workbench context: one dataset and gene list, set once in the
 * context bar and picked up by every analysis pane.
 *
 * The context only ever SEEDS a pane's own inputs. Each pane keeps its
 * own local state (so a pane can override -- correlate MYCN against MYC
 * while the context says MYCN), and the inputs are what get sent to the
 * API, exactly as before. Nothing here changes what an analysis computes.
 *
 * Rules:
 *  - Until the researcher sets a context, every pane keeps its own
 *    defaults (Survival's ETP-TF5 preset, CRISPR's MYB, ...).
 *  - When the context changes, every pane adopts it: single-gene inputs
 *    take the first gene, gene-set inputs take the whole list (a pane
 *    that needs at least N genes -- PCA -- only adopts a list that long),
 *    dataset pickers take the dataset if they offer it.
 *  - Editing a pane's input afterwards overrides it for that pane only,
 *    until the context changes again.
 *
 * Pages are rendered standalone in tests with no provider; every hook
 * here then behaves as plain local state with the page's own default.
 */
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from "react";

export interface WorkbenchState {
  genes: string[];
  datasetId: string;
}

interface WorkbenchValue extends WorkbenchState {
  setContext: (next: Partial<WorkbenchState>) => void;
}

const WorkbenchContext = createContext<WorkbenchValue | null>(null);

const STORAGE_KEY = "expression-explorer:workbench-context";

function readStored(): WorkbenchState {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    if (!raw) return { genes: [], datasetId: "" };
    const parsed = JSON.parse(raw) as Partial<WorkbenchState>;
    return {
      genes: Array.isArray(parsed.genes) ? parsed.genes.filter((g) => typeof g === "string") : [],
      datasetId: typeof parsed.datasetId === "string" ? parsed.datasetId : "",
    };
  } catch {
    return { genes: [], datasetId: "" };
  }
}

export function WorkbenchProvider({ children }: { children: ReactNode }) {
  const [state, setState] = useState<WorkbenchState>(readStored);

  const setContext = useCallback((next: Partial<WorkbenchState>) => {
    setState((prev) => {
      const merged = { ...prev, ...next };
      try {
        localStorage.setItem(STORAGE_KEY, JSON.stringify(merged));
      } catch {
        /* private window / blocked storage: the context still works for this session */
      }
      return merged;
    });
  }, []);

  const value = useMemo(() => ({ ...state, setContext }), [state, setContext]);
  return <WorkbenchContext.Provider value={value}>{children}</WorkbenchContext.Provider>;
}

export function useWorkbench(): WorkbenchValue | null {
  return useContext(WorkbenchContext);
}

/** What a pane input knows about its link to the context bar. */
export interface LinkState {
  /** The context's value for this input, or null when none is set. */
  contextValue: string | null;
  /** True when the input currently equals the context value. */
  linked: boolean;
  /** Puts the context value back into the input. */
  relink: () => void;
}

/**
 * Shared engine: local state seeded from, and re-seeded on every change
 * of, a context-derived value. `derive` maps the context to this input's
 * value, or null when this input should not adopt it.
 */
function useLinked<T>(
  fallback: T,
  derive: (ctx: WorkbenchState) => T | null,
  keyOf: (v: T) => string,
  onAdopt?: () => void,
): [T, (v: T) => void, T | null] {
  const ctx = useWorkbench();
  const derived = ctx ? derive(ctx) : null;
  const derivedKey = derived == null ? "" : keyOf(derived);

  const [value, setValue] = useState<T>(() => derived ?? fallback);
  // The context value this input last adopted. Seeded with the initial
  // adoption so StrictMode's double effect run doesn't re-adopt (and
  // re-fire onAdopt) on mount.
  const adopted = useRef(derivedKey);
  const onAdoptRef = useRef(onAdopt);
  onAdoptRef.current = onAdopt;

  useEffect(() => {
    if (!derivedKey || derivedKey === adopted.current || derived == null) return;
    adopted.current = derivedKey;
    setValue(derived);
    onAdoptRef.current?.();
    // eslint-disable-next-line react-hooks/exhaustive-deps -- keyed on derivedKey by design
  }, [derivedKey]);

  return [value, setValue, derived];
}

/** A single-gene input linked to the context's first gene. */
export function useLinkedGene(fallback: string): [string, (v: string) => void, LinkState] {
  const [value, setValue, derived] = useLinked<string>(
    fallback,
    (c) => c.genes[0] ?? null,
    (v) => v,
  );
  const link: LinkState = {
    contextValue: derived,
    linked: derived != null && value.trim().toUpperCase() === derived,
    relink: () => derived != null && setValue(derived),
  };
  return [value, setValue, link];
}

/**
 * A gene-set input linked to the context's gene list. `onAdopt` runs when
 * the context replaces the list (pages clear their "active preset" there,
 * since the list no longer is one).
 */
export function useLinkedGenes(
  fallback: string[],
  options: { minGenes?: number; onAdopt?: () => void } = {},
): [string[], (v: string[]) => void, LinkState] {
  const min = options.minGenes ?? 1;
  const [value, setValue, derived] = useLinked<string[]>(
    fallback,
    (c) => (c.genes.length >= min ? c.genes : null),
    (v) => v.join(","),
    options.onAdopt,
  );
  const link: LinkState = {
    contextValue: derived ? derived.join(", ") : null,
    linked: derived != null && value.join(",") === derived.join(","),
    relink: () => {
      if (derived) {
        setValue(derived);
        options.onAdopt?.();
      }
    },
  };
  return [value, setValue, link];
}

/** True when the context's gene list was adopted as this pane's initial value. */
export function useContextHasGenes(minGenes = 1): boolean {
  const ctx = useWorkbench();
  return (ctx?.genes.length ?? 0) >= minGenes;
}

/**
 * A dataset picker linked to the context dataset. `candidates` are the
 * datasets this pane can analyse; the context dataset is used only when
 * it is one of them, otherwise the pane's own first candidate (the
 * previous default) applies.
 */
export function useLinkedDataset(candidateIds: string[], defaultId?: string): [string, (v: string) => void, LinkState] {
  const ctx = useWorkbench();
  const ctxId = ctx?.datasetId ?? "";
  const ctxUsable = !!ctxId && candidateIds.includes(ctxId);
  const [datasetId, setDatasetId] = useState("");
  const adopted = useRef("");
  const candidatesKey = candidateIds.join(",");

  useEffect(() => {
    if (ctxUsable && ctxId !== adopted.current) {
      adopted.current = ctxId;
      setDatasetId(ctxId);
      return;
    }
    // Same default as before the context existed: the page's preferred
    // dataset, else the first candidate, whenever nothing valid is
    // selected yet.
    // Functional update: under StrictMode's double effect run, this
    // closure's datasetId can predate the context adoption just above,
    // and a plain set would overwrite it with the default.
    if (candidateIds.length > 0)
      setDatasetId((prev) =>
        candidateIds.includes(prev) ? prev : defaultId && candidateIds.includes(defaultId) ? defaultId : candidateIds[0],
      );
    // eslint-disable-next-line react-hooks/exhaustive-deps -- candidate list compared by key
  }, [ctxId, ctxUsable, candidatesKey, datasetId]);

  const link: LinkState = {
    contextValue: ctxUsable ? ctxId : null,
    linked: ctxUsable && datasetId === ctxId,
    relink: () => ctxUsable && setDatasetId(ctxId),
  };
  return [datasetId, setDatasetId, link];
}
