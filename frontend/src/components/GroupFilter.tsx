import { useMemo, useState } from "react";

interface GroupCount {
  group: string;
  count: number;
}

export function useGroupFilter(points: { group: string }[]) {
  const groupCounts = useMemo<GroupCount[]>(() => {
    const counts = new Map<string, number>();
    for (const p of points) counts.set(p.group, (counts.get(p.group) ?? 0) + 1);
    return Array.from(counts.entries())
      .map(([group, count]) => ({ group, count }))
      .sort((a, b) => b.count - a.count);
  }, [points]);

  const [selected, setSelected] = useState<Set<string> | null>(null);

  // Reset to "all selected" whenever the underlying group set changes
  // (new query), rather than carrying a stale selection across queries.
  const groupKey = groupCounts.map((g) => g.group).join("|");
  const [lastKey, setLastKey] = useState(groupKey);
  if (groupKey !== lastKey) {
    setLastKey(groupKey);
    setSelected(null);
  }

  const effectiveSelected = selected ?? new Set(groupCounts.map((g) => g.group));

  return { groupCounts, selected: effectiveSelected, setSelected };
}

export function GroupFilter({
  groupCounts,
  selected,
  onChange,
}: {
  groupCounts: { group: string; count: number }[];
  selected: Set<string>;
  onChange: (next: Set<string>) => void;
}) {
  if (groupCounts.length <= 6) return null; // not worth a filter UI for a handful of groups

  const toggle = (group: string) => {
    const next = new Set(selected);
    if (next.has(group)) next.delete(group);
    else next.add(group);
    onChange(next);
  };

  const allSelected = selected.size === groupCounts.length;

  return (
    <details className="rounded-[3px] border border-rule bg-surface">
      <summary className="flex cursor-pointer list-none items-center justify-between px-3 py-2 text-[12px] font-medium text-ink-soft hover:text-ink [&::-webkit-details-marker]:hidden">
        <span className="flex items-center gap-1.5">
          <svg width="9" height="9" viewBox="0 0 9 9" className="text-ink-mute" fill="none">
            <path d="M2.5 1L6.5 4.5L2.5 8" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" strokeLinejoin="round" />
          </svg>
          Filter groups shown in chart
        </span>
        <span className="font-mono text-[10.5px] text-ink-mute">
          {selected.size} / {groupCounts.length}
        </span>
      </summary>
      <div className="max-h-[220px] overflow-y-auto border-t border-rule px-3 py-2">
        <div className="mb-2 flex gap-3">
          <button
            type="button"
            onClick={() => onChange(new Set(groupCounts.map((g) => g.group)))}
            disabled={allSelected}
            className="font-mono text-[10.5px] text-accent hover:underline disabled:text-ink-mute disabled:no-underline"
          >
            Select all
          </button>
          <button
            type="button"
            onClick={() => onChange(new Set())}
            disabled={selected.size === 0}
            className="font-mono text-[10.5px] text-accent hover:underline disabled:text-ink-mute disabled:no-underline"
          >
            Clear
          </button>
          <button
            type="button"
            onClick={() => onChange(new Set(groupCounts.slice(0, 10).map((g) => g.group)))}
            className="font-mono text-[10.5px] text-accent hover:underline"
          >
            Top 10 by n
          </button>
        </div>
        <div className="flex flex-col gap-0.5">
          {groupCounts.map((g) => (
            <label
              key={g.group}
              className="flex cursor-pointer items-center gap-2 rounded-[3px] px-1.5 py-1 text-[12px] hover:bg-surface-hover"
            >
              <input
                type="checkbox"
                checked={selected.has(g.group)}
                onChange={() => toggle(g.group)}
                className="h-3.5 w-3.5 shrink-0 accent-accent"
              />
              <span className="min-w-0 flex-1 truncate text-ink-soft" title={g.group}>
                {g.group}
              </span>
              <span className="shrink-0 font-mono text-[10.5px] text-ink-mute">n={g.count}</span>
            </label>
          ))}
        </div>
      </div>
    </details>
  );
}
