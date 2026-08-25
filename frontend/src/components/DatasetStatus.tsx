import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";

// This is a connectivity/registry check, not a per-pane indicator -- it
// used to be labeled "Active dataset" and hardcoded to one dataset id,
// which read as if it tracked whatever the focused pane had selected. It
// never did: each pane picks its own dataset independently, and multiple
// panes open at once can each be on a different one. Listing every
// registered dataset here is what the sidebar can actually vouch for.
export function DatasetStatus() {
  const { data, isLoading, isError } = useQuery({
    queryKey: ["datasets"],
    queryFn: api.listDatasets,
    retry: false,
    staleTime: Infinity,
  });

  return (
    <div className="border-t border-rail-rule px-4 py-3">
      <p className="mb-1.5 font-mono text-[10px] uppercase tracking-wider text-rail-ink-mute">Datasets available</p>
      {isLoading && <p className="text-[12px] text-rail-ink-mute">Connecting…</p>}
      {isError && (
        <div className="flex items-start gap-1.5">
          <span className="mt-1 h-1.5 w-1.5 shrink-0 rounded-full bg-hot" />
          <p className="text-[12px] text-hot">API unreachable — is the backend running?</p>
        </div>
      )}
      {data && (
        <div className="flex flex-col gap-1.5">
          {data.datasets.map((d) => (
            <div key={d.dataset_id} className="flex items-start gap-1.5">
              <span className="mt-1 h-1.5 w-1.5 shrink-0 rounded-full bg-rail-accent" />
              <p className="text-[12.5px] font-medium text-rail-ink">{d.display_name}</p>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
