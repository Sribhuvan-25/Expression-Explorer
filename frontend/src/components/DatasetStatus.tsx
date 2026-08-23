import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";

export function DatasetStatus() {
  const { data, isLoading, isError } = useQuery({
    queryKey: ["dataset-info", "target_all_p2"],
    queryFn: () => api.datasetInfo("target_all_p2"),
    retry: false,
    staleTime: Infinity,
  });

  return (
    <div className="border-t border-rail-rule px-4 py-3">
      <p className="mb-1.5 font-mono text-[10px] uppercase tracking-wider text-rail-ink-mute">Active dataset</p>
      {isLoading && <p className="text-[12px] text-rail-ink-mute">Connecting…</p>}
      {isError && (
        <div className="flex items-start gap-1.5">
          <span className="mt-1 h-1.5 w-1.5 shrink-0 rounded-full bg-hot" />
          <p className="text-[12px] text-hot">API unreachable — is the backend running?</p>
        </div>
      )}
      {data && (
        <div className="flex items-start gap-1.5">
          <span className="mt-1 h-1.5 w-1.5 shrink-0 rounded-full bg-rail-accent" />
          <div>
            <p className="text-[12.5px] font-medium text-rail-ink">{data.display_name}</p>
            <p className="font-mono text-[10.5px] text-rail-ink-mute">{data.n_samples} samples</p>
          </div>
        </div>
      )}
    </div>
  );
}
