interface Source {
  name: string;
  kind: "paper" | "dataset";
  identifier: string;
  detail: string;
  url?: string;
}

// Static for now — becomes a real /sources endpoint once there are enough
// papers/datasets to justify the backend round trip (see redesign plan §7).
const SOURCES: Source[] = [
  {
    name: "Wang et al., J Exp Med 2025",
    kind: "paper",
    identifier: "10.1084/jem.20231349",
    detail: "seed paper",
    url: "https://doi.org/10.1084/jem.20231349",
  },
  {
    name: "TARGET-ALL-P2 (GDC)",
    kind: "dataset",
    identifier: "phs000218 / phs000463",
    detail: "469 samples",
  },
  {
    name: "DepMap 24Q4 Public",
    kind: "dataset",
    identifier: "via Figshare",
    detail: "186 lymphoid lines",
  },
  {
    name: "Liu et al., Nat Genet 2017",
    kind: "paper",
    identifier: "10.1038/ng.3909",
    detail: "ETP status classification",
    url: "https://doi.org/10.1038/ng.3909",
  },
];

export function Sources() {
  return (
    <details className="group border-t border-rule">
      <summary className="flex cursor-pointer list-none items-center gap-1.5 px-4 py-2.5 text-[12px] font-medium text-ink-soft hover:text-ink [&::-webkit-details-marker]:hidden">
        <svg
          width="9"
          height="9"
          viewBox="0 0 9 9"
          className="text-ink-mute transition-transform group-open:rotate-90"
          fill="none"
        >
          <path d="M2.5 1L6.5 4.5L2.5 8" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
        Sources
        <span className="font-mono text-[10px] text-ink-mute">({SOURCES.length})</span>
      </summary>
      <div className="flex flex-col gap-2.5 px-4 pb-3">
        {SOURCES.map((s) => (
          <div key={s.name} className="flex flex-col gap-0.5 border-t border-rule pt-2 first:border-t-0 first:pt-0">
            {s.url ? (
              <a href={s.url} target="_blank" rel="noreferrer" className="text-[12px] font-medium text-ink hover:text-accent">
                {s.name}
              </a>
            ) : (
              <span className="text-[12px] font-medium text-ink">{s.name}</span>
            )}
            <span className="font-mono text-[10px] text-ink-mute">
              {s.identifier} &middot; {s.detail}
            </span>
          </div>
        ))}
      </div>
    </details>
  );
}
