import { useQuery } from "@tanstack/react-query";
import { api } from "../lib/api";

const EXTERNAL_LINKS: { label: string; href: (info: { symbol: string; ensembl_gene_id: string | null }) => string }[] = [
  { label: "NCBI", href: (i) => `https://www.ncbi.nlm.nih.gov/gene/?term=${encodeURIComponent(i.symbol)}` },
  { label: "Ensembl", href: (i) => `https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=${i.ensembl_gene_id}` },
  { label: "GeneCards", href: (i) => `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${encodeURIComponent(i.symbol)}` },
];

// Standalone, dataset-independent annotation: a gene symbol means the same
// thing everywhere, so this doesn't take a dataset id -- unlike every other
// panel in this app, which is scoped to one dataset's own matrix.
export function GeneAnnotation({ gene }: { gene: string }) {
  const { data, isLoading, isError } = useQuery({
    queryKey: ["gene-info", gene],
    queryFn: () => api.geneInfo(gene),
    retry: false,
    staleTime: Infinity, // annotation is stable; no reason to refetch within a session
  });

  if (isLoading) {
    return <p className="text-[11.5px] text-ink-mute">Looking up {gene}…</p>;
  }
  // A gene with no external annotation (e.g. a symbol only meaningful on
  // one platform, or a typo) isn't an error state for the surrounding
  // page -- the comparison itself doesn't depend on this panel, so fail
  // quiet rather than showing a red error box under a successful query.
  if (isError || !data) {
    return null;
  }

  return (
    <details className="group rounded-[3px] border border-rule bg-surface" open>
      <summary className="flex cursor-pointer list-none items-center gap-1.5 px-3 py-2 text-[12px] font-medium text-ink-soft hover:text-ink [&::-webkit-details-marker]:hidden">
        <svg width="9" height="9" viewBox="0 0 9 9" className="text-ink-mute transition-transform group-open:rotate-90" fill="none">
          <path d="M2.5 1L6.5 4.5L2.5 8" stroke="currentColor" strokeWidth="1.4" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
        {data.symbol}
        {data.name && <span className="truncate font-normal text-ink-mute">— {data.name}</span>}
      </summary>
      <div className="flex flex-col gap-2 border-t border-rule px-3 py-2.5">
        {data.summary && <p className="text-[12px] leading-relaxed text-ink-soft">{data.summary}</p>}
        {data.aliases.length > 0 && (
          <p className="text-[11px] text-ink-mute">
            <span className="font-mono uppercase tracking-wider text-ink-mute">Aliases</span>{" "}
            {data.aliases.join(", ")}
          </p>
        )}
        <div className="flex flex-wrap items-center gap-x-3 gap-y-1 font-mono text-[10.5px] text-ink-mute">
          {data.ensembl_gene_id && <span>{data.ensembl_gene_id}</span>}
          {data.entrez_gene_id != null && <span>Entrez {data.entrez_gene_id}</span>}
          {EXTERNAL_LINKS.filter((l) => l.label !== "Ensembl" || data.ensembl_gene_id).map((l) => (
            <a
              key={l.label}
              href={l.href({ symbol: data.symbol, ensembl_gene_id: data.ensembl_gene_id })}
              target="_blank"
              rel="noreferrer"
              className="text-accent hover:underline"
            >
              {l.label}
            </a>
          ))}
        </div>
      </div>
    </details>
  );
}
