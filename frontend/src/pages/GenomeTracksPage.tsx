import { useEffect, useRef, useState } from "react";
import igv from "igv";
import { PageHeader, Panel, PresetButton } from "../components/ui";

// Public NCBI GEO FTP URLs for GSE225559 (Wang et al. 2025's own ChIP-seq
// and ATAC-seq deposit) -- bigWig signal tracks, RPM-normalized. No peak
// calls (BED/narrowPeak) were deposited, only continuous coverage, so this
// page shows raw signal for visual inspection, not called peaks.
const TRACKS = [
  {
    name: "ZMIZ1 ChIP-seq (LOUCY)",
    url: "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7050nnn/GSM7050319/suppl/GSM7050319_Loucy_Zmiz1.1m.bw",
    color: "#B4472E",
  },
  {
    name: "Input control (LOUCY)",
    url: "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7050nnn/GSM7050320/suppl/GSM7050320_Loucy_Input.1m.bw",
    color: "#767D82",
  },
  {
    name: "ATAC-seq (LOUCY)",
    url: "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7050nnn/GSM7050234/suppl/GSM7050234_LOUCY-ATAC_10000.1m.bw",
    color: "#2C6E9B",
  },
  {
    name: "ATAC-seq (THP-6)",
    url: "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7050nnn/GSM7050235/suppl/GSM7050235_THP6-ATAC_31000.1m.bw",
    color: "#1F6F5C",
  },
] as const;

// hg19 coordinates verified against UCSC's REST API (api.genome.ucsc.edu),
// not guessed. Padded beyond each gene body since the paper's putative
// enhancers sit at a distance from the gene itself (Fig. 3A-C: +140kb of
// MYB, +220kb of BCL2, +540kb of MYCN, -1kb of MEF2C) -- the padding here
// is generous enough to keep those enhancer sites in view, but exact
// enhancer coordinates were not re-derived from the paper's supplementary
// tables, so treat this as "the right neighborhood," not a precise locus.
const LOCI: Record<string, string> = {
  MYB: "chr6:135,400,000-135,750,000", // gene: 135,502,446-135,540,311; +140kb enhancer included
  MEF2C: "chr5:87,950,000-88,250,000", // gene: 88,012,934-88,200,074
  MYCN: "chr2:15,950,000-16,650,000", // gene: 16,080,672-16,087,129; +540kb enhancer included
  BCL2: "chr18:60,750,000-61,250,000", // gene: 60,790,579-60,987,361; +220kb enhancer included
  ZMIZ1: "chr10:80,750,000-81,150,000", // gene: 80,828,723-81,076,276
};

export function GenomeTracksPage() {
  const containerRef = useRef<HTMLDivElement>(null);
  const browserRef = useRef<any>(null);
  const [status, setStatus] = useState<"loading" | "ready" | "error">("loading");
  const [errorMsg, setErrorMsg] = useState("");
  const [locus, setLocus] = useState("MYB");

  useEffect(() => {
    if (!containerRef.current) return;
    let cancelled = false;

    igv
      .createBrowser(containerRef.current, {
        genome: "hg19",
        locus: LOCI[locus],
        tracks: TRACKS.map((t) => ({
          name: t.name,
          url: t.url,
          format: "bigwig",
          color: t.color,
          height: 60,
        })),
      })
      .then((browser: any) => {
        if (cancelled) return;
        // igv.js 3.8.5's createBrowser constructs browser.root correctly
        // but doesn't reliably append it to the given parent element in
        // this environment (confirmed: browser.parent === container, but
        // browser.root was never attached to it) -- attach it ourselves
        // as a defensive fallback so tracks actually render.
        if (browser.root && containerRef.current && !containerRef.current.contains(browser.root)) {
          containerRef.current.appendChild(browser.root);
        }
        browserRef.current = browser;
        setStatus("ready");
      })
      .catch((err: unknown) => {
        if (cancelled) return;
        setErrorMsg(err instanceof Error ? err.message : String(err));
        setStatus("error");
      });

    return () => {
      cancelled = true;
      if (browserRef.current) {
        igv.removeBrowser(browserRef.current);
        browserRef.current = null;
      }
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  useEffect(() => {
    if (browserRef.current && status === "ready") {
      browserRef.current.search(LOCI[locus]);
    }
  }, [locus, status]);

  return (
    <div className="mx-auto max-w-[1100px]">
      <PageHeader
        eyebrow="Genome Tracks"
        title="ZMIZ1 binding and chromatin accessibility"
        description="Raw ChIP-seq and ATAC-seq signal from Wang et al. 2025 (GSE225559), LOUCY and THP-6 cells, hg19. These are RPM-normalized coverage tracks — no peaks are called here, this is signal for visual inspection only."
      />

      <div className="flex flex-col gap-5 px-8 py-6">
        <Panel title="Locus">
          <div className="flex flex-wrap gap-1.5">
            {Object.keys(LOCI).map((gene) => (
              <PresetButton key={gene} label={gene} active={locus === gene} onClick={() => setLocus(gene)} />
            ))}
          </div>
          <p className="mt-2.5 text-[12px] text-ink-mute">
            Jumps the browser to the gene body plus flanking region (hg19). Pan and zoom freely once loaded — this
            is a live IGV instance, not a static image.
          </p>
        </Panel>

        <Panel>
          {status === "error" && (
            <div className="rounded-[3px] border-l-[3px] border-hot bg-hot-soft px-4 py-3">
              <p className="font-mono text-[10px] uppercase tracking-wider text-hot">Track load failed</p>
              <p className="mt-1 text-[13px] text-ink">{errorMsg}</p>
              <p className="mt-1 text-[12px] text-ink-mute">
                These files are hosted directly on NCBI's FTP server — this usually means a network or CORS issue
                reaching ftp.ncbi.nlm.nih.gov from the browser, not a bug in this page.
              </p>
            </div>
          )}
          {status === "loading" && <p className="py-8 text-center text-[13px] text-ink-mute">Loading tracks…</p>}
          <div ref={containerRef} className="igv-container" style={{ display: status === "error" ? "none" : "block" }} />
        </Panel>

        <Panel title="Tracks shown">
          <ul className="flex flex-col gap-1.5 text-[12.5px] text-ink-soft">
            {TRACKS.map((t) => (
              <li key={t.name} className="flex items-center gap-2">
                <span className="h-2.5 w-2.5 rounded-sm" style={{ background: t.color }} />
                {t.name}
              </li>
            ))}
          </ul>
        </Panel>
      </div>
    </div>
  );
}
