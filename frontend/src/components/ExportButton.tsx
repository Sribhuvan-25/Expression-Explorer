import { useState } from "react";
import { exportSvgAsPng } from "../lib/exportFigure";

export function ExportButton({
  svgRef,
  filename,
  title,
  subtitle,
  statLines,
}: {
  svgRef: React.RefObject<SVGSVGElement | null>;
  filename: string;
  title: string;
  subtitle?: string;
  statLines: string[];
}) {
  const [state, setState] = useState<"idle" | "exporting" | "error">("idle");

  const handleClick = async () => {
    if (!svgRef.current) return;
    setState("exporting");
    try {
      await exportSvgAsPng(svgRef.current, filename, { title, subtitle, statLines });
      setState("idle");
    } catch {
      setState("error");
      setTimeout(() => setState("idle"), 2000);
    }
  };

  return (
    <button
      type="button"
      onClick={handleClick}
      disabled={state === "exporting"}
      title="Download as PNG"
      className="flex items-center gap-1 rounded-[3px] border border-rule-firm bg-surface px-2 py-1 font-mono text-[10.5px] text-ink-mute transition-colors hover:border-accent hover:text-accent disabled:opacity-50"
    >
      <svg width="10" height="10" viewBox="0 0 10 10" fill="none">
        <path d="M5 1v6M2.5 5L5 7.5L7.5 5M1.5 8.5h7" stroke="currentColor" strokeWidth="1.2" strokeLinecap="round" strokeLinejoin="round" />
      </svg>
      {state === "error" ? "Failed" : state === "exporting" ? "Exporting…" : "PNG"}
    </button>
  );
}
