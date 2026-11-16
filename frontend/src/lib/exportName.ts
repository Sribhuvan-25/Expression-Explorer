/**
 * A short, stable tag for a gene list in an export filename. Two different
 * panels from one dataset were both "<dataset>-contrast-panel.png", so the
 * second overwrote the first.
 */
export function geneListTag(genes: string[]): string {
  const key = genes.join(",");
  let h = 5381;
  for (let i = 0; i < key.length; i++) h = ((h * 33) ^ key.charCodeAt(i)) >>> 0;
  return `${genes.length}genes-${h.toString(36).slice(0, 6)}`;
}
