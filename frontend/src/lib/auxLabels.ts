/**
 * What the CRISPR and drug results are called: title, y axis and PNG file
 * name. Each names both sides and the measurement, since the title read
 * "MYB vs MYB" (MYB expression against MYB's CRISPR gene effect) and an
 * exported figure cannot be read from that. Measurement names are the
 * layers' own (METHODS.md 8.3, 8.4): Chronos gene effect, and PRISM
 * viability log2 fold-change.
 */
export type AuxLayerKey = "crispr_gene_effect" | "drug_sensitivity";

export function auxResultTitle(layer: AuxLayerKey, gene: string, featureLabel: string): string {
  return layer === "crispr_gene_effect"
    ? `${gene} expression vs ${featureLabel} CRISPR gene effect`
    : `${gene} expression vs viability under ${featureLabel} (PRISM)`;
}

/** The y axis: which feature, and the measurement with its unit. */
export function auxAxisLabel(layer: AuxLayerKey, featureLabel: string): string {
  return layer === "crispr_gene_effect"
    ? `${featureLabel} CRISPR gene effect (Chronos)`
    : `Viability under ${featureLabel} (log2 fold-change)`;
}

/** File names keep only characters every file system accepts. */
export function auxExportFilename(args: {
  datasetId: string;
  layer: AuxLayerKey;
  gene: string;
  featureLabel: string;
  method: string;
}): string {
  const measure = args.layer === "crispr_gene_effect" ? "crispr-gene-effect" : "drug-viability";
  const safe = (s: string) => s.replace(/[^A-Za-z0-9._-]+/g, "_");
  return `${args.datasetId}-${safe(args.gene)}-expression-vs-${safe(args.featureLabel)}-${measure}-${args.method}`;
}
