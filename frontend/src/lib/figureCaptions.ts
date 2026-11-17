/**
 * Figure captions for the Correlation, Multi-omics and Signature Score
 * panes (METHODS.md 5.11). Each is one plain-text string, shown under the
 * figure with FigureCaption and handed unchanged to that figure's
 * ExportButton, so a saved PNG says what the on-screen figure says.
 *
 * A caption only describes what the backend computed (METHODS.md 2, 7.6,
 * 7.6a, 7.7, 8.2-8.5); it never adds a statistic.
 */
import { formatUnit } from "./units";

type CorrMethod = "pearson" | "spearman" | "kendall";

const METHOD_SENTENCE: Record<CorrMethod, string> = {
  pearson: "Pearson r is the linear correlation of the plotted values",
  spearman: "Spearman r is the rank correlation (rho), so the scale of the values does not affect it",
  kendall: "Kendall r is Kendall's tau, a rank concordance, so the scale of the values does not affect it",
};

const METHOD_NAME: Record<CorrMethod, string> = { pearson: "Pearson", spearman: "Spearman", kendall: "Kendall" };

/** Units the backend analyses as log2(x + 1) (`LINEAR_UNITS`, app/analysis/scale.py). */
const LINEAR_UNITS = new Set(["tpm", "fpkm", "cpm", "counts"]);

/** "log2(TPM+1)", plus how it was reached when the dataset stores linear values. */
function scalePhrase(shownUnit: string | undefined | null, datasetUnit: string | undefined | null): string {
  const shown = formatUnit(shownUnit ?? datasetUnit);
  if (datasetUnit && LINEAR_UNITS.has(datasetUnit) && shownUnit && shownUnit !== datasetUnit) {
    return `${shown} (the dataset's ${formatUnit(datasetUnit)} values, log2(x+1)-transformed before analysis)`;
  }
  return shown;
}

function listGenes(genes: string[], max = 12): string {
  if (genes.length <= max) return genes.join(", ");
  return `${genes.slice(0, max).join(", ")} and ${genes.length - max} more`;
}

const P_TWO_SIDED = "p is two-sided, against no correlation.";

export const POOLED_EXPERIMENTS_SENTENCE =
  "This dataset combines separate designed experiments (e.g. cell lines), so the result mostly reflects differences between them, not the perturbation.";

export function correlationCaption(args: {
  geneA: string;
  geneB: string;
  method: CorrMethod;
  n: number;
  valueUnit?: string | null;
  datasetUnit?: string | null;
  colorBy?: string;
  nColorGroups: number;
  /** The pooled-groups caveat, when the points span more than one colour group. */
  poolsNote?: string | null;
  pooledExperiments?: boolean;
}): string {
  const parts = [
    `Each point is one of ${args.n} samples: ${args.geneA} expression (x) against ${args.geneB} expression (y), on ${scalePhrase(args.valueUnit, args.datasetUnit)}.`,
    `${METHOD_SENTENCE[args.method]}, computed across all ${args.n} samples together; ${P_TWO_SIDED}`,
  ];
  if (args.colorBy && args.nColorGroups > 0) {
    parts.push(`Points are coloured by ${args.colorBy}; the colouring does not enter r.`);
  }
  if (args.poolsNote) parts.push(args.poolsNote);
  if (args.pooledExperiments) parts.push(POOLED_EXPERIMENTS_SENTENCE);
  return parts.join(" ");
}

export function networkCaption(args: {
  gene: string;
  method: "pearson" | "spearman";
  direction: "positive" | "negative";
  shown: number;
  n?: number | null;
  valueUnit?: string | null;
}): string {
  const across = args.n != null ? `across all ${args.n} samples` : "across the dataset's samples";
  const scale = args.valueUnit ? `, on ${formatUnit(args.valueUnit)}` : "";
  const end = args.direction === "positive" ? "highest (most positive)" : "most negative";
  const parts = [
    `Each row is one edge of ${args.gene}'s co-expression network: a gene whose expression correlates with ${args.gene} ${across}${scale}.`,
    `Every other gene that varies across these samples is scored by ${METHOD_NAME[args.method]} r with ${args.gene}, and the ${args.shown} with the ${end} r are listed in that order; there is no r threshold.`,
    "Ranked by r only, with no p-value or multiple-testing correction across the genes searched.",
  ];
  if (args.n != null && args.n < 20) parts.push(`With ${args.n} samples a high r is easily reached by chance.`);
  return parts.join(" ");
}

export function pcaCaption(args: {
  nSamples: number;
  genesUsed: string[];
  genesZeroVariance: string[];
  valueUnit?: string | null;
  datasetUnit?: string | null;
  colorBy?: string;
  pooledExperiments?: boolean;
}): string {
  const k = args.genesUsed.length;
  const parts = [
    `Each point is one of ${args.nSamples} samples, placed by its scores on the first two principal components of ${k} genes (${listGenes(args.genesUsed)}), on ${scalePhrase(args.valueUnit, args.datasetUnit)}.`,
    "Each gene was centred and scaled to unit variance first (PCA on the correlation matrix), so every gene weighs equally.",
    "Each axis gives that component's share of the total variance across these genes.",
  ];
  if (args.genesZeroVariance.length > 0) {
    parts.push(`Left out for having no variance: ${args.genesZeroVariance.join(", ")}.`);
  }
  if (args.colorBy) parts.push(`Coloured by ${args.colorBy}; the group labels were not used to compute the components.`);
  if (args.pooledExperiments) parts.push(POOLED_EXPERIMENTS_SENTENCE);
  return parts.join(" ");
}

export function mutationCaption(args: {
  gene: string;
  mutatedGene: string;
  expressionUnit?: string | null;
  nMutated: number;
  nWildtype: number;
  pValue: number | null;
  nDatasetTotal: number;
  nExcluded: number;
  exclusionReason: string | null;
}): string {
  const n = args.nMutated + args.nWildtype;
  const unit = args.expressionUnit ? ` (${formatUnit(args.expressionUnit)})` : "";
  const parts = [
    `${args.gene} expression${unit} in the ${n} samples with both expression and mutation data, split by whether the sample carries a non-silent ${args.mutatedGene} mutation (missense, nonsense, nonstop, frameshift or in-frame insertion/deletion, splice site, translation start): ${args.nMutated} mutated, ${args.nWildtype} wild-type.`,
    "Box: interquartile range; line: median; whiskers: minimum to maximum; dots: samples.",
    args.pValue != null
      ? "p: two-sided Mann-Whitney U test, which is rank-based, so the axis scale does not affect it."
      : "No test: a group has fewer than 2 samples.",
  ];
  if (args.nExcluded > 0) {
    parts.push(`${args.nExcluded} of ${args.nDatasetTotal} samples excluded (${args.exclusionReason ?? "no mutation data"}).`);
  }
  return parts.join(" ");
}

type AuxLayerKey = "crispr_gene_effect" | "drug_sensitivity";

export function auxCaption(args: {
  layer: AuxLayerKey;
  gene: string;
  featureLabel: string;
  featureId: string;
  moa?: string | null;
  method: "pearson" | "spearman";
  n: number;
  expressionUnit?: string | null;
  valueDescription: string;
  nDatasetTotal: number;
  nExcluded: number;
  exclusionReason: string | null;
}): string {
  const unit = args.expressionUnit ? ` (${formatUnit(args.expressionUnit)})` : "";
  const id = args.featureId !== args.featureLabel ? ` (${args.featureId})` : "";
  const y =
    args.layer === "crispr_gene_effect"
      ? `CRISPR knockout gene effect of ${args.featureLabel} (DepMap Chronos)`
      : `viability after exposure to ${args.featureLabel}${id} in the PRISM screen${args.moa ? `, mechanism ${args.moa}` : ""}`;
  const parts = [
    `Each point is one of ${args.n} samples: ${args.gene} expression${unit} (x) against ${y} (y).`,
    args.valueDescription,
    `${METHOD_SENTENCE[args.method]}, computed across the ${args.n} samples with both measurements; ${P_TWO_SIDED}`,
  ];
  if (args.nExcluded > 0) {
    parts.push(`${args.nExcluded} of ${args.nDatasetTotal} samples excluded: ${args.exclusionReason}.`);
  }
  return parts.join(" ");
}

/** Units the log2-mean score averages as they are (`LOG2_UNITS`, app/analysis/signature.py). */
const LOG2_SCORE_UNITS = new Set(["log2_tpm", "log2_intensity"]);

export function signatureCaption(args: {
  method: "auc" | "log2_mean";
  genes: string[];
  datasetUnit?: string | null;
  n: number;
  bins: number;
  shown: number;
}): string {
  const k = args.genes.length;
  const score =
    args.method === "auc"
      ? `AUCell scores of ${k} gene${k === 1 ? "" : "s"} (${listGenes(args.genes)}) for all ${args.n} samples. Within each sample genes are ranked by expression; the score is the area under the recovery curve of the signature genes within the top 25% of ranked genes, scaled so 1 means they are the sample's top-ranked genes and 0 means none reached the top 25%.`
      : `log2-mean scores of ${k} gene${k === 1 ? "" : "s"} (${listGenes(args.genes)}) for all ${args.n} samples: the mean across these genes of ${
          args.datasetUnit && LOG2_SCORE_UNITS.has(args.datasetUnit)
            ? `each sample's ${formatUnit(args.datasetUnit)} values, already on a log2 scale`
            : `log2(x+1) of each sample's ${formatUnit(args.datasetUnit)} values`
        }.`;
  return `Histogram of ${score} ${args.bins} equal-width bins from the lowest to the highest score; bar height is the number of samples. Shaded: the bins holding the top ${args.shown} samples listed below.`;
}
