/**
 * Expression units as the API reports them (`ExpressionUnit` in
 * backend/app/models/contract.py), formatted for a reader.
 */
export function formatUnit(unit: string | undefined | null): string {
  const map: Record<string, string> = {
    tpm: "TPM",
    fpkm: "FPKM",
    cpm: "CPM",
    counts: "counts",
    log2_tpm: "log2(TPM+1)",
    log2_intensity: "log2 intensity",
    delta_ct: "ΔCt",
  };
  if (!unit) return "expression";
  return map[unit] ?? unit.replace(/_/g, " ");
}

/**
 * True for units that are already on a log scale. Their zero is not "no
 * expression" (a log2 intensity of 0 is just a low reading), so a chart
 * must not anchor its axis at 0, and offering a further log transform
 * would plot log(log2(x)).
 */
export function isLogUnit(unit: string | undefined | null): boolean {
  return !!unit && (unit.startsWith("log2_") || unit === "delta_ct");
}
