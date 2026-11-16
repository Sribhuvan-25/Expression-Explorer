/**
 * A p-value as shown everywhere: scientific below 0.001, else 4 decimals.
 * Several stats printed every p in scientific notation, so p = 0.751
 * appeared as "7.51e-1".
 */
export function formatP(p: number | null | undefined): string {
  if (p == null || !Number.isFinite(p)) return "—";
  return p < 0.001 ? p.toExponential(2) : p.toFixed(4);
}
