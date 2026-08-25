/**
 * CSV export for tabular results. Sample Ranking's whole purpose is
 * picking candidates from a long ranked list -- up to 469 rows for
 * TARGET-ALL-P2 -- and without an export the only way to act on that list
 * is manually retyping from an on-screen table.
 */

function csvCell(value: unknown): string {
  const s = value == null ? "" : String(value);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

export function exportRowsAsCsv(
  filename: string,
  columns: { key: string; label: string }[],
  rows: Record<string, unknown>[],
): void {
  const lines = [
    columns.map((c) => csvCell(c.label)).join(","),
    ...rows.map((row) => columns.map((c) => csvCell(row[c.key])).join(",")),
  ];
  const blob = new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename.endsWith(".csv") ? filename : `${filename}.csv`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}
