/** Gene-list text parsing shared by the gene panel and the tag input. */

export interface ParsedGene {
  gene: string;
  // Text the user wrote in parentheses after the symbol, e.g. the former
  // name "KIAA1524" -- shown beside the row, never sent to the API.
  note: string | null;
  group: string | null;
}

/**
 * Parses "Heading: A, B (OLDNAME), C" lines. A line with a colon starts a
 * category heading; a line without one continues the current heading (or
 * has none, for plain free text -- which then keeps input order). Only
 * the symbol before any parenthesis is sent, so the preset sends CIP2A
 * and PABIR1, not their legacy names. Exported for tests.
 */
export function parseGeneInput(text: string): ParsedGene[] {
  const out: ParsedGene[] = [];
  const seen = new Set<string>();
  let group: string | null = null;
  for (const line of text.split(/\n/)) {
    let rest = line;
    const m = line.match(/^\s*([^:,()]+):(.*)$/);
    if (m) {
      group = m[1].trim() || null;
      rest = m[2];
    }
    for (const t of rest.matchAll(/([A-Za-z0-9][A-Za-z0-9._-]*)(?:\s*\(([^)]*)\))?/g)) {
      const gene = t[1];
      const key = gene.toUpperCase();
      if (seen.has(key)) continue;
      seen.add(key);
      out.push({ gene, note: t[2]?.trim() || null, group });
    }
  }
  return out;
}
