import { requestOpenPane } from "../lib/paneEvents";

interface TestSize {
  n_a: number;
  n_b: number;
  min_p_achievable?: number | null;
}

function joinAnd(parts: string[]): string {
  return parts.length <= 1 ? parts.join("") : `${parts.slice(0, -1).join(", ")} and ${parts[parts.length - 1]}`;
}

function formatMinP(p: number): string {
  return String(parseFloat(p.toFixed(3)));
}

/**
 * The tests whose group sizes make p < 0.05 impossible, one entry per
 * distinct size (6 pairwise tests on a 4 x 2-sample design are all
 * "2 vs 2" and should read as one statement, not six). Absent or null
 * min_p_achievable -- an older response -- is treated as "unknown", never
 * as "fine" or as an error.
 *
 * Exported for tests.
 */
export function underpoweredSizes(tests: TestSize[]): { label: string; minP: number }[] {
  const seen = new Map<string, number>();
  for (const t of tests) {
    const m = t.min_p_achievable;
    if (typeof m !== "number" || !Number.isFinite(m) || m <= 0.05) continue;
    // Unordered: "4 vs 3" and "3 vs 4" have the same floor and are one
    // statement (DepMap by subtype listed both).
    const [hi, lo] = t.n_a >= t.n_b ? [t.n_a, t.n_b] : [t.n_b, t.n_a];
    const label = `${hi} vs ${lo}`;
    if (!seen.has(label)) seen.set(label, m);
  }
  return Array.from(seen, ([label, minP]) => ({ label, minP })).sort((a, b) => b.minP - a.minP);
}

/**
 * An exact Mann-Whitney on 2 vs 2 samples cannot produce p below 0.333,
 * so a non-significant p there is uninformative, not evidence of no
 * difference. Without this, "p = 0.33" on a 2-vs-2 design reads as a
 * negative result.
 */
export function MinPNotice({
  tests,
  hasContrasts = false,
  partial = false,
  smallestShownP,
}: {
  tests: TestSize[];
  hasContrasts?: boolean;
  // True when only some of the tests on screen are affected, so the
  // notice does not imply every comparison is underpowered.
  partial?: boolean;
  // Smallest p actually displayed next to this notice. If it is below the
  // exact-test floor, those p-values cannot have come from an exact test
  // (the genome-wide ranking uses an approximation), and the notice must
  // say so rather than contradict the table beneath it.
  smallestShownP?: number | null;
}) {
  const sizes = underpoweredSizes(tests);
  if (sizes.length === 0) return null;

  const floor = Math.min(...sizes.map((s) => s.minP));
  const belowFloor =
    typeof smallestShownP === "number" && Number.isFinite(smallestShownP) && smallestShownP < floor - 1e-9;

  const statement =
    sizes.length === 1
      ? `With n = ${sizes[0].label}, the smallest p-value this test can produce is ${formatMinP(sizes[0].minP)}`
      : `With these group sizes, the smallest p-value this test can produce is ${joinAnd(
          sizes.map((s) => `${formatMinP(s.minP)} for n = ${s.label}`),
        )}`;

  return (
    <div
      role="note"
      className="mb-3 rounded-[3px] border-l-[3px] border-warn bg-warn-soft px-3 py-2 text-[12.5px] text-ink"
    >
      <p>
        {statement} — p ≥ 0.05 here does not mean 'no difference'.
        {partial && " Other comparisons shown have enough samples to reach p < 0.05."}
      </p>
      {belowFloor && (
        <p data-testid="below-floor" className="mt-1.5 text-[12px] text-ink-soft">
          Some p-values below are smaller than that limit, so they come from a normal approximation rather than the
          exact test; at these group sizes the approximation is unreliable — treat them with caution.
        </p>
      )}
      {hasContrasts && (
        <p className="mt-1.5 text-[12px] text-ink-soft">
          This dataset has precomputed count-based contrasts (moderated t, which borrows strength across
          genes) for its designed experiments.{" "}
          <button
            type="button"
            onClick={() => requestOpenPane("contrasts")}
            className="font-medium text-accent underline decoration-accent/40 underline-offset-2 hover:decoration-accent"
          >
            Open Perturbation Contrasts
          </button>
        </p>
      )}
    </div>
  );
}
