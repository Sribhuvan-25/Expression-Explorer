import type { ReactNode } from "react";

export function PageHeader({
  eyebrow,
  title,
  description,
}: {
  eyebrow: string;
  title: string;
  description: string;
}) {
  return (
    <header className="border-b border-rule px-8 py-6">
      <p className="font-mono text-[11px] uppercase tracking-widest text-accent">{eyebrow}</p>
      <h1 className="mt-1 font-display text-[26px] font-semibold tracking-tight text-ink">
        {title}
      </h1>
      <p className="mt-1.5 max-w-[60ch] text-[13.5px] text-ink-soft">{description}</p>
    </header>
  );
}

export function Panel({
  title,
  action,
  children,
  className = "",
}: {
  title?: string;
  action?: ReactNode;
  children: ReactNode;
  className?: string;
}) {
  return (
    <section className={`rounded-[4px] border border-rule bg-surface ${className}`}>
      {title && (
        <div className="flex items-center justify-between border-b border-rule px-4 py-2.5">
          <h2 className="text-[12.5px] font-semibold text-ink">{title}</h2>
          {action}
        </div>
      )}
      <div className="p-4">{children}</div>
    </section>
  );
}

export function GeneTagInput({
  genes,
  onChange,
  placeholder = "Type a gene symbol and press Enter",
}: {
  genes: string[];
  onChange: (genes: string[]) => void;
  placeholder?: string;
}) {
  return (
    <div>
      <div className="flex flex-wrap gap-1.5 rounded-[3px] border border-rule bg-ground p-2 focus-within:border-accent">
        {genes.map((g) => (
          <span
            key={g}
            className="flex items-center gap-1 rounded bg-accent-soft px-2 py-0.5 font-mono text-[12.5px] font-medium italic text-accent-ink"
          >
            {g}
            <button
              type="button"
              onClick={() => onChange(genes.filter((x) => x !== g))}
              className="text-accent-ink/60 hover:text-accent-ink"
              aria-label={`Remove ${g}`}
            >
              &times;
            </button>
          </span>
        ))}
        <input
          type="text"
          placeholder={genes.length === 0 ? placeholder : "Add another"}
          className="min-w-[10ch] flex-1 bg-transparent px-1 py-0.5 font-mono text-[12.5px] text-ink outline-none placeholder:text-ink-mute placeholder:font-sans placeholder:not-italic"
          onKeyDown={(e) => {
            const val = e.currentTarget.value.trim().toUpperCase();
            if ((e.key === "Enter" || e.key === ",") && val) {
              e.preventDefault();
              if (!genes.includes(val)) onChange([...genes, val]);
              e.currentTarget.value = "";
            } else if (e.key === "Backspace" && !e.currentTarget.value && genes.length) {
              onChange(genes.slice(0, -1));
            }
          }}
        />
      </div>
    </div>
  );
}

export function PresetButton({
  label,
  active,
  onClick,
}: {
  label: string;
  active?: boolean;
  onClick: () => void;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      className={[
        "whitespace-nowrap rounded-[3px] border px-3 py-1 font-mono text-[11px] transition-colors",
        active
          ? "border-accent bg-accent-soft text-accent-ink"
          : "border-rule text-ink-mute hover:border-rule-firm hover:text-ink-soft",
      ].join(" ")}
    >
      {label}
    </button>
  );
}

export function PrimaryButton({
  children,
  onClick,
  disabled,
  loading,
}: {
  children: ReactNode;
  onClick: () => void;
  disabled?: boolean;
  loading?: boolean;
}) {
  return (
    <button
      type="button"
      onClick={onClick}
      disabled={disabled || loading}
      className="flex items-center justify-center gap-2 whitespace-nowrap rounded-[3px] bg-accent px-4 py-2 text-[13px] font-semibold text-white transition-colors hover:bg-accent-strong disabled:cursor-not-allowed disabled:opacity-50"
    >
      {loading && (
        <span className="h-3 w-3 animate-spin rounded-full border-2 border-white/40 border-t-white" />
      )}
      {children}
    </button>
  );
}

export function Stat({
  label,
  value,
  tone = "neutral",
}: {
  label: string;
  value: string;
  tone?: "neutral" | "hot" | "cool" | "accent";
}) {
  const toneClass = {
    neutral: "text-ink",
    hot: "text-hot",
    cool: "text-cool",
    accent: "text-accent",
  }[tone];
  return (
    <div className="flex flex-col gap-0.5">
      <span className="font-mono text-[10px] uppercase tracking-wider text-ink-mute">{label}</span>
      <span className={`font-mono text-[15px] font-semibold tabular-nums ${toneClass}`}>{value}</span>
    </div>
  );
}

export function EmptyState({ icon, title, description }: { icon: ReactNode; title: string; description: string }) {
  return (
    <div className="flex flex-col items-center justify-center gap-2 rounded-[4px] border border-dashed border-rule-firm py-16 text-center">
      <div className="text-ink-mute">{icon}</div>
      <p className="text-[13.5px] font-medium text-ink-soft">{title}</p>
      <p className="max-w-[40ch] text-[12.5px] text-ink-mute">{description}</p>
    </div>
  );
}

export function ErrorState({ message }: { message: string }) {
  return (
    <div className="rounded-[3px] border-l-[3px] border-hot bg-hot-soft px-4 py-3">
      <p className="font-mono text-[10px] uppercase tracking-wider text-hot">Request failed</p>
      <p className="mt-1 text-[13px] text-ink">{message}</p>
    </div>
  );
}

export function DataTable({
  columns,
  rows,
}: {
  columns: { key: string; label: string; align?: "left" | "right" }[];
  rows: Record<string, ReactNode>[];
}) {
  return (
    <div className="overflow-x-auto rounded-[3px] border border-rule">
      <table className="w-full min-w-[28rem] border-collapse text-[12.5px]">
        <thead>
          <tr>
            {columns.map((c) => (
              <th
                key={c.key}
                className={`sticky top-0 border-b border-rule bg-surface-sunk px-3 py-1.5 font-mono text-[10px] uppercase tracking-wider text-ink-mute ${
                  c.align === "right" ? "text-right" : "text-left"
                }`}
              >
                {c.label}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row, i) => (
            <tr key={i} className="border-b border-rule last:border-none">
              {columns.map((c) => (
                <td
                  key={c.key}
                  className={`px-3 py-1.5 tabular-nums ${c.align === "right" ? "text-right" : "text-left"}`}
                >
                  {row[c.key]}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
