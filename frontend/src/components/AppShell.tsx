import { NavLink, Outlet } from "react-router-dom";
import type { ReactNode } from "react";

const NAV_ITEMS = [
  { to: "/compare", label: "Expression Compare", hint: "Gene across groups" },
  { to: "/signature", label: "Signature Score", hint: "Gene set, ranked" },
  { to: "/rank", label: "Sample Ranking", hint: "Pick candidates" },
  { to: "/survival", label: "Survival", hint: "Kaplan–Meier + Cox" },
] as const;

function MarkGlyph() {
  return (
    <svg width="22" height="22" viewBox="0 0 22 22" fill="none" aria-hidden="true">
      <path
        d="M4 18C4 18 6 4 11 4C16 4 18 18 18 18"
        stroke="var(--accent)"
        strokeWidth="1.6"
        strokeLinecap="round"
      />
      <path
        d="M4 11C4 11 7 15 11 11C15 7 18 11 18 11"
        stroke="var(--accent)"
        strokeWidth="1.6"
        strokeLinecap="round"
        opacity="0.55"
      />
    </svg>
  );
}

function NavRow({ to, label, hint }: { to: string; label: string; hint: string }) {
  return (
    <NavLink
      to={to}
      className={({ isActive }) =>
        [
          "group flex flex-col gap-0.5 rounded-md px-3 py-2.5 transition-colors",
          isActive
            ? "bg-accent-soft text-accent-ink"
            : "text-ink-soft hover:bg-surface-hover hover:text-ink",
        ].join(" ")
      }
    >
      {({ isActive }) => (
        <>
          <span className="flex items-center gap-2 text-[13.5px] font-medium">
            <span
              className={[
                "h-1.5 w-1.5 rounded-full transition-colors",
                isActive ? "bg-accent" : "bg-rule-firm group-hover:bg-ink-mute",
              ].join(" ")}
            />
            {label}
          </span>
          <span className="pl-3.5 font-mono text-[10.5px] uppercase tracking-wider text-ink-mute">
            {hint}
          </span>
        </>
      )}
    </NavLink>
  );
}

export function AppShell({ rail }: { rail?: ReactNode }) {
  return (
    <div className="flex h-screen overflow-hidden bg-ground text-ink">
      <aside className="flex w-[248px] shrink-0 flex-col border-r border-rule bg-surface">
        <div className="flex items-center gap-2.5 border-b border-rule px-4 py-4">
          <MarkGlyph />
          <div className="flex flex-col leading-tight">
            <span className="font-display text-[15px] font-semibold tracking-tight">
              Expression Explorer
            </span>
            <span className="font-mono text-[10px] uppercase tracking-wider text-ink-mute">
              ETP-ALL &middot; v0
            </span>
          </div>
        </div>

        <nav className="flex flex-col gap-1 px-2.5 py-3">
          <span className="px-3 pb-1 pt-1 font-mono text-[10px] uppercase tracking-wider text-ink-mute">
            Analyses
          </span>
          {NAV_ITEMS.map((item) => (
            <NavRow key={item.to} {...item} />
          ))}
        </nav>

        <div className="flex-1" />

        {rail}

        <div className="border-t border-rule px-4 py-3">
          <p className="font-mono text-[10px] leading-relaxed text-ink-mute">
            Data: TARGET-ALL-P2 (GDC), DepMap 24Q4
            <br />
            Method: Wang et al., J Exp Med 2025
          </p>
        </div>
      </aside>

      <main className="flex-1 overflow-y-auto">
        <Outlet />
      </main>
    </div>
  );
}
