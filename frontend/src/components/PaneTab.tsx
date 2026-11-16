import { useEffect, useRef, useState } from "react";
import type { IDockviewPanelHeaderProps } from "dockview-react";

/**
 * Pane tab without a button nested inside the focusable tab.
 *
 * dockview's default tab puts a close <button> inside its role="tab"
 * element; axe reports that as nested-interactive (assistive tech still
 * focuses the inner control, and a tab cannot contain one). Here the ×
 * is a pointer-only control hidden from assistive tech, and a keyboard
 * user closes the focused tab with Delete or Backspace -- stated in the
 * tab's accessible name.
 */
export function PaneTab({ api }: IDockviewPanelHeaderProps) {
  const ref = useRef<HTMLSpanElement>(null);
  const [title, setTitle] = useState(api.title ?? "");

  useEffect(() => {
    const sub = api.onDidTitleChange((e) => setTitle(e.title));
    return () => sub.dispose();
  }, [api]);

  useEffect(() => {
    const tab = ref.current?.closest<HTMLElement>(".dv-tab");
    if (!tab) return;
    tab.setAttribute("aria-label", `${title} — press Delete to close`);
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Delete" || e.key === "Backspace") {
        e.preventDefault();
        api.close();
      }
    };
    tab.addEventListener("keydown", onKey);
    return () => tab.removeEventListener("keydown", onKey);
  }, [api, title]);

  return (
    <span ref={ref} className="dv-default-tab">
      <span className="dv-default-tab-content">{title}</span>
      <span
        aria-hidden="true"
        className="dv-default-tab-action"
        title={`Close ${title}`}
        onPointerDown={(e) => e.stopPropagation()}
        onClick={(e) => {
          e.stopPropagation();
          api.close();
        }}
      >
        <svg height="11" width="11" viewBox="0 0 28 28" aria-hidden="true" focusable="false" className="dv-svg">
          <path d="M2.1 27.3L0 25.2L11.55 13.65L0 2.1L2.1 0L13.65 11.55L25.2 0L27.3 2.1L15.75 13.65L27.3 25.2L25.2 27.3L13.65 15.75L2.1 27.3Z" />
        </svg>
      </span>
    </span>
  );
}
