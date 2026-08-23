import type { IDockviewPanelProps } from "dockview-react";
import { ComparePage } from "../pages/ComparePage";
import { GenomeTracksPage } from "../pages/GenomeTracksPage";
import { RankPage } from "../pages/RankPage";
import { SignaturePage } from "../pages/SignaturePage";
import { SurvivalPage } from "../pages/SurvivalPage";

// Each existing page component becomes a pane's content, unmodified —
// panes are just a different frame around the same analysis logic. The
// wrapper below gives every pane its own scroll region and strips the
// page-level width centering (panes are already width-constrained by the
// dockview grid), without forking each page component.
//
// @container: page components wrap their query fields with `sm:flex-row`
// etc, which key off *viewport* width — fine for a single full-width
// page, but wrong once a page is squeezed into a ~300px pane sitting
// inside a wide window. Establishing a container context here lets the
// page components' field-wrap rules key off @container variants (see
// ui.tsx / page field rows) that respond to the pane's actual rendered
// width instead of the viewport.
function paneFrame(Component: React.ComponentType) {
  return function PaneContent(_props: IDockviewPanelProps) {
    return (
      <div className="@container h-full overflow-y-auto [&>div]:mx-0 [&>div]:max-w-none">
        <Component />
      </div>
    );
  };
}

export const PANE_COMPONENTS = {
  compare: paneFrame(ComparePage),
  signature: paneFrame(SignaturePage),
  rank: paneFrame(RankPage),
  survival: paneFrame(SurvivalPage),
  "genome-tracks": paneFrame(GenomeTracksPage),
};

export type PaneType = keyof typeof PANE_COMPONENTS;

export const PANE_LABELS: Record<PaneType, string> = {
  compare: "Expression Compare",
  signature: "Signature Score",
  rank: "Sample Ranking",
  survival: "Survival",
  "genome-tracks": "Genome Tracks",
};
