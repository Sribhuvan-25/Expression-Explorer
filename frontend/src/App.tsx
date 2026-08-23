import { BrowserRouter, Routes, Route, Navigate } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { DatasetStatus } from "./components/DatasetStatus";
import { ComparePage } from "./pages/ComparePage";
import { GenomeTracksPage } from "./pages/GenomeTracksPage";
import { RankPage } from "./pages/RankPage";
import { SignaturePage } from "./pages/SignaturePage";
import { SurvivalPage } from "./pages/SurvivalPage";

const queryClient = new QueryClient({
  defaultOptions: { queries: { refetchOnWindowFocus: false, retry: 1 } },
});

export default function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <BrowserRouter>
        <Routes>
          <Route element={<AppShell rail={<DatasetStatus />} />}>
            <Route index element={<Navigate to="/compare" replace />} />
            <Route path="/compare" element={<ComparePage />} />
            <Route path="/signature" element={<SignaturePage />} />
            <Route path="/rank" element={<RankPage />} />
            <Route path="/survival" element={<SurvivalPage />} />
            <Route path="/genome-tracks" element={<GenomeTracksPage />} />
          </Route>
        </Routes>
      </BrowserRouter>
    </QueryClientProvider>
  );
}
