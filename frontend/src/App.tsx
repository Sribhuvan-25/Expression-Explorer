import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppShell } from "./components/AppShell";
import { ApiError } from "./lib/api";

const queryClient = new QueryClient({
  defaultOptions: {
    queries: {
      refetchOnWindowFocus: false,
      // A retry only makes sense for a transient failure (a dropped
      // connection, a 5xx). ApiError carries the real HTTP status, so a
      // 4xx here means the request itself was invalid -- e.g. "gene not
      // found" -- and retrying sends the exact same request that will
      // fail the exact same way. Retrying those was also silently making
      // queries un-resolvable in some conditions: TanStack Query's retry
      // step gates on focusManager.isFocused(), which falls back to
      // document.visibilityState and only latches to a fresh reading on
      // its own visibilitychange listener -- confirmed via query-core's
      // source that once that's false, a query sits retrying forever with
      // no error and no data, and no page-visible symptom to explain why.
      // Not retrying 4xx sidesteps that failure mode entirely for the
      // whole class of error users actually hit (typos, unsupported
      // genes) instead of depending on a browser signal this app doesn't
      // otherwise care about.
      retry: (failureCount, error) => {
        if (error instanceof ApiError && error.status >= 400 && error.status < 500) return false;
        return failureCount < 1;
      },
    },
  },
});

export default function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <AppShell />
    </QueryClientProvider>
  );
}
