# Canonical repository

- This checkout belongs exclusively to `https://github.com/proteininnovation/OpenAntigens`.
- Before any fetch, push, pull request, merge, release, or deployment, verify that `origin` resolves to `proteininnovation/OpenAntigens` and that the target branch is in that repository.
- Never push, open pull requests, merge, release, or deploy from a personal or private fork. If the remote differs from the canonical repository, stop and correct it before doing any work.
- Changes intended for `main` must go through a pull request in `proteininnovation/OpenAntigens` unless the user explicitly requests a different workflow.
- Never create or use a branch whose name begins with `codex`. Use an informative prefix such as `feat/`, `fix/`, `perf/`, `docs/`, `test/`, `refactor/`, or `chore/`.
- Keep `core.hooksPath` set to `.githooks`; the tracked pre-push hook rejects every non-canonical destination.

# Human and mouse portal parity

- Every change to shared page layout, navigation, copy, citations, agent guidance, documentation, or interactive behavior MUST apply to both human and mouse portals. A change is incomplete until both portals are checked. Do not update only the portal named in a request when the change also applies to the other.
- Differences MUST have a species-specific reason: target identity and provenance, available data and homolog mappings, species labels and theme, relative links, or supported features (for example, human disease context). Preserve those differences; do not copy human data or claims into mouse pages to make them look alike.
- Reuse shared renderers and content for common behavior. Inspect both build paths and their overrides before editing; remove stale overrides that prevent shared changes from reaching one portal.
- Verify the actual human and mouse build paths, not only the shared renderer or a custom preview. Refresh both local previews and inspect affected pages, copy, navigation, and links in each portal. Run relevant existing checks and add a regression check when it protects against divergence.
- Before reporting completion, publishing, or deploying, confirm that both portal artifacts contain the changes and state any intentional species-specific differences. A successful human build or preview alone does not establish mouse parity.
