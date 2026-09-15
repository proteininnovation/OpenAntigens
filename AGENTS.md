# Canonical repository

- This checkout belongs exclusively to `https://github.com/proteininnovation/OpenAntigens`.
- Before any fetch, push, pull request, merge, release, or deployment, verify that `origin` resolves to `proteininnovation/OpenAntigens` and that the target branch is in that repository.
- Never push, open pull requests, merge, release, or deploy from a personal or private fork. If the remote differs from the canonical repository, stop and correct it before doing any work.
- Changes intended for `main` must go through a pull request in `proteininnovation/OpenAntigens` unless the user explicitly requests a different workflow.
- Keep `core.hooksPath` set to `.githooks`; the tracked pre-push hook rejects every non-canonical destination.
