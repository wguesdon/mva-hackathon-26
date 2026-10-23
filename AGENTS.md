# Repository rules

This repository contains reusable research code and synthetic tests.
Keep patient records, findings, submission files, prompts from real runs and credentials outside this repository.
Do not copy history from another repository.

Use `uv` for Python and Podman for containers.
Keep generated outputs in the ignored output directories.
Read `docs/repository_scope.md` before you add a file.
Run `make harness-tests` and `make privacy-scan-repo` before a push.
Use Conventional Commits and do not add an assistant attribution to a commit.
