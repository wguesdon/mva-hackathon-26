# Local harness walkthrough

This example creates an invented board and previews its review prompts.
It uses the existing registry and board commands without a model call or cloud operation.
The example was added on October 24, 2026, after the owner's competition submissions.
It demonstrates the coordination mechanism and does not replay a historical research board.

## Preparation

Run these commands from the repository root in a fresh checkout.
Install the locked Python environment first.

```bash
uv sync --frozen
```

The example assignments use placeholder model identifiers.
The dry run validates their syntax but does not contact a provider or establish that a model exists.
Use the [assignment file](examples/reviewer_assignments.json) and [board goal](examples/board_goal.md) as written for this preview.

## Register the two tasks

Each task has its own question and acceptance criterion.
Both tasks share the same track, target, control and suite identity.
The identifiers describe an invented documentation exercise.

```bash
uv run python harness/research.py task \
  --task-id example_methods_review --phase 4 --track t2_drug \
  --assignee claude --data-class A --target synthetic_graph_review \
  --control example_baseline --suite-hash example_fixture_v1 \
  --question "Which inputs and outputs define the synthetic graph example?" \
  --criterion "Name source paths and state reproduction limits."
```

```bash
uv run python harness/research.py task \
  --task-id example_controls_review --phase 4 --track t2_drug \
  --assignee codex --data-class A --target synthetic_graph_review \
  --control example_baseline --suite-hash example_fixture_v1 \
  --question "Which checks separate software behavior from biological evidence?" \
  --criterion "Propose controls and identify unsupported conclusions."
```

The commands create local records in `tracking/research.json`.
Git ignores that directory.
A task contract is immutable. Reusing either task identifier will fail if the task already exists.

## Inspect the context and board

Preview the methods task context.

```bash
uv run python harness/research.py context --task-id example_methods_review --recipient claude
```

Create the board with two reviewers, three rounds and a six-call budget.
This command creates local board state but does not dispatch a reviewer.

```bash
uv run python harness/blackboard.py create \
  --board example_graph_board \
  --assignments docs/examples/reviewer_assignments.json \
  --goal-file docs/examples/board_goal.md --rounds 3 --max-calls 6
```

Preview the exact pending prompts.

```bash
uv run python harness/blackboard.py run --board example_graph_board --dry-run
```

Inspect the board status.

```bash
uv run python harness/blackboard.py status --board example_graph_board
```

The status remains `READY` with zero calls after the preview.
The local board views appear under `reports/blackboard/example_graph_board/`.
Git ignores these generated files.

## Optional live review

A live review requires owner authorization, configured provider credentials and installed CLIs.
Select approved model identifiers and create new tasks and a new board with the revised assignments.
The board freezes assignments at creation, so editing the example assignment file does not update an existing board.
Review the exact prompts and data class before dispatch.
Removing `--dry-run` from a board run dispatches model calls and can incur charges.

The configured backend restrictions are implemented in [`headless.py`](../harness/headless.py).
The public tests mock provider calls. Confirm the restrictions of an installed CLI before a live review.

## Other local checks

The synthetic graph demonstration uses invented identifiers.

```bash
make demo-track2
```

The repository checks cover software fixtures and tracked-file privacy rules.

```bash
make harness-tests
make privacy-scan-repo
```

The check results do not establish a treatment effect or reproduce the submitted scientific reports.
