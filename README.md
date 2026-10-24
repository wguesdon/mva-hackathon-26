# MVA Hackathon 2026 code

This repository contains the agent software and analysis code for the MVA Hackathon 2026.
The code coordinates research tasks, ranks variant hypotheses and evaluates literature graphs for drug repurposing.
The tests use synthetic inputs and public reference examples.

## Agent setup and research workflow

The [workflow guide](docs/agent_workflow.md) describes the owner's setup role, agent roles, task contracts, discussion board and analysis records.
It distinguishes implemented controls from components that remained in the original design.
The [local walkthrough](docs/harness_walkthrough.md) demonstrates task registration and prompt previews without model calls.
The [review protocol](docs/review_protocol.md) explains decision gates, manuscript checks and the limits of model review.

[`AGENTS.md`](AGENTS.md) and [`CLAUDE.md`](CLAUDE.md) contain the same public operating rules.
The [privacy guide](docs/llm_data_privacy.md) explains the historical data classes and the exported guard behavior.
The [review rubric](prompts/jury_review.md) and [blank record templates](docs/templates/) preserve reusable parts of the research process.
These documents were added on October 24, 2026, after the owner's competition submissions.
They describe the workflow without distributing its private conversations or completed research records.

## Code map

| Component | Location | Function |
|---|---|---|
| Agent coordination | `harness/` | The task registry, discussion board and CLI adapters record assignments, responses and decisions. |
| Variant analysis | `pipelines/` | Nextflow workflows annotate variants, construct candidate pairs, examine structural variants and assess phase. |
| Submission CSV construction | `scripts/merge_track1_submission.py` | The script combines ranked small and structural variant rows within the row budget. |
| Public reference benchmarks | `scripts/track1_*.py` | The scripts evaluate variant calls and candidate ranks against caller-supplied reference inputs. |
| Literature graph analysis | `scripts/track2_run_disease.py` | The runner retrieves public relations, constructs a graph, ranks compounds and runs statistical controls. |
| Document format | `writeups/` | The Quarto extensions and blank example support report creation. |

The [method guide](docs/methods.md) explains the inputs, commands and outputs.
The [repository scope](docs/repository_scope.md) defines which files can be committed.

## Run the tests

Python 3.11 or later and `uv` are required.

```bash
uv sync --frozen
make harness-tests
make privacy-scan-repo
```

The tests do not call a model or submit an entry to the competition.
They create temporary inputs and check the code against expected outputs.
They do not establish clinical accuracy.

## Run a synthetic graph example

```bash
make demo-track2
```

The example creates invented genes, diseases, compounds and publication identifiers.
The ranker writes its results under `reports/`, which Git ignores.
No downloaded data is required for this example.

## Agent workflow

The software records each research task with its question, data class, control and acceptance criterion.
The discussion board sends bounded context to configured reviewers and records their responses.
The command adapters support Claude Code, Codex and OpenCode.
The adapters require separately installed CLIs and account configuration.
The tests mock these calls, so a test pass does not verify a live provider connection.

The code retains the scientific scoring logic of the research project.
This export removes machine-specific credential searches and private storage addresses.
Historical decision records and completed analysis outputs are not included.
The registry records proposals but does not implement the general measurement gate described in the original project design.
Analysis-specific programs and their records supported scientific claims outside that registry.
Some analysis commands require separately supplied records, tables or registered plans.
The [input requirements](docs/methods.md#inputs-not-distributed) distinguish those commands from the synthetic tests.

## Attribution

Copyright 2026 William Guesdon. The code and documentation use CC BY 4.0, as described in [LICENSE](LICENSE).
Downloaded resources retain their source licenses. [DATA_LICENSE](DATA_LICENSE) records the resource attribution.

> "This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership
> with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment
> Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to
> the child and their family who generously contributed their data and their story to advance
> research into this rare disease. We acknowledge their trust in making this Hackathon possible."

The hackathon dataset citation is Synapse `syn76251147`, DOI `10.7303/syn76251147`.
This repository does not distribute that dataset or individual findings.
