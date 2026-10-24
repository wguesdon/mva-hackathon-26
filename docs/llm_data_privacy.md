# Privacy rules for the public harness

This document describes competition-approved model use, the project's data boundaries and the exported software controls.
It was added on October 24, 2026, after the owner's competition submissions.

## Competition-approved model use

The organizers' September 28 clarification permits consumer subscription accounts to process patient-specific or derived information under its stated conditions.
Those conditions require training opt-out, no feedback on outputs and no safety-review flag on the content.
The methods must record the provider, product tier and training setting.
The project's Claude Code sessions used Anthropic Max, and its Codex sessions used ChatGPT Pro.
Both used training opt-out with no feedback provided.
Use of these services under the organizer's conditions is permitted processing.

The source is the organizers' corrected guidance in [discussion 2](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2).

## Public repository boundary

Use public references, code and synthetic examples in this checkout.
Keep source genomic data, clinical documents, individual findings, completed prompts, transcripts and credentials outside the public repository.
The [repository scope](repository_scope.md) applies even when the competition permits publication of a limited finding elsewhere.
Do not copy an actual patient record into a synthetic example.

## Data classes

| Class | Meaning in the research project |
|---|---|
| A | The material contains public code, references, literature or synthetic examples. |
| B1 | A derived summary contains no individual genotype or phenotype. The project treated eligible B1 summaries as class A. |
| B2 | A limited derived finding contains individual-level information. It requires an approved processing route. |
| C | The material contains source reads, alignments, genome-scale variant or genotype records, or the full clinical document. |

The project prohibited class C content from every model route.
Approved B2 processing routes included Amazon Bedrock and subscription sessions that met the organizer's conditions above.
Kimi and OpenRouter routes remained at class A.

## Implemented dispatcher ceilings

The exported guard accepts class labels `A`, `B` and `C`. It does not implement separate `B1` and `B2` labels.
Its dispatch rules are more restrictive than the approved interactive-session policy.

| Backend | Model selection | Implemented ceiling |
|---|---|---|
| `claude` | The assignment supplies an explicit model. | A |
| `codex` | The assignment supplies an explicit model. | A |
| `opencode` | The model starts with `amazon-bedrock/`. | B |
| `opencode` | Any other model selection applies. | A |

The provider and plan strings in `arm_profile` describe intended configurations.
They do not inspect account settings, verify contractual retention or establish that an account has the declared subscription.
Do not change a task label to evade a class ceiling.
Use class A for every public walkthrough in this repository.

## Content and path checks

The scanner detects selected raw genomic record formats and local deny tokens.
The scanner also examines JSON-escaped tabs and newlines.
`.env.example` names the local deny-token setting. Keep its completed value in the ignored `.env` file.
Pattern matching cannot identify every sensitive sentence or prove that a document has public provenance.

The registry rejects source paths outside the project, symlinks, hidden paths and configured private directories.
The structured publication check validates allowed fields before pipeline outputs reach a less restricted surface.
The command guard checks selected restricted paths and shell forms.
These are controls at specific entry points. They do not replace operating-system isolation or authorize an arbitrary interactive command.

## Model calls and records

The board calls the swarm runner, which checks the outgoing prompt before dispatch.
The runner checks returned text before it registers a proposed finding.
The headless adapter removes unrelated environment credentials and redacts known credentials from completed transcripts.

The guard records its decision, the declared class and the prompt hash for each call through that entry point.
The `SENT` value records guard acceptance before the provider process starts. It does not prove that the provider received the prompt.
Call receipts and terminal events must establish whether execution completed.

## Software tests and retention

The synthetic tests check the implemented content, publication and command rules using synthetic inputs and public examples.
The [repository scope](repository_scope.md) defines which reusable files this export includes.

Keep new runtime material outside Git and apply the retention terms of its source.
The hackathon required deletion of source and genome-scale derived data within its stated retention period.
The owner completed a separate deletion process. This code export supplies no copy of that data.
