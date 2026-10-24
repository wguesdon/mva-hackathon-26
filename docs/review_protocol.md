# Research review protocol

This document records the review method used during the project and its limits.
The documentation was added on October 24, 2026, after the owner's competition submissions.

## Research discussion and manuscript review

The research board and the adversarial manuscript review served different purposes.
The board proposed analysis plans through the restricted review adapters.
Separate model sessions reviewed manuscripts against a jury-style rubric and produced defect lists.
Historical manuscript review invocations included `codex --yolo` and separate Claude Code reviews in prepared workspaces.
Those historical invocations were not the restricted board route in `headless.py`.
This document does not prescribe unrestricted execution for a new review.

The reusable [rubric](../prompts/jury_review.md) preserves the review questions and output structure.
Its path instructions are adapted to this public repository, which has no submitted manuscript or completed research receipts.
The rubric is a prompt template. It is not a transcript or an official jury judgment.

## Review cycle

1. Define the claim, control and acceptance criterion before an analysis.
2. Record the relevant input hashes, code revision, parameters and tool versions.
3. Run the analysis and its declared controls.
4. Check each reported quantity against its owning output or receipt.
5. Give the reviewer the authorized manuscript, methods, code and evidence.
6. Require a defect list with file references and reasons.
7. Repair each defect or record why it remains unresolved.
8. Repeat affected checks after a repair.
9. Record the reviewed revision and the remaining limits.

A change to a decision gate after a result requires a dated amendment.
Keep the original criterion and state why it changed.
A missing record stays missing unless a documented reconstruction method supports it.
Label a reconstructed receipt and identify fields that the original run did not record.

## Distinct forms of evidence

| Check | What it establishes | What it does not establish |
|---|---|---|
| Synthetic software test | The code meets the tested expectation on the fixture. | It does not establish clinical accuracy. |
| Public benchmark | The method performs as measured on the declared benchmark. | It does not establish performance on every disease or variant class. |
| Reference resolution | A cited identifier resolves to a source. | It does not prove that the source supports the claim. |
| Manuscript quantity check | A reported value agrees with the checked record. | It does not establish the validity of the experiment or assumptions. |
| Model review | The reviewer identifies defects under the supplied rubric. | Agreement between models does not establish independent biological validation. |
| Human review | A person checked the declared evidence or decision. | A submission approval does not establish that every scientific judgment received human review. |

## Model attribution

Record the requested model separately from identifiers emitted by the service.
Record the CLI version, provider, product tier and applicable data-handling settings when known.
An initialization identifier does not prove which model produced every later response.
Do not reconstruct a model identity from a filename or an alias alone.
State the number of human-checked judgments where the analysis records that number.
State when an interactive session has no record that attributes an individual step.

## Limits of this export

The public repository includes the analysis code, synthetic tests, instructions and reusable templates.
It excludes completed conversations, manuscripts, private review outputs and historical receipts.
The command examples demonstrate setup and prompt construction without claiming to replay the original project.
The general registry measurement gate remains unimplemented, as documented in [the workflow](agent_workflow.md#decision-gates-and-evidence).
