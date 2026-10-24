# Repository scope

This repository contains reusable source code, public resource configuration and synthetic test fixtures.
It starts with new Git history. It does not inherit the research repository's commits or tags.

## Included files

The export includes agent coordination, variant analysis, literature retrieval, graph analysis, CSV construction and blank report formatting.
Public gene panels describe literature and ontology evidence. They do not describe an individual sample.
Synthetic tests define their own expected variants and graph relationships.
Generic parameter names such as `patient_vcf` identify inputs that a caller must supply.
The repository also includes public operating instructions, a workflow description, reusable review prompts and blank decision and receipt templates.
The local walkthrough uses invented tasks and public documentation.

## Excluded files

The export excludes individual findings, clinical descriptions, submitted predictions, completed methods forms, papers, slides and video materials.
It also excludes model conversations, review outputs, session summaries, run receipts, private infrastructure configuration and credentials.
Reusable prompt templates are permitted. Completed prompts from research runs remain excluded.
The original repository and its history remain private.

## Check a change

1. Review each new file and each changed line.
2. Keep runtime outputs in an ignored directory.
3. Run `make harness-tests`.
4. Run `make privacy-scan-repo`.
5. Inspect `git diff --cached` before you commit.

The scope check rejects tracked output directories, binary artifacts and common credential patterns.
The genomic record check permits only named synthetic fixtures and public source examples.
These checks cannot determine the origin of arbitrary prose. Review remains necessary before publication.

## Reproduction limits

The synthetic tests run without the research dataset.
Full genomic runs require separately supplied reads, references and execution resources.
Live literature queries require network access and can change when the source services change.
The code does not reconstruct the private submission from this repository alone.
