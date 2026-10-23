# Methods and execution

## Agent coordination

`harness/research.py` stores task state and validates project paths.
`harness/blackboard.py` manages reviewer exchanges and a configured call budget.
`harness/swarm.py` executes assignments through the CLI adapters in `harness/headless.py`.
`harness/privacy_guard.py` checks outgoing content against declared data classes and configured deny tokens.

Use each command's help output to inspect its arguments.

```bash
uv run python harness/research.py --help
uv run python harness/blackboard.py --help
uv run python harness/swarm.py --help
```

The example limits are in `config.yaml`.
Copy `.env.example` to `.env` for local deny-token configuration.
Do not commit the completed file or model transcripts.

## Track 1 code

The calling configuration example uses nf-core/raredisease 3.1.2 with GRCh38 and DeepVariant.
The caller must supply reads, a sample sheet and reference paths.
Nextflow workflows under `pipelines/` cover annotation, forced genotypes, structural variants, local reassembly, phase and benchmark comparison.
Those workflows require Nextflow and their specified containers.

`pipelines/track1_b0/bin/b0_rank.py` scores annotated variants from consequence, population frequency and clinical classification.
It constructs same-gene compound heterozygous hypotheses and homozygous hypotheses.
Separate workflows assess phase and structural evidence.
The merge script preserves rank order and separates tied confidence values within the submission row budget.

```bash
uv run python pipelines/track1_b0/bin/b0_rank.py --help
uv run python scripts/merge_track1_submission.py --help
uv run python scripts/check_merged_submission.py --help
```

The source contains default reference filenames from the original analysis.
Supply current local input paths through the command arguments.
Generated CSVs, genotypes and run records must remain outside Git.

## Track 2 code

The graph runner resolves a disease identifier, retrieves public relationships and constructs candidate paths through genes and pathways.
The ranker applies path scores and drug eligibility rules.
The withheld-indication control removes known indication edges and tests whether the ranker recovers their drugs.
The permutation control compares the observed ranking with graphs whose node degrees are preserved during edge changes.

The resource manifest is `literature/track2_runner_resources.yaml`.
The cost model contains coefficients fitted on public disease graphs and does not contain individual run records.
The runner checks downloaded reference checksums before a full run.

```bash
uv run python scripts/track2_run_disease.py --help
uv run python scripts/track2_controls.py --help
uv run python scripts/literature_pack.py --help
```

`make demo-track2` exercises the graph interface on invented inputs.
A successful synthetic run confirms software behavior on that fixture. It does not validate a treatment hypothesis.

## Report creation

`writeups/_extensions/` contains the paper and REVTeX formats.
`writeups/Containerfile` specifies Quarto 1.8.27 and its PDF dependencies.
The blank document at `writeups/example/report_template.qmd` contains no research results.

```bash
podman build -t localhost/mva-writeups:quarto-1.8.27 -f writeups/Containerfile writeups
podman run --rm --userns=keep-id -v "$PWD/writeups:/work:Z" -e HOME=/tmp localhost/mva-writeups:quarto-1.8.27 render example/report_template.qmd --to paper-pdf
```

Git ignores the generated PDF. Completed reports and submission packages belong in the private repository.

## Inputs not distributed

The synthetic tests and graph example run from this checkout.
The commands below require additional inputs. The export does not claim to reproduce the original reports without them.

| Command or component | Required inputs and preparation |
|---|---|
| Nextflow genomic workflows | Supply authorized sample files, references and container execution resources. |
| `build_s6_targets.py` | Supply ClinVar and a resolved public literature allele table with `--alleles`. The default dated output is not included. |
| `build_spikein_suite.py` | Supply public benchmark variants and annotation resources. The default cached HG002 call file is not included. |
| `summarize_s6_patient.py` | Supply the requested sample results and a control summary. The default historical control receipt is not included. |
| `track2_run_disease.py` | Fetch the public reference releases listed in the resource manifest. The run then queries public APIs. |
| `t2_regulator_screen_v1.py` | Build a graph first for the graph-presence comparison, or supply its expected `reports/track2/graph/nodes.tsv` input locally. |
| `t2_endpoint_search_complete_2026_10_22.py` | Create and commit a new local registration document before the query. Pass its path and the new commit through `--registration` and `--registration-commit`. No private commit is required. |
| `t2_preclinical_protocol.py` | Supply a protocol parameter YAML, evidence TSV, bibliography check status and supporting experimental-design records. Completed protocol configurations are excluded. |
| `t2_nad_pharmacokinetics.py` and `t2_nad_exposure_compartments_v1.py` | Supply the evidence tables described by their command help. Completed tables are excluded. |
| `track2_functional_compare.py` | Supply the graph layers and reference control receipts. Historical control outputs are excluded. |
| `check_report_claims.py` | Supply a manuscript, its claim configuration and its referenced receipts. This export includes no completed manuscript. |

Historical source comments can name tests or records that belong to the private research archive.
`make harness-tests` runs only the test scripts present in this checkout.
The gene-panel consistency test mentioned in some source comments is not part of this exported suite.
The query registration check accepts a newly committed plan in this repository; it does not require the research repository's history.
