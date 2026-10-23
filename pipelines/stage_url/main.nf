#!/usr/bin/env nextflow
// Download one public URL on an AWS Batch worker and publish it to S3. Class A reference data only.

nextflow.enable.dsl = 2

params.url        = null
params.outdir     = null
params.sha256_url = null

process DOWNLOAD {
    // python:3.12-slim has no ps, which Nextflow needs for task metrics.
    container 'docker.io/library/python:3.12-bookworm'
    cpus 2
    memory '4 GB'
    disk '80 GB'
    publishDir params.outdir, mode: 'copy'

    input:
    val url

    output:
    path "*", includeInputs: false

    script:
    def name = url.tokenize('/')[-1]
    """
    python3 - <<'PY'
import hashlib, shutil, urllib.request
name = "${name}"
with urllib.request.urlopen("${url}", timeout=120) as response, open(name, "wb") as handle:
    shutil.copyfileobj(response, handle, length=16 * 1024 * 1024)
digest = hashlib.sha256()
with open(name, "rb") as handle:
    for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
        digest.update(block)
open(name + ".sha256", "w").write(digest.hexdigest() + "  " + name + "\\n")
PY
    """
}

workflow {
    if (!params.url || !params.outdir) {
        error "Set url and outdir."
    }
    DOWNLOAD(params.url)
}
