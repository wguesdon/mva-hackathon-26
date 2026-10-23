"""Send one class A question to a second and a third model, then record both answers.

The script exists so that a decision does not rest on one model. It sends the same prompt to Amazon
Bedrock and to OpenRouter, and it writes both answers to `reports/second_opinion/`.

The prompt must be class A. The privacy guard scans the prompt file first, and the script stops when
the guard finds a variant record, a genotype or a read. OpenRouter is approved for class A only, see
`docs/llm_data_privacy.md`.

Usage:
    uv run python scripts/ask_second_opinion.py --prompt prompts/class_a/question.md --slug llm_terms
    uv run python scripts/ask_second_opinion.py --prompt p.md --slug s --skip openrouter
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import httpx

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "harness"))
import privacy_guard  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
# The `us.` inference profile keeps the request inside the United States. `global.` does not.
BEDROCK_MODELS = {"bedrock_claude": "us.anthropic.claude-opus-5", "bedrock_gptoss": "openai.gpt-oss-120b-1:0"}
BEDROCK_REGION = "us-east-1"
OPENROUTER_MODEL = "openai/gpt-5.1"
OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"
OPENROUTER_KEY_PARAMETER = "/openclaw/openrouter-api-key"
# Claude Opus 5 returns a reasoningContent block beside the text block and both count against
# maxTokens, so a low cap truncates the answer while the token count still looks healthy.
MAX_TOKENS = 32000
READ_TIMEOUT_SECONDS = 900
PROVIDERS = ("bedrock_claude", "bedrock_gptoss", "openrouter")


def guard(prompt: str) -> None:
    """Stop the run when the prompt holds class B or class C content.

    Args:
        prompt: The prompt text.

    Raises:
        SystemExit: When the privacy guard reports a problem.
    """
    problems = privacy_guard.scan(prompt)
    if problems:
        for problem in problems:
            print(f"privacy guard: {problem}", file=sys.stderr)
        raise SystemExit("The prompt is not class A. Nothing was sent.")


def ssm_parameter(name: str) -> str:
    """Read one SecureString parameter.

    Args:
        name: The parameter name.

    Returns:
        The decrypted value.
    """
    result = subprocess.run(["aws", "ssm", "get-parameter", "--name", name, "--with-decryption",
                             "--query", "Parameter.Value", "--output", "text"],
                            capture_output=True, text=True, check=True)
    return result.stdout.strip()


def ask_bedrock(prompt: str, model: str) -> str:
    """Ask a Bedrock model through the AWS CLI.

    Args:
        prompt: The prompt text.
        model: The Bedrock model identifier or inference profile.

    Returns:
        The answer text.
    """
    # Claude Opus 5 rejects `temperature`, so only the token cap is set.
    body = {"messages": [{"role": "user", "content": [{"text": prompt}]}],
            "inferenceConfig": {"maxTokens": MAX_TOKENS}}
    with tempfile.TemporaryDirectory() as folder:
        # converse writes the answer to stdout; it takes no output file argument.
        request = Path(folder) / "request.json"
        request.write_text(json.dumps(body))
        # The CLI read timeout defaults to 60 seconds, which a long answer exceeds.
        result = subprocess.run(["aws", "bedrock-runtime", "converse", "--region", BEDROCK_REGION,
                                 "--cli-read-timeout", str(READ_TIMEOUT_SECONDS),
                                 "--model-id", model, "--cli-input-json", f"file://{request}"],
                                capture_output=True, text=True, check=True)
        answer = json.loads(result.stdout)
    blocks = answer["output"]["message"]["content"]
    if answer.get("stopReason") == "max_tokens":
        print(f"warning: {model} stopped at the token cap, the answer is truncated", file=sys.stderr)
    return "\n".join(block.get("text", "") for block in blocks).strip()


def ask_openrouter(prompt: str, model: str) -> str:
    """Ask a model through OpenRouter with logging off.

    Args:
        prompt: The prompt text.
        model: The OpenRouter model identifier.

    Returns:
        The answer text.
    """
    key = os.environ.get("OPENROUTER_API_KEY") or ssm_parameter(OPENROUTER_KEY_PARAMETER)
    payload = {"model": model, "max_tokens": MAX_TOKENS, "temperature": 0.2,
               "messages": [{"role": "user", "content": prompt}]}
    headers = {"Authorization": f"Bearer {key}", "Content-Type": "application/json",
               "X-Title": "mva-hackathon-2026"}
    with httpx.Client(timeout=300) as client:
        response = client.post(OPENROUTER_URL, headers=headers, json=payload)
        response.raise_for_status()
        return response.json()["choices"][0]["message"]["content"].strip()


def run_one(name: str, prompt: str) -> tuple[str, str]:
    """Call one provider and keep the error text on a failure.

    Args:
        name: A key of BEDROCK_MODELS, or `openrouter`.
        prompt: The prompt text.

    Returns:
        The model identifier and the answer or the error text.
    """
    caller = ask_openrouter if name == "openrouter" else ask_bedrock
    model = OPENROUTER_MODEL if name == "openrouter" else BEDROCK_MODELS[name]
    try:
        return model, caller(prompt, model)
    except subprocess.CalledProcessError as error:
        return model, f"CALL FAILED: {error.stderr.strip()[:2000]}"
    except httpx.HTTPError as error:
        return model, f"CALL FAILED: {error}"


def main() -> int:
    """Ask both providers and write the report.

    Returns:
        Zero when at least one provider answered.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--prompt", type=Path, required=True)
    parser.add_argument("--slug", required=True)
    parser.add_argument("--skip", action="append", default=[], choices=PROVIDERS)
    parser.add_argument("--attach", action="append", default=[], type=Path,
                        help="A class A file to append to the prompt. Repeat the flag for more than one.")
    args = parser.parse_args()
    parts = [args.prompt.read_text()]
    for path in args.attach:
        parts.append(f"\n\n# Attached file: {path.name}\n\n{path.read_text()}")
    prompt = "".join(parts)
    guard(prompt)
    print(f"prompt characters: {len(prompt)} from {1 + len(args.attach)} file(s)")
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    sections, answered = [], False
    for name in PROVIDERS:
        if name in args.skip:
            continue
        model, answer = run_one(name, prompt)
        answered = answered or not answer.startswith("CALL FAILED")
        sections.append(f"## {name}, `{model}`\n\n{answer}\n")
        print(f"{name}: {len(answer)} characters")
    out_dir = ROOT / "reports" / "second_opinion"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{stamp}_{args.slug}.md"
    header = [f"# Second opinion: {args.slug}", "",
              f"Asked on {stamp}. Data class A. The prompt is `{args.prompt}`.", "",
              "The answers below come from a model. Check every quotation against the primary source.", ""]
    path.write_text("\n".join(header) + "\n".join(sections) + "\n")
    print(f"Wrote {path}")
    return 0 if answered else 1


if __name__ == "__main__":
    sys.exit(main())
