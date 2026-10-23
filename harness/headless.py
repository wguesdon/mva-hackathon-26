"""Run read-only reviewers through headless coding-agent CLIs.

Ported from kaggriculture_harness/harness/headless.py. The MVA port keeps three backends:
Claude Code and Codex on the owner's subscriptions, and OpenCode for Amazon Bedrock models.
A reviewer reads files inside this repository only. No backend has a write or shell tool.
"""

from __future__ import annotations

import asyncio
import json
import os
import shutil
import signal
import uuid
from pathlib import Path

if __package__:
    from . import paths
else:
    import paths

BACKENDS = ("claude", "codex", "opencode")


def command(backend: str, model: str, prompt: str) -> tuple[list[str], dict[str, str]]:
    """Build a noninteractive review command without a shell.

    Args:
        backend: The CLI name.
        model: The explicit model identifier accepted by that CLI.
        prompt: The complete review prompt.

    Returns:
        The argument vector and environment overrides.
    """
    if not model.strip() or model.startswith("-"):
        raise ValueError("A headless assignment requires an explicit model.")
    if len(prompt.encode()) > 100_000:
        raise ValueError("The headless prompt exceeds 100000 bytes. Use a smaller task context.")
    if backend == "claude":
        # Restricted mode confines the file tools to the working directory.
        # The repository holds no patient data, so a reviewer gets no --add-dir.
        return ["claude", "-p", prompt, "--model", model, "--output-format", "stream-json", "--verbose",
                "--safe-mode", "--restricted", "--strict-mcp-config", "--permission-mode", "dontAsk",
                "--permission-prompts", "none", "--allowedTools", "Read,Glob,Grep",
                "--tools", "Read,Glob,Grep", "--no-session-persistence"], {}
    if backend == "codex":
        return ["codex", "exec", "--model", model, "--sandbox", "read-only", "--json",
                "--ephemeral", "--ignore-user-config", "--ignore-rules", "--skip-git-repo-check",
                "-c", "approval_policy=\"never\"", prompt], {}
    if backend == "opencode":
        # OpenCode recursively merges agent maps. A fresh name prevents inherited
        # agent-specific allow rules from surviving our deny-by-default policy.
        reviewer = "research_reviewer_" + uuid.uuid4().hex
        # The reviewer reads this repository only.
        permissions = {"*": "deny", "read": "allow", "grep": "allow", "glob": "allow", "list": "allow",
                       "external_directory": "deny"}
        settings = {"permission": permissions, "share": "disabled", "agent": {
            reviewer: {"mode": "primary", "permission": permissions,
                       "prompt": "Review the assigned research task. Use only read-only file tools."}}}
        return ["opencode", "run", "--pure", "--agent", reviewer, "--model", model, "--format", "json", prompt], {
            "OPENCODE_CONFIG_CONTENT": json.dumps(settings),
            "OPENCODE_PERMISSION": json.dumps(permissions),
            "OPENCODE_DISABLE_PROJECT_CONFIG": "true",
        }
    raise ValueError("The headless backend is unsupported.")


def parse_output(backend: str, output: str, exit_code: int) -> dict:
    """Require a final answer and a successful terminal event.

    Args:
        backend: The CLI name.
        output: The captured JSON or JSONL output.
        exit_code: The process exit code.

    Returns:
        The answer, completion status, and failure reason when needed.
    """
    if exit_code != 0:
        return {"status": "FAILED", "answer": "", "reason": f"The process exited with code {exit_code}."}
    try:
        # Claude also accepts a historical single-result JSON fixture.
        if backend == "claude" and output.lstrip().startswith("{"):
            try:
                events = [json.loads(output)]
            except json.JSONDecodeError:
                events = [json.loads(line) for line in output.splitlines() if line.strip()]
        else:
            events = [json.loads(line) for line in output.splitlines() if line.strip()]
    except (json.JSONDecodeError, TypeError):
        return {"status": "INCOMPLETE", "answer": "", "reason": "The CLI returned invalid JSON."}
    answer, complete, failed = [], False, False
    opencode_step = []
    models = []
    for event in events:
        if not isinstance(event, dict):
            return {"status": "INCOMPLETE", "answer": "", "reason": "The CLI event is not an object."}
        kind = event.get("type")
        if kind in ("error", "turn.failed") or event.get("is_error"):
            failed = True
        if backend == "claude":
            if kind == "system" and event.get("model"):
                models.append(event["model"])
            if kind == "result":
                complete = event.get("subtype") == "success" and not event.get("is_error")
                answer.append(event.get("result", ""))
                models.extend(event.get("modelUsage", {}))
        elif backend == "codex":
            complete = complete or kind == "turn.completed"
            item = event.get("item", {})
            if kind == "item.completed" and item.get("type") == "agent_message":
                # Codex emits progress as complete agent messages too. The last
                # message is the answer; the transcript retains earlier messages.
                answer = [item.get("text", "")]
        elif backend == "opencode":
            part = event.get("part", {})
            if kind == "text":
                opencode_step.append(part.get("text", ""))
            if kind == "step_finish":
                if part.get("reason") == "stop":
                    answer = list(opencode_step)
                    complete = True
                opencode_step = []
    text = "\n".join(value for value in answer if isinstance(value, str)).strip()
    status = "FAILED" if failed else "COMPLETE" if complete and text else "INCOMPLETE"
    return {"status": status, "answer": text, "reported_models": sorted(set(models))}


def parse_transcript(backend: str, path: Path, exit_code: int) -> dict:
    """Read terminal events without retaining verbose tool output in memory.

    Args:
        backend: The CLI name.
        path: The JSONL transcript path.
        exit_code: The process exit code.

    Returns:
        The parsed completion record.
    """
    retained = []
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                return {"status": "INCOMPLETE", "answer": "", "reason": "The transcript contains invalid JSON."}
            if not isinstance(event, dict):
                continue
            kind = event.get("type")
            if kind in ("result", "error", "turn.completed", "turn.failed", "message_end", "text", "step_finish"):
                retained.append(line)
            elif kind == "system" and event.get("model"):
                retained.append(json.dumps({"type": "system", "model": event["model"]}))
            elif kind == "item.completed" and event.get("item", {}).get("type") == "agent_message":
                retained.append(line)
    return parse_output(backend, "\n".join(retained), exit_code)


def credential_leaves(value: object) -> list[str]:
    """Collect credential fields from a CLI authentication document.

    Args:
        value: A decoded authentication value.

    Returns:
        Credential strings, without account identifiers or expiry fields.
    """
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            if isinstance(item, str) and len(item) >= 8 and (key.lower() in ("access", "refresh") or any(
                word in key.lower() for word in ("key", "token", "password", "secret")
            )):
                found.append(item)
            elif isinstance(item, (dict, list)):
                found.extend(credential_leaves(item))
    elif isinstance(value, list):
        for item in value:
            found.extend(credential_leaves(item))
    return found


def secret_values() -> list[str]:
    """Read known credentials only to prevent their persistence in review artifacts.

    Returns:
        Nonempty credential values from the environment and configured CLI stores.
    """
    values = credential_leaves(dict(os.environ))
    stores = (Path.home() / ".local/share/opencode/auth.json",
              Path.home() / ".claude/.credentials.json",
              Path(os.environ.get("CODEX_HOME", str(Path.home() / ".codex"))) / "auth.json")
    for path in stores:
        if path.is_file():
            values.extend(credential_leaves(json.loads(path.read_text())))
    for path in paths.ENV_FILES:
        if path.is_file():
            for line in path.read_text().splitlines():
                if "=" in line and not line.lstrip().startswith("#"):
                    key, value = line.split("=", 1)
                    values.extend(credential_leaves({key: value.strip().strip("\"'")}))
    return sorted(set(values))


def reviewer_environment(backend: str, model: str) -> dict[str, str]:
    """Remove unrelated credentials while preserving the selected provider's authentication.

    Args:
        backend: The selected CLI.
        model: The selected model identifier.

    Returns:
        A filtered copy of the process environment.
    """
    provider = {"claude": "ANTHROPIC", "codex": "OPENAI"}.get(backend, model.split("/")[0].upper().replace("-", "_"))
    bedrock = backend == "opencode" and model.startswith("amazon-bedrock/")
    environment = dict(os.environ)
    for key in list(environment):
        sensitive = any(word in key.upper() for word in ("KEY", "TOKEN", "SECRET", "PASSWORD", "CREDENTIAL"))
        selected_auth = key.startswith(provider + "_") or backend == "claude" and key == "CLAUDE_CODE_OAUTH_TOKEN"
        # A Bedrock reviewer authenticates through the AWS profile files. It keeps the profile
        # name and the region, but never a key in the environment.
        bedrock_setting = bedrock and key in ("AWS_PROFILE", "AWS_REGION", "AWS_DEFAULT_REGION")
        if bedrock_setting:
            continue
        if key.startswith(("KAGGLE", "AWS_", "HF_", "HUGGING")) or sensitive and not selected_auth:
            del environment[key]
    if bedrock:
        environment.setdefault("AWS_REGION", "us-east-1")
    return environment


def sanitize_transcript(path: Path, secrets: list[str]) -> bool:
    """Redact known credentials from a completed log and report any exposure.

    Args:
        path: The process log path.
        secrets: Known credential values.

    Returns:
        True if the transcript contained a credential.
    """
    text = path.read_text(errors="replace")
    found = False
    for secret in secrets:
        if secret and secret in text:
            text = text.replace(secret, "[REDACTED]")
            found = True
    if found:
        path.write_text(text)
    return found


async def run(backend: str, model: str, prompt: str, root: Path, output: Path, deadline: float) -> dict:
    """Run one reviewer with a wall-clock timeout and captured process logs.

    Args:
        backend: The CLI name.
        model: The explicit model identifier.
        prompt: The complete prompt.
        root: The review working directory.
        output: The output filename stem.
        deadline: The process deadline in seconds.

    Returns:
        A completion record with the answer and log references.
    """
    argv, overrides = command(backend, model, prompt)
    if shutil.which(argv[0]) is None:
        return {"status": "FAILED", "answer": "", "reason": f"The {backend} executable is absent."}
    secrets = secret_values()
    environment = reviewer_environment(backend, model)
    environment.update(overrides)
    stdout, stderr = output.with_suffix(".jsonl"), output.with_suffix(".err")
    interrupted = None
    with stdout.open("wb") as out, stderr.open("wb") as err:
        process = await asyncio.create_subprocess_exec(
            *argv, cwd=root, env=environment, stdin=asyncio.subprocess.DEVNULL,
            stdout=out, stderr=err, start_new_session=True,
        )
        try:
            await asyncio.wait_for(process.wait(), timeout=deadline)
        except (TimeoutError, asyncio.CancelledError) as error:
            interrupted = "The reviewer exceeded its deadline." if isinstance(error, TimeoutError) else "The reviewer received a stop request."
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                await asyncio.wait_for(process.wait(), timeout=10)
            except (TimeoutError, asyncio.CancelledError):
                return {"status": "FAILED", "answer": "", "pid": process.pid,
                        "reason": "The reviewer may remain live after SIGTERM. Keep the campaign lock and inspect its private logs."}
    secrets = sorted(set(secrets + secret_values()))
    exposed = sanitize_transcript(stdout, secrets)
    exposed = sanitize_transcript(stderr, secrets) or exposed
    if exposed:
        return {"status": "FAILED", "answer": "", "reason": "A known credential appeared in the reviewer output. The log was redacted."}
    if interrupted:
        return {"status": "INCOMPLETE", "answer": "", "reason": interrupted,
                "stdout": str(stdout), "stderr": str(stderr)}
    result = parse_transcript(backend, stdout, process.returncode)
    result.update(stdout=str(stdout), stderr=str(stderr))
    return result
