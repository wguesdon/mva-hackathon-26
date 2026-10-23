"""Measure the peak resident memory of a child process and every descendant of it.

No run of this project recorded its memory, so a reader could not tell whether the method needs a
laptop or a server. `resource.getrusage` reports a high water mark that never falls, so a value
taken after one stage also covers every earlier stage of the same process tree. This module samples
`/proc/<pid>/status` instead, which gives the peak of one stage alone.

The sampler reads `VmRSS` of the stage process and of every descendant at a fixed interval and
keeps the largest sum. A process that lives for less than one interval is not seen, so the value is
a lower bound on a very short stage. Every stage of the Track 2 pipeline that allocates memory runs
for more than one second, which is five intervals.

The module works on Linux, because `/proc` is a Linux interface. On a system with no `/proc` the
sampler reports `None` and the caller records that the measurement did not run.

Usage:
    sampler = start_sampler(process.pid)
    ...
    peak_kb = stop_sampler(sampler)
"""

from __future__ import annotations

import threading
import time
from pathlib import Path

PROC = Path("/proc")
SAMPLE_INTERVAL_S = 0.2


def available() -> bool:
    """Return whether this system exposes the `/proc` interface.

    Returns:
        True when `/proc/self/status` is readable.
    """
    return (PROC / "self" / "status").is_file()


def read_rss_kb(pid: int) -> int:
    """Read the resident set size of one process.

    Args:
        pid: The process identifier.

    Returns:
        The resident set size in kibibytes, or 0 when the process has already exited.
    """
    try:
        for line in (PROC / str(pid) / "status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        return 0
    return 0


def children_of(pid: int) -> list[int]:
    """Read the direct children of one process.

    Args:
        pid: The process identifier.

    Returns:
        The process identifiers of the direct children.
    """
    found: list[int] = []
    try:
        for task in (PROC / str(pid) / "task").iterdir():
            try:
                text = (task / "children").read_text()
            except OSError:
                continue
            found.extend(int(value) for value in text.split())
    except OSError:
        return found
    return found


def tree_rss_kb(pid: int) -> int:
    """Sum the resident set size of one process and every descendant of it.

    Args:
        pid: The process identifier of the root of the tree.

    Returns:
        The sum in kibibytes.
    """
    total = 0
    pending = [pid]
    seen: set[int] = set()
    while pending:
        current = pending.pop()
        if current in seen:
            continue
        seen.add(current)
        total += read_rss_kb(current)
        pending.extend(children_of(current))
    return total


def sample_until_stopped(pid: int, stop: threading.Event, result: dict) -> None:
    """Sample the memory of a process tree until the caller sets the stop event.

    Args:
        pid: The process identifier of the root of the tree.
        stop: The event that ends the sampling.
        result: The mapping that receives the key `peak_kb`.
    """
    peak = 0
    while not stop.is_set():
        peak = max(peak, tree_rss_kb(pid))
        stop.wait(SAMPLE_INTERVAL_S)
    peak = max(peak, tree_rss_kb(pid))
    result["peak_kb"] = peak


def start_sampler(pid: int) -> dict:
    """Start the sampler of one process tree.

    Args:
        pid: The process identifier of the root of the tree.

    Returns:
        The sampler handle, which holds the thread, the stop event and the result mapping. The
        handle carries `supported` set to False when this system has no `/proc`.
    """
    if not available():
        return {"supported": False, "result": {}}
    stop = threading.Event()
    result: dict = {}
    thread = threading.Thread(target=sample_until_stopped, args=(pid, stop, result), daemon=True)
    thread.start()
    return {"supported": True, "thread": thread, "stop": stop, "result": result}


def stop_sampler(handle: dict) -> float | None:
    """Stop the sampler and return the peak memory in mebibytes.

    Args:
        handle: The handle that `start_sampler` returned.

    Returns:
        The peak resident memory of the process tree in mebibytes, or None when the measurement did
        not run on this system.
    """
    if not handle.get("supported"):
        return None
    handle["stop"].set()
    handle["thread"].join(timeout=2.0 + SAMPLE_INTERVAL_S)
    peak_kb = handle["result"].get("peak_kb", 0)
    return round(peak_kb / 1024.0, 1)


def main() -> int:
    """Print the peak memory of a short command, as a check of the sampler.

    Returns:
        0 on success.
    """
    import subprocess  # noqa: PLC0415

    command = "import time; x = bytearray(200 * 1024 * 1024); time.sleep(1.5); del x"
    process = subprocess.Popen(["uv", "run", "python", "-c", command])
    handle = start_sampler(process.pid)
    process.wait()
    print(f"peak of the process tree: {stop_sampler(handle)} MiB, expected above 200 MiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
