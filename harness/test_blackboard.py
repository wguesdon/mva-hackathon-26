"""Verify forum delivery, the privacy guard, context packing and budgets without model calls."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import blackboard
import privacy_guard
import research


def response_for(prompt: str) -> str:
    """Build a synthetic reviewer response that replies to peer posts.

    Args:
        prompt: The exact dispatched forum prompt.

    Returns:
        A valid JSON response.
    """
    context = json.loads(prompt.split("\n\nForum context:\n", 1)[1])
    participant = context["participant"]
    peer_posts = [post for post in context["posts"] if post["sender"] != participant]
    parent = peer_posts[-1] if peer_posts else None
    last = context["round"] == context["rounds"]
    kind = "proposal" if context["round"] == 1 else "decision" if last else "reply"
    body = f"The {participant} reviewer proposes a spike-in benchmark. " + "Detail. " * 300
    to = ["all"]
    reply_to = parent["post_id"] if parent and not last else None
    return json.dumps({"messages": [{"thread": "strategy", "to": to, "reply_to": reply_to,
                                     "kind": kind, "body": body[:3000]}]})


async def fake_review(backend: str, model: str, prompt: str, root: Path, stem: Path, deadline: float) -> dict:
    """Return a valid synthetic final answer through the real swarm adapter.

    Args:
        backend: The configured CLI backend.
        model: The configured model.
        prompt: The exact dispatched prompt.
        root: The synthetic campaign root.
        stem: The output file stem.
        deadline: The configured response limit.

    Returns:
        A successful synthetic CLI response.
    """
    return {"status": "COMPLETE", "answer": response_for(prompt)}


async def synthetic_restricted_record_review(backend: str, model: str, prompt: str, root: Path, stem: Path, deadline: float) -> dict:
    """Return an answer that holds a VCF record.

    Args:
        backend: The configured CLI backend.
        model: The configured model.
        prompt: The exact dispatched prompt.
        root: The synthetic campaign root.
        stem: The output file stem.
        deadline: The configured response limit.

    Returns:
        A response that the privacy scan must refuse.
    """
    response = json.loads(response_for(prompt))
    response["messages"][0]["body"] = "chr15\t40161100\t.\tG\tA\t50\tPASS\t.\tGT:AD\t0/1:1,2\n"
    return {"status": "COMPLETE", "answer": json.dumps(response)}


class BlackboardTest(unittest.TestCase):
    """Exercise the forum against the real registry and swarm runner."""

    def setUp(self) -> None:
        """Create a temporary root with two phase four tasks."""
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        for task_id, assignee in (("calling", "claude"), ("ranking", "codex")):
            research.create_task(self.root, {
                "task_id": task_id, "phase": 4, "track": "t1_strategy", "assignee": assignee,
                "question": f"Propose the {task_id} strategy.", "criterion": "A testable plan.",
                "control": "none", "target": "task1_strategy", "suite_hash": "none", "data_class": "A"})
        self.assignments = self.root / "assignments.json"
        self.assignments.write_text(json.dumps([
            {"task_id": "calling", "backend": "claude", "model": "test_claude"},
            {"task_id": "ranking", "backend": "codex", "model": "test_codex"},
        ]))
        self.goal = self.root / "goal.md"
        self.goal.write_text("Define a Task 1 strategy from public evidence only.")
        for target, value in ((blackboard.headless, "secret_values"), (privacy_guard, "deny_tokens")):
            patch = mock.patch.object(target, value, return_value=[])
            patch.start()
            self.addCleanup(patch.stop)

    def run_forum(self) -> dict:
        """Execute or resume the forum under the campaign lock.

        Returns:
            The resulting forum state.
        """
        with research.CampaignLock(self.root, "test_forum") as lock:
            return blackboard.run_forum(self.root, "strategy_board", lock)

    def test_full_exchange_logs_every_call(self) -> None:
        """A complete forum records one SENT exposure row per call and two decision findings."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 4, 8, {})
        with mock.patch.object(blackboard.headless, "run", side_effect=fake_review) as call:
            board = self.run_forum()
        self.assertEqual(board["status"], "COMPLETE")
        self.assertEqual(call.call_count, 8)
        self.assertEqual(len(board["findings"]), 2)
        with (self.root / privacy_guard.CALL_LOG).open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["verdict"] == "SENT" and row["data_class"] == "A" for row in rows))

    def test_budget_limit_from_config(self) -> None:
        """A board above the configured call limit is refused."""
        with self.assertRaises(ValueError):
            blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 4, 51,
                                    {"board_max_calls": 50})

    def test_class_b_task_refused_for_subscription_arm(self) -> None:
        """The exported board enforces its configured class A ceiling for Claude Code."""
        research.create_task(self.root, {
            "task_id": "private_review", "phase": 4, "track": "t1_strategy", "assignee": "claude",
            "question": "Review candidates.", "criterion": "A ranked list.", "control": "none",
            "target": "task1_strategy", "suite_hash": "none", "data_class": "B"})
        path = self.root / "b_assignments.json"
        path.write_text(json.dumps([{"task_id": "private_review", "backend": "claude", "model": "m"},
                                    {"task_id": "ranking", "backend": "codex", "model": "m"}]))
        with self.assertRaises(ValueError):
            blackboard.create_board(self.root, "b_board", path, self.goal, 3, 6, {})

    def test_synthetic_restricted_record_is_not_posted(self) -> None:
        """An answer with a VCF record fails and stops the forum for resume."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 3, 6, {})
        with mock.patch.object(blackboard.headless, "run", side_effect=synthetic_restricted_record_review):
            board = self.run_forum()
        self.assertEqual(board["status"], "NEEDS_RESUME")
        self.assertEqual(board["posts"], {})

    def test_refused_prompt_is_logged_and_not_dispatched(self) -> None:
        """A goal with a private identifier is refused before any dispatch."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 3, 6, {})
        with mock.patch.object(privacy_guard, "deny_tokens", return_value=["Task 1 strategy"]), \
                mock.patch.object(blackboard.headless, "run", side_effect=fake_review) as call:
            board = self.run_forum()
        self.assertEqual(call.call_count, 0)
        self.assertEqual(board["status"], "NEEDS_RESUME")
        with (self.root / privacy_guard.CALL_LOG).open() as handle:
            self.assertTrue(all(row["verdict"] == "REFUSED" for row in csv.DictReader(handle, delimiter="\t")))

    def test_body_limit(self) -> None:
        """A 6000-byte body passes and a 9000-byte body fails."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 3, 6, {})
        board = blackboard.get_board(self.root, "strategy_board")
        attempt = {"round": 1, "task_id": "calling", "delivered_ids": []}
        for size, ok in ((6000, True), (9000, False)):
            text = json.dumps({"messages": [{"thread": "t", "to": ["all"], "reply_to": None,
                                             "kind": "proposal", "body": "a" * size}]})
            if ok:
                self.assertEqual(len(blackboard.parse_messages(board, attempt, text)), 1)
            else:
                with self.assertRaises(ValueError):
                    blackboard.parse_messages(board, attempt, text)

    def test_long_decision_passes(self) -> None:
        """A 12000-byte decision passes and a 12000-byte proposal fails."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 3, 6, {})
        board = blackboard.get_board(self.root, "strategy_board")
        for kind, round_number, ok in (("decision", 3, True), ("proposal", 1, False)):
            attempt = {"round": round_number, "task_id": "calling", "delivered_ids": []}
            text = json.dumps({"messages": [{"thread": "t", "to": ["all"], "reply_to": None,
                                             "kind": kind, "body": "a" * 12000}]})
            if ok:
                self.assertEqual(len(blackboard.parse_messages(board, attempt, text)), 1)
            else:
                with self.assertRaises(ValueError):
                    blackboard.parse_messages(board, attempt, text)

    def test_resume_accepts_one_answer_per_round(self) -> None:
        """Two saved failing answers for one round yield one post after the limit is fixed."""
        blackboard.create_board(self.root, "strategy_board", self.assignments, self.goal, 3, 12, {})
        with mock.patch.object(blackboard, "MAX_BODY", 100), \
                mock.patch.object(blackboard.headless, "run", side_effect=fake_review):
            self.run_forum()
            self.run_forum()
        state = blackboard.get_board(self.root, "strategy_board")
        self.assertEqual(state["posts"], {})
        self.assertEqual(state["calls"], 4)
        with mock.patch.object(blackboard.headless, "run", side_effect=fake_review) as call:
            board = self.run_forum()
        self.assertEqual(board["status"], "COMPLETE")
        round_one = [post for post in board["posts"].values() if post["round"] == 1]
        self.assertEqual(len(round_one), 2)
        self.assertEqual(call.call_count, 4)
        statuses = sorted(item["status"] for item in board["attempts"].values() if item["round"] == 1)
        self.assertEqual(statuses, ["COMPLETE", "COMPLETE", "SUPERSEDED", "SUPERSEDED"])

    def test_pack_posts_keeps_unread_and_stubs_old(self) -> None:
        """A small budget keeps unread bodies in full and shortens older posts."""
        board = {"board_id": "b"}
        posts = [{"post_id": f"post_{index}", "thread": "t", "sender": "calling", "to": ["all"], "kind": "reply",
                  "round": index, "reply_to": None, "body": "x " * 1500, "created_utc": str(index),
                  "run_id": "r", "artifacts": []} for index in range(6)]
        packed = blackboard.pack_posts(board, posts, {"post_5"}, 6000)
        self.assertEqual(packed[5]["body"], posts[5]["body"])
        self.assertIn("full_text_path", packed[0])
        self.assertLessEqual(len(json.dumps(packed).encode()), 6000)


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=1).result
    sys.exit(0 if result.wasSuccessful() else 1)
