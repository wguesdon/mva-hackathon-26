PY := uv run python

.PHONY: harness-tests privacy-scan-repo demo-track2

harness-tests:
	cd harness && $(PY) run_tests.py

privacy-scan-repo:
	$(PY) scripts/check_repository_scope.py
	@git ls-files -z | xargs -0 -n 50 $(PY) harness/privacy_guard.py --allowlist literature/privacy_guard_allowlist.yaml --quiet-pass

demo-track2:
	$(PY) scripts/track2_make_graph_fixture.py --out-dir reports/track2/fixture_graph
	$(PY) scripts/track2_rank.py --nodes reports/track2/fixture_graph/nodes.tsv --edges reports/track2/fixture_graph/edges.tsv --drug-evidence reports/track2/fixture_graph/drug_evidence.tsv --seed DISEASE:SYNDIS00 --defect loss_of_function --defect-gene GENE:SYNG01 --out-dir reports/track2/fixture_ranking
