"""Verify the code-only exact-arithmetic package.

Default: validate the manifest and supplied finite summaries, and recompute
the final rational comparisons. --quick adds the short arithmetic audits.
--replay-tree independently reconstructs and replays the supplied certificates.
--regenerate reruns the finite enumerations and regenerates all check results.
No manuscript, network connection, or TeX installation is required.
"""
from fractions import Fraction
from pathlib import Path
import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "supplementary"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(path, *args):
    print("Checking " + path.relative_to(ROOT).as_posix(), flush=True)
    logs = ROOT / "verification/run_logs"
    logs.mkdir(parents=True, exist_ok=True)
    name = path.relative_to(SUPP).as_posix().replace("/", "__") + ".log"
    with (logs / name).open("wb") as stream:
        subprocess.run([sys.executable, str(path), *map(str, args)], cwd=path.parent,
                       stdout=stream, stderr=subprocess.STDOUT, check=True)


def manifest_files():
    allowed = {".py", ".json", ".jsonl", ".csv", ".txt", ".md"}
    for path in sorted(ROOT.rglob("*")):
        relative = path.relative_to(ROOT)
        if not path.is_file() or any(part in (".git", "__pycache__", "run_logs", ".venv", "venv")
                                     for part in relative.parts):
            continue
        if relative.as_posix() == "SHA256SUMS.json":
            continue
        if path.suffix in allowed or relative.as_posix() in (".gitignore", ".gitattributes"):
            yield path


def refresh_manifest():
    manifest = {p.relative_to(ROOT).as_posix(): sha256(p) for p in manifest_files()}
    (ROOT / "SHA256SUMS.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    parser.add_argument("--replay-tree", action="store_true")
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    manifest = json.loads((ROOT / "SHA256SUMS.json").read_text(encoding="utf-8"))
    for name, expected in manifest.items():
        target = (ROOT / name).resolve()
        assert ROOT in target.parents and target.is_file(), name
        assert sha256(target) == expected, "Integrity mismatch: " + name

    if args.regenerate:
        run(SUPP / "edge/verification/verify_all.py")
        for name in ("equal_length_checks.py", "link_sum_checks.py", "cusp_defect_checks.py",
                     "coordinate_proof_checks.py"):
            run(SUPP / "edge/verification" / name)
        run(SUPP / "vertex/verification/verify_available.py")
        run(SUPP / "tree/tree_audit.py")
    else:
        run(SUPP / "vertex/verification/verify_available.py", "--cached")

    if args.regenerate or args.quick:
        run(SUPP / "shared_moment_audit.py")
        run(SUPP / "edge/verification/square_link_geometry.py")
        run(SUPP / "edge/verification/chapter4_readability_checks.py", "--report",
            SUPP / "edge/logs/edge_examples_checks.json")

    if args.regenerate or args.replay_tree:
        run(SUPP / "tree/verify_certificates.py")
    if args.regenerate or args.replay_tree or args.quick:
        run(SUPP / "vertex_tree_certificate_audit.py")

    terminal = SUPP / "vertex/verification/terminal_audit.py"
    if args.regenerate:
        run(terminal)
    else:
        run(terminal, "--check")

    edge = json.loads((SUPP / "edge/logs/verification_summary.json").read_text(encoding="utf-8"))
    assert edge["status"] == "passed"
    assert edge["results"]["local"]["counts"] == {"I": 1781, "II": 728, "III": 2761}
    assert edge["results"]["local"]["minima"] == {"I": "3", "II": "2", "III": "1"}
    assert edge["results"]["marking_counts"]["candidate_entries"] == 61
    assert edge["results"]["type_IV"]["statuses"] == {
        "nonpositive cross deletion": 76, "incompatible shared labels": 3, "survivor": 13}
    generation = json.loads((SUPP / "tree/summary.json").read_text(encoding="utf-8"))
    replay = json.loads((SUPP / "tree/replay_summary.json").read_text(encoding="utf-8"))
    expected = {"negative": 39389, "short": 11167, "failed": 28557}
    assert generation["counts"] == replay["counts"] == expected
    assert generation["initial_chains"] == 450 and generation["all_records_certified"]
    assert generation["unresolved_count"] == 0
    assert json.loads((SUPP / "tree/unresolved.json").read_text(encoding="utf-8")) == []
    assert replay["all_certificates_replayed"]
    assert replay["distinct_records"] == 79113 and replay["distinct_initial_chains"] == 450
    assert Fraction(replay["minimum_determinant"]) == Fraction(1, 32400)
    assert sha256(SUPP / "tree/records.jsonl") == replay["record_sha256"]

    report = {
        "status": "passed",
        "mode": "regenerated_all_computations" if args.regenerate else
                "validated_supplied_results_and_quick_arithmetic" if args.quick else
                "validated_supplied_exact_results",
        "tree_replayed_in_this_run": bool(args.regenerate or args.replay_tree),
        "tree_certificate_arithmetic_checked_in_this_run": bool(args.regenerate or args.replay_tree or args.quick),
        "terminal_comparisons_recomputed_in_this_run": True,
        "tree_records": 79113, "tree_counts": replay["counts"],
        "classification_reductions": 61, "three_face_markings": 5270,
        "planar_cusp_lists": 2770038,
        "python": platform.python_version(), "platform": platform.platform(),
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "scope": "Exact finite arithmetic and supplied-result validation; geometric completeness and arbitrary-rank arguments are mathematical inputs, not formalized here."
    }
    (ROOT / "verification/verification_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    # Some result files include timing metadata; record the package after this run.
    refresh_manifest()
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
