"""Independent checks of the finite tree-certificate arithmetic and paths.

This reader does not import the generator. Matrix reconstruction is checked
separately by the independent replay program; this audit verifies the exact
2x2 certificate conditions and every recorded unchanged short path.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
import json

ROOT = Path(__file__).resolve().parent
counts = Counter()
minimum = None
for line in (ROOT / "tree/records.jsonl").open(encoding="utf-8"):
    row = json.loads(line)
    certificate = row["certificate"]
    counts[row["outcome"]] += 1
    if row["outcome"] == "negative":
        a, b = map(Fraction, certificate["matrix"][0])
        c, d = map(Fraction, certificate["matrix"][1])
        det = a*d-b*c
        assert b == c and a < 0 and det > 0
        assert det == Fraction(certificate["determinant"])
        assert Fraction(certificate["x"]) > 0 and Fraction(certificate["y"]) > 0
        minimum = det if minimum is None else min(minimum, det)
    elif row["outcome"] == "failed":
        if certificate["condition"] == "S+mark PD":
            assert Fraction(certificate["pivot"]) <= 0
        else:
            a, b = map(Fraction, certificate["matrix"][0])
            c, d = map(Fraction, certificate["matrix"][1])
            assert b == c
            assert a < 0 or d < 0 or a*d-b*c < 0
    elif row["outcome"] == "short":
        P, distance, Q = row["P"], row["ell"], row["Q"]
        last = P+distance+Q
        a, b = P, P+distance
        edges = {frozenset((i, i+1)) for i in range(last)}
        leaf = last+1
        for endpoint, kind in zip((0, last), row["endtype"]):
            for _ in range({"A": 0, "B": 1, "D": 2}[kind]):
                edges.add(frozenset((endpoint, leaf)))
                leaf += 1
        for name in ("c", "d"):
            for site, label in row[f"support_{name}"]:
                assert label in (3, 4)
                edges.add(frozenset((name, site)))
        if row["variant_c"][0]:
            edges.add(frozenset(("c", a)))
        if row["variant_d"][0]:
            edges.add(frozenset(("d", b)))
        changing = set()
        for variable, left, right, base in ((P, 2, a-2, 9),
                                          (distance, a+2, b-2, 10),
                                          (Q, b+2, last-2, 9)):
            if variable == base:
                changing.update(frozenset((i, i+1)) for i in range(left, right))
        path = certificate["path"]
        pair = certificate["pair"]
        assert path[0] == pair[0] and path[-1] == pair[1]
        assert len(path)-1 == certificate["distance"] <= 7
        assert len(path) == len(set(path))
        forbidden = {a, b, "c", "d"} - set(pair)
        assert not set(path) & forbidden
        for u, v in zip(path, path[1:]):
            edge = frozenset((u, v))
            assert edge in edges and edge not in changing
        assert certificate["avoids_all_variable_segments"] is True
    else:
        raise AssertionError(row["outcome"])

assert counts == {"negative": 39389, "short": 11167, "failed": 28557}
assert minimum == Fraction(1, 32400)
result = {
    "records_checked": sum(counts.values()),
    "outcomes": dict(counts),
    "all_negative_certificate_matrices_exactly_negative_definite": True,
    "all_witness_coordinates_positive": True,
    "all_short_paths_are_valid_adjacent_angle_paths": True,
    "all_short_paths_avoid_all_variable_internal_segments": True,
    "all_failed_certificates_violate_the_stated_2x2_or_scalar_condition": True,
    "minimum_negative_witness_determinant": str(minimum),
    "scope_limit": "Does not independently reconstruct Schur matrices; the independent replay handles reconstruction. The support classification and arbitrary-length theorem remain mathematical inputs.",
}
(ROOT / "vertex_tree_certificate_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
