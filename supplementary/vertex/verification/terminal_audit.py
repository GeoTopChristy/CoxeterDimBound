"""Recompute the final rational comparisons from the supplied exact budgets.

This does not rerun the cusp census. With --check, require the stored table to
match the recomputed comparisons; otherwise regenerate terminal_results.json.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent


def calculate():
    budgets = json.loads((HERE / "budget_results.json").read_text(encoding="utf-8"))
    cusps = json.loads((HERE / "cusp_results.json").read_text(encoding="utf-8"))
    finite = {row["n"]: F(row["bound"]) for row in budgets["finite_DP"]}
    ideal = {row["n"]: F(row["bound"]) for row in cusps["cusp_table"]}
    assert set(finite) == set(ideal) == set(range(82, 90))
    rows = []
    for n in range(82, 90):
        threshold = F(n - 7, 8) + (F(3, 4 * n) if n % 2 else 0)
        bound = max(finite[n], ideal[n])
        gap = threshold - bound
        assert gap > 0, (n, threshold, bound)
        rows.append({"n": n, "T": str(threshold), "C": str(bound), "gap": str(gap)})

    # The comparison at dimension 81 uses the five-factor list (16,16,16,16,16).
    # Reconstruct its product h-polynomial and first g-moment directly.
    dimensions = (16,) * 5
    h = [1]
    for d in dimensions:
        new = [0] * (len(h) + d)
        for j, coefficient in enumerate(h):
            for k in range(d + 1):
                new[j + k] += coefficient
        h = new
    rank = sum(dimensions)
    vertices = sum(h)
    first = sum(F(j * (rank + 1 - 2 * j) * (h[j] - (h[j-1] if j else 0)), vertices)
                for j in range(rank // 2 + 1))
    affine = {row["d"]: max(F(value) for value in row["scores"].values())
              for row in budgets["affine"]}
    candidate = 5 * affine[16] - F(3, 2) * first
    threshold81 = F(81 - 7, 8) + F(3, 4 * 81)
    assert candidate == F(53148337, 5679428)
    assert threshold81 == F(250, 27)
    assert candidate > threshold81
    return {"terminal": rows, "dimension81_candidate": str(candidate),
            "dimension81_threshold": str(threshold81)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = calculate()
    output = HERE / "terminal_results.json"
    if args.check:
        assert json.loads(output.read_text(encoding="utf-8")) == result
    else:
        output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "passed", "terminal_rows": len(result["terminal"]),
                      "mode": "checked_stored_table" if args.check else "regenerated_table"}, indent=2))


if __name__ == "__main__":
    main()
