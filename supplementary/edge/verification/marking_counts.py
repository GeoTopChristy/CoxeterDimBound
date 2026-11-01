"""Reconstruct the marked-record counts in the local lower-bound proof using exact deletions.

Run: python verification/marking_counts.py
This counts markings before the buffer search. It also compares every marking
with the records in tables/local_checks.json; no floating-point tests are used.
"""
from collections import Counter
from itertools import combinations
from math import comb
import csv
import json

from diagram_checks import ROOT, generate_cores, positivity, subgraph, valid


def choose(n, k):
    return comb(n, k) if n >= k else 0


def check_counts(candidates=None):
    if candidates is None:
        candidates = generate_cores()
    rows, expected_records = [], set()
    for name, G in candidates:
        r = G[0]
        status = [positivity(subgraph(G, (v,))) for v in range(r)]
        assert all(s >= 1 for s in status), name
        assert all(positivity(subgraph(G, D)) == 2
                   for D in combinations(range(r), 2)), name
        positive = [v for v, s in enumerate(status) if s == 2]
        singular = [v for v, s in enumerate(status) if s == 1]
        p = len(positive)
        formula = dict(I=choose(r, 4), II=choose(p, 3),
                       III=(r-2)*choose(p, 2))
        enumerated = Counter()

        # Enumerate the deletion conditions separately from the formulas.
        for marks in combinations(range(r), 4):
            conditions = tuple([((a,), 1) for a in marks]
                               + [(D, 2) for D in combinations(marks, 2)])
            if valid(G, conditions):
                expected_records.add((name, 'I', marks))
                enumerated['I'] += 1
        for marks in combinations(range(r), 3):
            conditions = tuple(((a,), 2) for a in marks)
            if valid(G, conditions):
                expected_records.add((name, 'II', marks))
                enumerated['II'] += 1
            for a in marks:
                b, c = [v for v in marks if v != a]
                conditions = tuple(sorted((((a,), 1), ((b,), 2), ((c,), 2))))
                if valid(G, conditions):
                    expected_records.add((name, 'III', (a, b, c)))
                    enumerated['III'] += 1
        assert all(enumerated[k] == formula[k] for k in formula), name
        rows.append(dict(name=name, r=r, edges=G[1], p=p,
                         positive_deletions=positive,
                         singular_deletions=singular,
                         counts=formula,
                         type_III_singular_a=(r-p)*choose(p, 2)))

    # Compare all 5270 identifiers, not just the three totals.
    records = json.loads((ROOT/'tables/local_checks.json').read_text())
    stored = [(rec['core'], rec['type'], tuple(rec['marks'])) for rec in records]
    assert len(stored) == len(set(stored))
    assert set(stored) == expected_records

    summary = []
    for r in range(3, 11):
        group = [row for row in rows if row['r'] == r]
        summary.append(dict(r=r, candidates=len(group),
                            p_distribution=dict(sorted(Counter(row['p'] for row in group).items())),
                            counts={k: sum(row['counts'][k] for row in group)
                                    for k in ('I', 'II', 'III')}))
    totals = {k: sum(row['counts'][k] for row in rows) for k in ('I', 'II', 'III')}
    singular_a = sum(row['type_III_singular_a'] for row in rows)
    assert len(rows) == 61
    assert totals == dict(I=1781, II=728, III=2761)
    assert totals['III'] == 3*totals['II'] + singular_a
    assert singular_a == 577
    example = next(row for row in rows if row['name'] == 'Q9-E')
    assert example['positive_deletions'] == [0, 1, 2, 3, 4, 5, 6, 8]
    assert example['singular_deletions'] == [7]
    assert example['counts'] == dict(I=126, II=56, III=196)

    result = dict(method='Exact rational congruence and deletion tests; no quotient by marked automorphisms.',
                  candidate_entries=len(rows), summary=summary, totals=totals,
                  type_III_singular_a=singular_a,
                  comparison_with_local_checks='All 5270 marked identifiers agree.',
                  candidates=rows)
    (ROOT/'tables/marking_counts.json').write_text(json.dumps(result, indent=2)+'\n')
    with (ROOT/'tables/marking_counts.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['diagram', 'r', 'p', 'positive_deletions', 'singular_deletions',
                         'I', 'II', 'III', 'III_singular_a'])
        for row in rows:
            writer.writerow([row['name'], row['r'], row['p'],
                             ' '.join(map(str, row['positive_deletions'])),
                             ' '.join(map(str, row['singular_deletions'])),
                             *[row['counts'][k] for k in ('I', 'II', 'III')],
                             row['type_III_singular_a']])
    return dict(candidate_entries=len(rows), counts=totals, type_III_singular_a=singular_a,
                comparison_with_local_checks=result['comparison_with_local_checks'])


if __name__ == '__main__':
    print(json.dumps(check_counts(), indent=2))
