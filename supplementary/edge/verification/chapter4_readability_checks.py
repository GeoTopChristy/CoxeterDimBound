"""Exact checks for the examples and algebra in the rewritten Chapter 4.

This checks the exposition changes, not the full Coxeter classification.
Run with Python and SymPy; use --report PATH to save the JSON report.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import json

import sympy as sp


def check():
    C, delta = F(73, 12), F(23, 4)
    # Full-product sum and the nonuniformly selected subset.
    losses = [2*C-1, 2*C-1, 2*C]
    vertices = list(product(range(3), range(2)))
    mean = sum(losses)/3
    assert Counter(a for a,b in vertices) == {0:2, 1:2, 2:2}
    assert sum(losses[a] for a,b in vertices) == 6*mean
    selected = [(0,0), (0,1), (1,0)]
    assert sum(losses[a] for a,b in selected) == 6*C-3 < 3*mean == 6*C-2
    # Match squared lengths in two identical right isosceles triangles.
    pts = [(0,0), (1,0), (0,1)]
    edges = {(i,j): sum((a-b)**2 for a,b in zip(pts[i],pts[j]))
             for i,j in combinations(range(3), 2)}
    counts = Counter(edges.values())
    assert counts == {1:2, 2:1}
    matches = sum(a == b for a in edges.values() for b in edges.values())
    assert matches == sum(x*x for x in counts.values()) == 5
    assert F(matches, 9) == F(5,9)
    S = [sum(x*x for x in Counter(t for ij,t in edges.items() if v in ij).values())
         for v in range(3)]
    assert S == [4,2,2]
    assert 2*C-1 < C*4/2 == (2*C-1)+1
    assert 2*C-2 < 4*C*1/2
    # Rectangles of a product of three intervals, counted by their corners.
    cube = list(product(range(2), repeat=3))
    faces = []
    for factors in combinations(range(3), 2):
        fixed = next(i for i in range(3) if i not in factors)
        for value in range(2):
            faces.append([v for v in cube if v[fixed] == value])
    assert len(faces) == 6 and all(len(face) == 4 for face in faces)
    assert all(sum(v in face for face in faces) == 3 for v in cube)
    assert len(cube)*3 == len(faces)*4 == 24
    # Small-pair alternating-color bounds: exhaust the finite arrays.
    color_maxima = {}
    for rows,cols in [(2,3),(3,3)]:
        maximum = 0
        for coloring in product(range(2), repeat=rows*cols):
            count = 0
            for a,b in combinations(range(rows),2):
                for c,d in combinations(range(cols),2):
                    p,q,r,s = [coloring[x*cols+y] for x,y in [(a,c),(a,d),(b,c),(b,d)]]
                    count += p == s and q == r and p != q
            maximum = max(maximum, count)
        color_maxima[f'{rows}x{cols}'] = maximum
    assert color_maxima == {'2x3':2, '3x3':4}
    corrected_lower = {n: n*(C-F(n-1,2))-int(n==2) for n in range(1,9)}
    assert all(x >= C for x in corrected_lower.values())
    assert min(corrected_lower[n]/n for n in corrected_lower) == F(31,12)
    small_pairs = 0
    for nj in range(3,9):
        assert corrected_lower[nj] >= F(61,4)
        for ni in range(1,nj+1):
            assert corrected_lower[ni]/ni >= C-F(nj-1,2)
            assert ni*corrected_lower[nj]+nj*corrected_lower[ni] >= F(61,2)*ni > 24*C*ni/5
            small_pairs += 1
    for (ni,nj), normalized in [((1,1),F(1,4)), ((1,2),F(1,3)), ((2,2),F(4,9))]:
        assert ni*corrected_lower[nj]+nj*corrected_lower[ni] >= 8*C*normalized
        small_pairs += 1
    assert F(17,2)+9*F(31,12) == F(127,4) > 4*C
    gaps = [delta+2-C*(1+F(2,21)), 2*delta+2-C*(2+F(2,21)),
            delta+3-C*(1+F(4,21)), 2*delta+3-C*(2+F(8,21))]
    assert min(gaps) == F(1,63)
    # Formal local-to-global identities, without choosing a dimension.
    n,c,d,B,f,fu,k = sp.symbols('n C delta B f f_u correction')
    lam = 2/(n-2)*(c-d/(n-1))
    assert sp.cancel((c*(n-1)-d)/((n-1)*(n-2)/2)-lam) == 0
    I = 8*c*B/(n-1)+f*(c-d)+f*k
    J = 4*c*B/(n-1)+fu*(c-d)-fu*k
    combined = 12*c*B/(n-1)+(f+fu)*(c-d)+(f-fu)*k
    assert sp.expand(I+J-combined) == 0
    local = 12*c/(n-1)+12*(c-d)/((n-1)*(n-2))
    assert sp.cancel(local-6*lam) == 0
    assert sp.simplify(local*3/2-9*lam) == 0
    # Endpoint multiplicities: I counts ideal endpoints; J adds the mixed edge.
    multiplicities = {}
    for endpoints in [(0,0), (0,1), (1,1)]:
        ideal_ends = sum(endpoints)
        I_count = ideal_ends
        J_count = int(ideal_ends == 1)
        multiplicities[str(endpoints)] = I_count+J_count
        assert I_count+J_count == (2 if ideal_ends else 0)
    return {'status':'passed', 'arithmetic':'exact rational and symbolic',
            'examples':['triangle times interval', 'selected subset below the full mean',
                        'right-isosceles matching counts and vertex counts',
                        'equilateral correction', 'cube rectangle incidences'],
            'alternating_color_maxima':color_maxima,
            'small_dimension_pairs_checked':small_pairs,
            'stable_pointwise_minimum_at_21':str(min(gaps)),
            'symbolic_identities':['normalized edge coefficient', 'sum of I and J with correction remainder',
                                   'local coefficient 6 lambda_n', 'global coefficient 9 lambda_n'],
            'edge_endpoint_multiplicities':multiplicities,
            'scope':'New examples and rewritten algebra only; full classification suite not rerun.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    report = json.dumps(check(), indent=2)+'\n'
    if args.report:
        args.report.write_text(report, encoding='utf-8')
    print(report)
