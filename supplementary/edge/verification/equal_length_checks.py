"""Exact local verification of the edge-length multiplicity calculation.

The infinite-family proof is the counting formula, not the finite check.
For exceptional types, reconstruct a simplex using perpendicular tree edges,
verify its facet-angle Gram matrix, and list every squared length and pair.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
import sympy as sp

ROOT = Path(__file__).absolute().parents[1]
OUT = ROOT / 'tables'


def tree_simplex(name, edges):
    nodes = ['O']
    pos = {'O': sp.zeros(len(edges), 1)}
    for k, (parent, child, weight, label) in enumerate(edges):
        assert parent in pos and child not in pos
        pos[child] = pos[parent].copy()
        pos[child][k] += 1
        nodes.append(child)
    weights = [sp.Rational(e[2]) for e in edges]
    metric = sp.diag(*weights)
    basis = sp.Matrix.hstack(*(pos[x] for x in nodes[1:]))
    covectors = basis.inv()
    gradients = sp.Matrix.vstack(-sp.ones(1, len(edges))*covectors, covectors)
    gram = gradients * metric.inv() * gradients.T
    expected = {frozenset((u,v)):m for u,v,w,m in edges}
    for i,j in combinations(range(len(nodes)), 2):
        m = expected.get(frozenset((nodes[i], nodes[j])), 2)
        assert gram[i,j] <= 0
        assert sp.simplify(gram[i,j]**2/(gram[i,i]*gram[j,j]) - sp.cos(sp.pi/m)**2) == 0, (name,i,j,m)
    pairs = {}
    for u,v in combinations(nodes, 2):
        diff = pos[u]-pos[v]
        square = (diff.T*metric*diff)[0]
        pairs.setdefault(square, []).append([u,v])
    assert sum(map(len, pairs.values())) == len(nodes)*(len(nodes)-1)//2
    assert max(map(len, pairs.values())) <= sp.Rational(6,5)*len(nodes)
    return {'type':name, 'dimension':len(edges),
            'orthogonal_tree_edges':edges,
            'squared_normal_norms':list(map(str,gram.diagonal())),
            'squared_lengths':[{'squared_length':str(t),'multiplicity':len(pairs[t]),'pairs':pairs[t]} for t in sorted(pairs)],
            'maximum_multiplicity':max(map(len,pairs.values()))}


def star_edges(arms):
    edges = []
    for letter, ws in zip('ABC', arms):
        prev = 'O'
        for j,w in enumerate(ws,1):
            node = letter+str(j)
            edges.append((prev,node,w,3))
            prev = node
    return edges


def classical(kind, n):
    p = lambda k: sp.Matrix([int(j<k) for j in range(n)])
    if kind == 'C':
        vertices = [p(i) for i in range(n+1)]
    elif kind == 'B':
        vertices = [p(0),2*p(1)]+[p(i) for i in range(2,n+1)]
    else:
        end = p(n-1)
        axis = p(n)-p(n-1)
        vertices = [p(0),2*p(1)]+[p(i) for i in range(2,n-1)]+[end-axis,end+axis]
    counts = Counter(int((v-u).dot(v-u)) for u,v in combinations(vertices,2))
    for k in range(1,max(n,4)+1):
        if kind == 'C': formula=max(n+1-k,0)
        elif kind == 'B': formula=max(n-1-k,0)+2*(2<=k<=n)+(k==4)
        else: formula=max(n-3-k,0)+4*(2<=k<=n-2)+2*(k==4)+4*(k==n)
        assert counts[k] == formula, (kind,n,k,counts[k],formula)
    expected = n if kind=='C' else n-1 if kind=='B' else {4:6,5:4,6:6}.get(n,n-1)
    assert max(counts.values()) == expected
    assert max(counts.values()) <= sp.Rational(6,5)*(n+1)
    return {str(k):v for k,v in sorted(counts.items())}


def main():
    OUT.mkdir(exist_ok=True)
    data = [tree_simplex('E6',star_edges([[1,3],[1,3],[1,3]])),
            tree_simplex('E7',star_edges([[3],[2,4,12],[2,4,12]])),
            tree_simplex('E8',star_edges([[20],[15,45],[12,18,30,60,180]])),
            tree_simplex('F4',[('O','A',3,3),('A','B',1,4),('B','C',2,3),('C','D',6,3)]),
            tree_simplex('G2',[('O','A',1,6),('A','B',3,3)])]
    small = {}
    for kind, minimum in [('C',2),('B',3),('D',4)]:
        for n in range(minimum,21):
            counts=classical(kind,n)
            if n <= 8: small[kind+str(n)]=counts
    # Independent comparison with the original Gram-matrix reconstruction.
    import sys
    sys.path.insert(0,str(ROOT / 'verification'))
    from metric_checks import metric_data
    for item in data:
        name=item['type']; n=item['dimension']
        ss,mm,mx=metric_data(name[0],n)
        counts=[x['multiplicity'] for x in item['squared_lengths']]
        assert mx == max(counts)
        assert sp.Rational(sum(x*x for x in counts),(n+1)**2) == sp.Rational(mm.numerator,mm.denominator)
    result={'arithmetic':'exact','exceptional_gram_checks':'passed','classical_formula_checks_through_dimension':20,
            'original_gram_reconstruction_cross_check':'passed',
            'exceptional_types':data,'small_classical_counts':small}
    (OUT/'equal_length_counts.json').write_text(json.dumps(result,indent=2)+'\n')
    for x in data: print(x['type'],x['maximum_multiplicity'],len(x['squared_lengths']))
    print('All exact checks passed.')


if __name__ == '__main__':
    main()
