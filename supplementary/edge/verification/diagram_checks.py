"""Exact finite checks for the single weight used in the dimension-162 paper."""
from fractions import Fraction as F
from collections import Counter
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json
import networkx as nx

ROOT = Path(__file__).absolute().parents[1]
ALPHA, DELTA = F(73, 12), F(23, 4)

def graph(n, edges):
    return n, tuple(sorted((min(a,b), max(a,b), m) for a,b,m in edges))

def subgraph(G, deleted):
    n, edges = G
    kept = [i for i in range(n) if i not in deleted]
    index = {v:i for i,v in enumerate(kept)}
    return graph(len(kept), [(index[a],index[b],m) for a,b,m in edges
                             if a in index and b in index])

def adj(G, deleted=()):
    n, edges = G
    out = {i:{} for i in range(n) if i not in deleted}
    for a,b,m in edges:
        if a in out and b in out:
            out[a][b] = m
            out[b][a] = m
    return out

def rational_matrix(G):
    """A positive diagonal congruence of 2G; None implies not PSD."""
    n, edges = G
    A = adj(G)
    parity = {}
    for root in range(n):
        if root in parity:
            continue
        parity[root] = 0
        stack = [root]
        while stack:
            u = stack.pop()
            for v,m in A[u].items():
                wanted = parity[u] ^ (m==4)
                if v in parity:
                    if parity[v] != wanted:
                        return None
                else:
                    parity[v] = wanted
                    stack.append(v)
    out = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        out[i][i] = F(2 << parity[i])
    for a,b,m in edges:
        out[a][b] = out[b][a] = F(-2 if m==4 or parity[a] else -1)
    return out

@lru_cache(maxsize=300000)
def positivity(G):
    """Return 2 (PD), 1 (singular PSD), or 0 (not PSD), exactly."""
    A = rational_matrix(G)
    if A is None:
        return 0
    singular = False
    while A:
        if any(A[i][i] < 0 for i in range(len(A))):
            return 0
        if A[0][0] == 0:
            if any(A[0][1:]):
                return 0
            singular = True
            A = [r[1:] for r in A[1:]]
        else:
            p = A[0][0]
            A = [[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,len(A))]
                 for i in range(1,len(A))]
    return 1 if singular else 2

@lru_cache(maxsize=300000)
def valid(G, conditions):
    return all(positivity(subgraph(G, D)) >= level for D,level in conditions)

def path_data(G, a, b, deleted=()):
    A = adj(G, deleted)
    stack = [(a, [a])]
    seen = {a}
    while stack:
        u, p = stack.pop()
        if u == b:
            special = any(len(A[z]) >= 3 for z in p) or any(
                A[x][y] >= 4 for x,y in zip(p,p[1:]))
            return p, special
        for v in A[u]:
            if v not in seen:
                seen.add(v)
                stack.append((v,p+[v]))
    return None, False

def base_weight(d, special):
    if d <= 6: return F(1)
    if d == 7: return F(1) if special else F(1,2)
    if d <= 12: return F(1,3)
    if d == 13: return F(1,3) if special else F(1,6)
    if d == 14: return F(1,6) if special else F(1,12)
    return F(0)

def has_arm(A, u, prev, steps):
    for _ in range(steps):
        onward = [z for z in A[u] if z != prev]
        if len(onward) != 1 or A[u][onward[0]] != 3:
            return False
        prev, u = u, onward[0]
    return True

def weight(G,a,b):
    p,special = path_data(G,a,b)
    if p is None: return F(0)
    d = len(p)-1
    if not special and 3 <= d <= 7:
        A = adj(G)
        if has_arm(A,p[0],p[1],9-d) and has_arm(A,p[-1],p[-2],9-d):
            return F(1,2) if d<=5 else F(1,3)
    return base_weight(d,special)

def total_weight(G):
    return sum((weight(G,a,b) for a,b in combinations(range(G[0]),2)), F(0))

@lru_cache(maxsize=200000)
def possible_arm(G, u, prev, remaining, core_n, deleted, conditions):
    if remaining == 0:
        return True
    A = adj(G,deleted)
    onward = [z for z in A[u] if z != prev]
    if len(onward)>1:
        return False
    if onward:
        z = onward[0]
        return A[u][z]==3 and possible_arm(G,z,u,remaining-1,core_n,deleted,conditions)
    n,edges = G
    candidates = [n]
    if u >= core_n:
        component = {u}
        stack = [u]
        while stack:
            for z in A[stack.pop()]:
                if z not in component:
                    component.add(z)
                    stack.append(z)
        candidates += [z for z in A if z not in component
                       and (remaining==1 or len(A[z])<=1)]
    for z in candidates:
        nxt = graph(n+(z==n),edges+((min(u,z),max(u,z),3),))
        if valid(nxt,conditions) and possible_arm(nxt,z,u,remaining-1,core_n,deleted,conditions):
            return True
    return False

def safe_weight(G, a,b,deleted,conditions):
    p,special = path_data(G,a,b,deleted)
    if p is None: return F(0), 'disconnected'
    d = len(p)-1
    if not special and 3<=d<=7:
        left = possible_arm(G,p[0],p[1],9-d,G[0],deleted,conditions)
        right = possible_arm(G,p[-1],p[-2],9-d,G[0],deleted,conditions)
        if left and right:
            return (F(1,2) if d<=5 else F(1,3)), 'relaxed'
        return base_weight(d,special), 'buffer excluded'
    return base_weight(d,special), 'path fixed'

def generate_cores():
    """Read one labelled representative of each classification reduction.

    Every label at least four is reduced to four, then labelled isomorphic
    diagrams are identified. All permitted marked-node choices are retained
    separately by check_local; markings are not quotiented by automorphisms.
    The input includes its classification sources. A reduced diagram need
    not remain indefinite, so definiteness is imposed only on deletions.
    """
    entries = json.loads((ROOT/'verification/reduced_cores.json').read_text())
    assert Counter(item['n'] for item in entries)=={3:5,4:17,5:13,6:12,7:3,8:4,9:4,10:3}
    assert len({item['name'] for item in entries}) == len(entries) == 61
    seen = []
    cores = []
    edge_match = nx.algorithms.isomorphism.categorical_edge_match('label',3)
    for item in entries:
        n = item['n']
        assert item['classification_sources'], item['name']
        assert all(0 <= a < b < n and m in (3,4) for a,b,m in item['edges'])
        C = graph(n,item['edges'])
        X = nx.Graph()
        X.add_nodes_from(range(n))
        X.add_edges_from((a,b,{'label':m}) for a,b,m in C[1])
        assert nx.is_connected(X)
        assert not any(nx.is_isomorphic(X,Y,edge_match=edge_match) for Y in seen)
        seen.append(X)
        assert all(positivity(subgraph(C,(v,))) >= 1 for v in range(n))
        assert all(positivity(subgraph(C,D)) == 2 for D in combinations(range(n),2))
        cores.append((item['name'],C))
    return cores

def check_local(cores):
    rows = []
    counts = dict(I=0,II=0,III=0)
    minima = dict(I=F(100),II=F(100),III=F(100))
    for name,G in cores:
        n = G[0]
        for marks in combinations(range(n),4):
            cond = tuple([((a,),1) for a in marks]+[(D,2) for D in combinations(marks,2)])
            if valid(G,cond):
                pairs = [(a,b,tuple(x for x in marks if x not in (a,b))) for a,b in combinations(marks,2)]
                record_local(rows,counts,minima,name,G,'I',marks,cond,pairs)
        for marks in combinations(range(n),3):
            cond = tuple(((a,),2) for a in marks)
            if valid(G,cond):
                pairs = [(a,b,tuple(x for x in marks if x not in (a,b))) for a,b in combinations(marks,2)]
                record_local(rows,counts,minima,name,G,'II',marks,cond,pairs)
            for a in marks:
                b,c = [z for z in marks if z!=a]
                cond = tuple(sorted((((a,),1),((b,),2),((c,),2))))
                if valid(G,cond):
                    record_local(rows,counts,minima,name,G,'III',(a,b,c),cond,[(a,b,(c,)),(a,c,(b,))])
        print('core',name,'complete',flush=True)
    expected = dict(I=F(3),II=F(2),III=F(1))
    assert minima == expected, minima
    assert counts==dict(I=1781,II=728,III=2761), counts
    (ROOT/'tables/local_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    return dict(counts=counts,minima={k:str(v) for k,v in minima.items()})

def record_local(rows,counts,minima,name,G,kind,marks,cond,pairs):
    terms=[]
    score=F(0)
    for a,b,D in pairs:
        value,reason = safe_weight(G,a,b,D,cond)
        score+=value
        terms.append(dict(pair=[a,b],deleted=D,lower_bound=str(value),reason=reason))
    counts[kind]+=1
    minima[kind]=min(minima[kind],score)
    if score < {'I':3,'II':2,'III':1}[kind]:
        raise AssertionError((name,kind,marks,score,terms))
    rows.append(dict(core=name,type=kind,marks=marks,conditions=cond,terms=terms,total=str(score)))

if __name__=='__main__':
    cores=generate_cores()
    (ROOT/'tables/minimal_core_catalog.json').write_text(json.dumps([
        dict(name=name,n=G[0],edges=G[1]) for name,G in cores],indent=2)+'\n')
    print(check_local(cores))
