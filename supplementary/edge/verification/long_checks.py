"""Enumerate every critical position and every long type-IV overlap."""
from itertools import combinations, combinations_with_replacement, permutations
import json
import networkx as nx
import sympy as sp
from diagram_checks import ROOT,F,graph,subgraph,positivity,rational_matrix,adj,total_weight,weight

def ordinary(G):
    X=nx.Graph()
    X.add_nodes_from(range(G[0]))
    X.add_edges_from((a,b) for a,b,m in G[1])
    return X

def inv_diagonal(G,u):
    A=sp.Matrix(rational_matrix(G))
    return F(A[u,u]*A.inv()[u,u])

def deleted_inv(G,x,u):
    return inv_diagonal(subgraph(G,(x,)),u-(u>x))

def glue(G,K,identifications):
    index={v:u for u,v in identifications}
    n=G[0]
    for v in range(K[0]):
        if v not in index:
            index[v]=n
            n+=1
    edges={(a,b):m for a,b,m in G[1]}
    compatible=True
    for a,b,m in K[1]:
        p=tuple(sorted((index[a],index[b])))
        if p in edges and edges[p]!=m:
            compatible=False
        edges[p]=m
    return graph(n,[(a,b,m) for (a,b),m in edges.items()]),index,compatible

def run():
    source=json.loads((ROOT/'verification/high_cores.json').read_text())
    records=[]
    critical=[]
    for entry in source:
        G=graph(entry['n'],entry['edges'])
        X=ordinary(G)
        assert nx.diameter(X)<=8
        distances=dict(nx.all_pairs_shortest_path_length(X))
        for marks in combinations(range(G[0]),2):
            if not all(positivity(subgraph(G,(x,)))==2 for x in marks):
                continue
            d={u:min(distances[u][x] for x in marks) for u in X}
            assert max(d.values())<=7
            positions=[u for u in X if d[u]>=6]
            if not positions:
                continue
            assert G[0]>=9
            assert sum(d[u]==7 for u in positions)<=1
            record=dict(name=entry['name'],marks=marks,positions=positions,d=d,G=G)
            records.append(record)
            for u in positions:
                zeta=inv_diagonal(G,u)
                rho=max(deleted_inv(G,x,u) for x in marks)
                assert zeta<=0
                if d[u]==7: assert zeta==0
                if zeta<0: assert rho>=F(16,5)
                critical.append(dict(core=entry['name'],marks=marks,node=u,distance=d[u],zeta=str(zeta),rho=str(rho)))
    assert len(critical)==21 and len(records)==16
    candidates=[]
    for i,j in combinations_with_replacement(range(len(records)),2):
        L,M=records[i],records[j]
        for u in L['positions']:
            for v in M['positions']:
                if L['d'][u]+M['d'][v]>=13:
                    candidates.append((i,j,((u,v),)))
        for u1,u2 in combinations(L['positions'],2):
            if u2 not in adj(L['G'])[u1]: continue
            for v1,v2 in permutations(M['positions'],2):
                if v2 not in adj(M['G'])[v1]: continue
                if L['d'][u1]+M['d'][v1]<13 or L['d'][u2]+M['d'][v2]<13: continue
                ident=tuple(sorted(((u1,v1),(u2,v2))))
                inverse=tuple(sorted(((v1,u1),(v2,u2))))
                if i==j and ident>inverse: continue
                candidates.append((i,j,ident))
    outcomes=[]
    survivor_count=0
    statuses={}
    for i,j,ident in candidates:
        L,M=records[i],records[j]
        G,index,compatible=glue(L['G'],M['G'],ident)
        X=L['marks']
        Y=tuple(index[v] for v in M['marks'])
        row=dict(left=dict(core=L['name'],marks=X),right=dict(core=M['name'],marks=M['marks']),
                 identification=ident,n=G[0],edges=G[1])
        if not compatible:
            status='incompatible shared labels'
        else:
            terms=[]
            failures=[]
            for x in X:
                for y in Y:
                    deleted=tuple(a for a in X+Y if a not in (x,y))
                    kept=[a for a in range(G[0]) if a not in deleted]
                    T=subgraph(G,deleted)
                    pos=positivity(T)
                    if pos!=2: failures.append(dict(deleted=deleted,positivity=pos))
                    else:
                        xx,yy=kept.index(x),kept.index(y)
                        terms.append(dict(pair=[x,y],deleted=deleted,weight=str(weight(T,xx,yy))))
            if failures:
                status='nonpositive cross deletion'
                row['failures']=failures
            else:
                status='survivor'
                score=sum((F(t['weight']) for t in terms),F(0))
                assert score>=F(1,3)
                survivor_count+=1
                row.update(terms=terms,total=str(score))
        statuses[status]=statuses.get(status,0)+1
        row['status']=status
        outcomes.append(row)
    print('long candidates',len(candidates),'statuses',statuses,flush=True)
    assert len(candidates)==92 and survivor_count==13
    assert statuses=={'nonpositive cross deletion':76,'incompatible shared labels':3,'survivor':13}
    (ROOT/'tables/critical_positions.json').write_text(json.dumps(critical,indent=2)+'\n')
    (ROOT/'tables/long_candidates.json').write_text(json.dumps(outcomes,indent=2)+'\n')
    return dict(critical_positions=len(critical),marked_pairs=len(records),candidates=len(candidates),statuses=statuses)

if __name__=='__main__':
    print(run())
