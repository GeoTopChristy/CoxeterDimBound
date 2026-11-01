"""Independent exact reconstruction of21tails and92quadrilateral overlaps.

Only the standard diagram edge catalogue is shared with the existing162
paper; its numerical evaluators are not imported.
"""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement,permutations
from functools import lru_cache
from collections import deque,Counter
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent
ROOT=OUT

def graph(n,edges):
    return n,tuple(sorted((min(a,b),max(a,b),m) for a,b,m in edges))

def adj(G):
    n,edges=G;A=[{} for _ in range(n)]
    for a,b,m in edges:A[a][b]=A[b][a]=m
    return A

def deletion(G,removed):
    n,edges=G;kept=[i for i in range(n) if i not in removed];ix={u:i for i,u in enumerate(kept)}
    return graph(len(kept),[(ix[a],ix[b],m) for a,b,m in edges if a in ix and b in ix])

def matrix(G):
    n,edges=G;A=adj(G);par={}
    for start in range(n):
        if start in par:continue
        par[start]=0;stack=[start]
        while stack:
            u=stack.pop()
            for v,m in A[u].items():
                want=par[u]^(m==4)
                if v in par and par[v]!=want:return None
                if v not in par:par[v]=want;stack.append(v)
    M=[[F(0) for _ in range(n)] for _ in range(n)]
    for u in range(n):M[u][u]=F(2<<par[u])
    for a,b,m in edges:M[a][b]=M[b][a]=F(-2 if m==4 or par[a] else -1)
    return M

@lru_cache(maxsize=None)
def positivity(G):
    A=matrix(G)
    if A is None:return 0
    singular=False
    while A:
        if A[0][0]<0:return 0
        if A[0][0]==0:
            if any(A[0][1:]):return 0
            singular=True;A=[row[1:] for row in A[1:]];continue
        p=A[0][0]
        A=[[A[i][j]-A[i][0]*A[0][j]/p for j in range(1,len(A))] for i in range(1,len(A))]
    return 1 if singular else 2

def inverse_diag(G,u):
    A=matrix(G);scale=A[u][u];n=len(A)
    B=[row+[F(i==u)] for i,row in enumerate(A)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if B[i][k])
        B[k],B[pivot]=B[pivot],B[k]
        p=B[k][k];B[k]=[x/p for x in B[k]]
        for i in range(n):
            if i==k:continue
            t=B[i][k]
            if t:B[i]=[x-t*y for x,y in zip(B[i],B[k])]
    return scale*B[u][-1]

def distances(G,start):
    A=adj(G);d={start:0};q=deque([start])
    while q:
        u=q.popleft()
        for v in A[u]:
            if v not in d:d[v]=d[u]+1;q.append(v)
    return d

def glue(L,M,pairs):
    ix={v:u for u,v in pairs};n=L[0]
    for v in range(M[0]):
        if v not in ix:ix[v]=n;n+=1
    edges={(a,b):m for a,b,m in L[1]};okay=True
    for a,b,m in M[1]:
        p=tuple(sorted((ix[a],ix[b])))
        if p in edges and edges[p]!=m:okay=False
        edges[p]=m
    return graph(n,[(a,b,m) for (a,b),m in edges.items()]),ix,okay

def preliminary(G,a,b):
    A=adj(G);parents={a:None};q=deque([a])
    while q:
        u=q.popleft()
        for v in A[u]:
            if v not in parents:parents[v]=u;q.append(v)
    if b not in parents:return F(0)
    p=[b]
    while p[-1]!=a:p.append(parents[p[-1]])
    p.reverse();d=len(p)-1
    special=any(len(A[u])>=3 for u in p) or any(A[u][v]>=4 for u,v in zip(p,p[1:]))
    return F(1) if d<=12 or special and d<=14 else F(1,2) if d==13 else F(1,4) if d==14 else F(0)

def run():
    source=json.loads((OUT/'high_cores.json').read_text())
    records=[];tails=[]
    for item in source:
        G=graph(item['n'],item['edges']);ds=[distances(G,u) for u in range(G[0])]
        assert max(max(d.values()) for d in ds)<=8
        for marks in combinations(range(G[0]),2):
            if any(positivity(deletion(G,{x}))!=2 for x in marks):continue
            near={u:min(ds[u][x] for x in marks) for u in range(G[0])}
            assert max(near.values())<=7
            pos=[u for u in near if near[u]>=6]
            if not pos:continue
            assert sum(near[u]==7 for u in pos)<=1
            if len(pos)>1 and any(near[u]==7 for u in pos):assert pos[1] in adj(G)[pos[0]]
            records.append({'name':item['name'],'marks':marks,'pos':pos,'near':near,'G':G})
            for u in pos:
                gamma=inverse_diag(G,u)
                rho=max(inverse_diag(deletion(G,{x}),u-(u>x)) for x in marks)
                assert gamma<=0
                if near[u]==7:assert gamma==0
                if gamma<0:assert rho>=F(16,5)
                tails.append({'core':item['name'],'marks':marks,'node':u,'distance':near[u],'gamma':str(gamma),'rho':str(rho)})
    assert len(tails)==21 and len(records)==16
    # Compare the current literal21-row data; same values as earlier162 table.
    old=json.loads((OUT/'expected_tail_positions.json').read_text())
    for new,expected in zip(tails,old):
        assert new['core']==expected['core'] and list(new['marks'])==expected['marks'] and new['node']==expected['node']
        assert new['distance']==expected['distance'] and new['gamma']==expected['zeta'] and new['rho']==expected['rho']
    candidates=[]
    for i,j in combinations_with_replacement(range(len(records)),2):
        L,M=records[i],records[j]
        for u in L['pos']:
            for v in M['pos']:
                if L['near'][u]+M['near'][v]>=13:candidates.append((i,j,((u,v),)))
        for u1,u2 in combinations(L['pos'],2):
            for v1,v2 in permutations(M['pos'],2):
                if u2 not in adj(L['G'])[u1] or v2 not in adj(M['G'])[v1]:continue
                if L['near'][u1]+M['near'][v1]<13 or L['near'][u2]+M['near'][v2]<13:continue
                ids=tuple(sorted(((u1,v1),(u2,v2))));inverse=tuple(sorted(((v1,u1),(v2,u2))))
                if i==j and ids>inverse:continue
                candidates.append((i,j,ids))
    assert len(candidates)==92
    rows=[];counts=Counter();scores=[]
    for i,j,ids in candidates:
        L,M=records[i],records[j];G,ix,okay=glue(L['G'],M['G'],ids)
        X=L['marks'];Y=tuple(ix[u] for u in M['marks']);allmarks=set(X+Y)
        row={'left':[L['name'],X],'right':[M['name'],M['marks']],'identification':ids}
        if not okay:status='conflicting_labels'
        elif positivity(deletion(G,allmarks))!=2:status='S_not_PD'
        else:
            terms=[];fail=False
            for x,y in product2(X,Y):
                removed=allmarks-{x,y};T=deletion(G,removed)
                if positivity(T)==0:fail=True;break
                kept=[u for u in range(G[0]) if u not in removed]
                terms.append(preliminary(T,kept.index(x),kept.index(y)))
            if fail:status='angle_not_PSD'
            else:
                status='survivor';score=sum(terms,F(0));assert score>=1;scores.append(score)
                row.update(score=str(score),pair_weights=[str(t) for t in terms])
        row['status']=status;counts[status]+=1;rows.append(row)
    assert counts=={'conflicting_labels':3,'S_not_PD':65,'angle_not_PSD':3,'survivor':21},counts
    assert min(scores)==1
    result={'tails':tails,'marked_records':len(records),'candidates':len(candidates),'counts':dict(counts),'minimum_score':str(min(scores)),'records':rows}
    (OUT/'tail_gluing_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','tails')},indent=2),flush=True)

def product2(X,Y):
    for x in X:
        for y in Y:yield x,y

if __name__=='__main__':run()
