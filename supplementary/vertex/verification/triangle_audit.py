"""Exact1265triangle marked subsets from the standard high-order catalogue."""
from fractions import Fraction as F
from itertools import combinations
from collections import deque,Counter
from pathlib import Path
import json
from tail_gluing_audit import graph,deletion,adj,ROOT,OUT

def pair(T,a,b):
    A=adj(T);parents={a:None};q=deque([a])
    while q:
        u=q.popleft()
        for v in A[u]:
            if v not in parents:parents[v]=u;q.append(v)
    if b not in parents:return None,F(0)
    p=[b]
    while p[-1]!=a:p.append(parents[p[-1]])
    p.reverse();d=len(p)-1
    special=any(len(A[u])>=3 for u in p) or any(A[u][v]>=4 for u,v in zip(p,p[1:]))
    q=F(1) if d<=3 and special else F(1,2) if d<=4 or special and d<=7 else F(0)
    return d,q

def run():
    source=json.loads((OUT/'high_cores.json').read_text())
    counts=Counter();records=[];minimum=F(99);minshort=99
    for item in source:
        G=graph(item['n'],item['edges'])
        for marks in combinations(range(G[0]),3):
            terms=[]
            for removed in marks:
                others=[u for u in marks if u!=removed];T=deletion(G,{removed});kept=[u for u in range(G[0]) if u!=removed]
                d,q=pair(T,kept.index(others[0]),kept.index(others[1]));terms.append((d,q))
            short=sum(d is not None and d<=7 for d,q in terms);score=sum((q for d,q in terms),F(0))
            assert short>=2 and score>=1,(item['name'],marks,terms)
            counts[G[0]]+=1;minimum=min(minimum,score);minshort=min(minshort,short)
            records.append({'diagram':item['name'],'marks':marks,'distances':[d for d,q in terms],'q':[str(q) for d,q in terms]})
    assert len(records)==1265
    result={'count':len(records),'counts_by_order':dict(counts),'minimum_auxiliary_sum':str(minimum),'minimum_short_pairs':minshort,'records':records}
    (OUT/'triangle_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2),flush=True)

if __name__=='__main__':run()
