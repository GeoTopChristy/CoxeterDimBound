"""Independent local reconstruction of Bound82 point-pair budgets (2026-10-07).

No original Bound82 verification code is imported. Integer units are used for
the pair weights: w4=4w, q4=4q, j16=16J=4w4+q4.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json

OUT=Path(__file__).absolute().parent

def diagram(kind,r):
    if kind in ('A','B'):
        return r,[(i,i+1,4 if kind=='B' and i==0 else 3) for i in range(r-1)]
    if kind=='D':
        if r<=3: return diagram('A',r)
        return r,[(i,i+1,3) for i in range(r-2)]+[(r-3,r-1,3)]
    if kind=='a':
        if r==2: return r,[(0,1,999)]
        return r,[(i,i+1,3) for i in range(r-1)]+[(0,r-1,3)]
    if kind=='b':
        return r,[(i,i+1,4 if i==r-3 else 3) for i in range(r-2)]+[(1,r-1,3)]
    if kind=='c':
        return r,[(i,i+1,4 if i in (0,r-2) else 3) for i in range(r-1)]
    if kind=='d':
        return r,[(i,i+1,3) for i in range(r-3)]+[(1,r-2,3),(r-4,r-1,3)]
    if kind.startswith('E'):
        arms={'E6':(2,2,2),'E7':(1,3,3),'E8':(1,2,5)}[kind]
        edges=[]; node=1
        for arm in arms:
            p=0
            for _ in range(arm):
                edges.append((p,node,3));p=node;node+=1
        return node,edges
    if kind=='F4': return 5,[(i,i+1,4 if i==1 else 3) for i in range(4)]
    if kind=='G2': return 3,[(0,1,6),(1,2,3)]
    raise ValueError(kind)

def budget(kind,r):
    n,edges=diagram(kind,r)
    adj=[{} for _ in range(n)]
    for a,b,m in edges: adj[a][b]=adj[b][a]=m
    cyclic=len(edges)==n and all(len(z)==2 for z in adj)
    W4=Q4=0
    def arm(u,prev):
        for _ in range(5):
            onward=[v for v in adj[u] if v!=prev]
            if len(onward)!=1 or adj[u][onward[0]]!=3: return False
            prev,u=u,onward[0]
        return True
    for a in range(n):
        parents={a:None}; queue=deque([a])
        while queue:
            u=queue.popleft()
            for v in adj[u]:
                if v not in parents: parents[v]=u;queue.append(v)
        for b in range(a+1,n):
            p=[b]
            while p[-1]!=a: p.append(parents[p[-1]])
            p.reverse();dist=len(p)-1
            special=any(len(adj[u])>=3 for u in p) or any(adj[u][v]>=4 for u,v in zip(p,p[1:]))
            w4=4 if dist<=12 else (4 if special else 2) if dist==13 else (4 if special else 1) if dist==14 else 0
            q4=(4 if special else 2) if dist<=3 else 2 if dist==4 or special and dist<=7 else 0
            if cyclic and dist>=9: w4=0
            elif not cyclic and dist>=8 and arm(p[0],p[1]) and arm(p[-1],p[-2]): w4=0
            W4+=w4;Q4+=q4
    return F(W4,4),F(Q4,4),F(4*W4+Q4,16)

def run():
    formula={'A':(F(59,2),-5,F(113,4)), 'B':(F(73,2),-2,F(36)), 'D':(F(183,4),5,F(47)),
             'c':(F(87,2),1,F(175,4)), 'b':(F(211,4),8,F(219,4)), 'd':(F(62),15,F(263,4))}
    checks=0
    for r in list(range(28,50))+[64,80,100,150]:
        for kind,(wc,qc,jc) in formula.items():
            actual=budget(kind,r)
            expected=7*r+wc,2*r+qc,F(15*r,2)+jc
            assert actual==expected,(kind,r,actual,expected)
            checks+=1
    # Small budgets in Appendix D, reconstructing all three finite families.
    v={r:F(5,4)*F(r*(r-1),2) for r in range(1,9)}
    expected_v=['333/8','205/4','247/4','293/4','343/4','397/4','455/4','517/4','573/4','313/2','679/4','183','781/4','413/2','867/4','226','937/4','242','499/2']
    finite=[]
    for r in range(9,90):
        scores={kind:budget(kind,r)[2] for kind in ('A','B','D')}
        v[r]=max(scores.values())
        if r<=27: assert v[r]==F(expected_v[r-9]),(r,v[r])
        else: assert v[r]==F(15*r,2)+47
        finite.append({'r':r,'scores':{k:str(a) for k,a in scores.items()},'bound':str(v[r])})
    psi={d:F(5*d,8) for d in range(6)}
    expected_psi=['103/28','271/64','43/9','85/16','511/88','605/96','353/52','407/56','31/4','527/64','593/68','647/72','705/76','761/80','407/42','867/88','229/23','961/96','501/50','1039/104','268/27']
    affine=[]
    for d in range(2,100):
        r=d+1
        scores={'a':budget('a',r)[2]/r,'c':budget('c',r)[2]/r}
        if d>=3: scores['b']=budget('b',r)[2]/r
        if d>=4: scores['d']=budget('d',r)[2]/r
        if d>=6:
            psi[d]=max(scores.values())
            if d<=26: assert psi[d]==F(expected_psi[d-6]),(d,psi[d])
            else: assert psi[d]==max(F(17,2),F(15,2)+F(263,4*(d+1)))
        assert scores['a']==min(F(d,2),8)+F(1,4)*min(F(d,4),2)
        affine.append({'d':d,'scores':{k:str(a) for k,a in scores.items()}})
    expected_exc={'E6':F(51,14),'E7':F(133,32),'E8':F(337,72),'F4':F(19,8),'G2':F(29,24)}
    exceptional={}
    for kind,val in expected_exc.items():
        r=diagram(kind,0)[0]
        actual=budget(kind,0)[2]/r
        assert actual==val,(kind,actual,val)
        assert actual<=psi[r-1]
        exceptional[kind]=str(actual)
    # Exact finite-vertex dynamic program.
    D=[F(0)]; winners=[[]]
    for n in range(1,90):
        pairs=[(D[n-r]+v[r],r) for r in range(1,n+1)]
        val,r=max(pairs)
        D.append(val);winners.append(winners[n-r]+[r])
    expected_fin=['1513/164','3075/332','781/84','3169/340','1607/172','3259/348','413/44','3345/356']
    dp=[]
    for n in range(82,90):
        assert D[n]/n==F(expected_fin[n-82]),(n,D[n]/n)
        dp.append({'n':n,'bound':str(D[n]/n),'witness':winners[n]})
    report={'stable_formula_checks':checks,'classical_affine_checks':len(affine),'exceptional':exceptional,'finite':finite,'affine':affine,'finite_DP':dp}
    (OUT/'budget_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'stable_formula_checks':checks,'finite_DP':dp,'exceptional':exceptional},indent=2),flush=True)
    return psi,v

if __name__=='__main__':run()
