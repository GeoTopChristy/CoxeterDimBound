"""Reconstruct the 192 cyclic Schur witnesses from Appendix C."""
from fractions import Fraction as F
from pathlib import Path
from collections import deque
from itertools import product
import json

OUT=Path(__file__).absolute().parent
MODES=('O','E+','E-','H+','H-','N+','N-','B')
SHORT={(9,9):{('E+','N-'),('E-','N+'),('N+','E-'),('N-','E+'),('N+','N-'),('N-','N+')},(9,10):{('N+','N-')},(10,10):set()}

def support(mode,z,N):
    if mode=='O':return {},2,1
    if mode=='B':return {(z-1)%N:-1,(z+1)%N:-1},2,0
    step=(1 if mode[0] in ('E','H') else 2)*(1 if mode[1]=='+' else -1)
    mag=2 if mode[0]=='H' else 1
    return {(z+step)%N:-mag},4 if mag==2 else 2,0

def schur(L,M,mc,md):
    N=L+M;a=0;b=L;c=N;d=N+1
    ret=(a,b,c,d)
    A=[[F(0) for _ in range(N+2)] for _ in range(N+2)]
    for i in range(N):
        A[i][i]=F(2);j=(i+1)%N;A[i][j]=A[j][i]=F(-1)
    for new,z,mode in ((c,a,mc),(d,b,md)):
        s,diag,opp=support(mode,z,N);A[new][new]=F(diag)
        for old,val in s.items():A[new][old]=A[old][new]=F(val)
        A[new][z]=A[z][new]=F(-opp)
    un=[i for i in range(N) if i not in (a,b)]
    conn=[{j:A[i][j] for j in un if A[i][j]} for i in ret]
    def inverse(i,j):
        if 0<i<L and 0<j<L:
            return F(min(i,j)*(L-max(i,j)),L)
        if L<i<N and L<j<N:
            ii=i-L;jj=j-L
            return F(min(ii,jj)*(M-max(ii,jj)),M)
        return F(0)
    R=[[A[ret[i]][ret[j]]-sum(x*inverse(u,v)*y for u,x in conn[i].items() for v,y in conn[j].items()) for j in range(4)] for i in range(4)]
    return A,R,ret

def restricted(R,x,y):
    p=R[0][0]+2*x*R[0][2]+x*x*R[2][2]
    q=R[1][1]+2*y*R[1][3]+y*y*R[3][3]
    r=R[0][1]+y*R[0][3]+x*R[2][1]+x*y*R[2][3]
    return p,q,r,p*q-r*r

def paths(A,ret):
    pairs=((0,1),(1,2),(2,3),(3,0))
    out=[]
    for i,j in pairs:
        keep=set(range(len(A)))-(set(ret)-{ret[i],ret[j]})
        adj={u:[v for v in keep if v!=u and A[u][v]] for u in keep}
        src=ret[i];dst=ret[j]
        queue=deque([src]);parents={src:None}
        while queue:
            u=queue.popleft()
            for v in adj[u]:
                if v not in parents:parents[v]=u;queue.append(v)
        if dst not in parents:out.append(None);continue
        p=[dst]
        while p[-1]!=src:p.append(parents[p[-1]])
        p.reverse()
        comp=set(parents)
        cycle=all(len(adj[u])==2 for u in comp)
        out.append({'marks':[i,j],'distance':len(p)-1,'noncyclic':not cycle,'path':p,'degrees':{u:len(adj[u]) for u in p}})
    return out

def run():
    rows=[];summaries=[]
    for L,M in ((10,10),(9,9),(9,10)):
        initials=[];diags=[];dets=[]
        for mc,md in product(MODES,repeat=2):
            A,R,ret=schur(L,M,mc,md)
            row={'L':L,'M':M,'modes':[mc,md]}
            if (mc,md) in SHORT[L,M]:
                eligible=[p for p in paths(A,ret) if p and p['noncyclic'] and p['distance']<=9]
                assert eligible,(L,M,mc,md)
                row.update(outcome='short',paths=eligible)
            else:
                x=y=F(3,5) if (L,M)==(10,10) else F(1,2)
                if (L,M)==(9,9) and (mc,md) in {('N+','N+'),('N-','N-')}:x,y=F(1,2),F(3,5)
                if (L,M)==(9,10):
                    if (mc,md)==('E+','N-'):x,y=F(1,2),F(3,5)
                    elif (mc,md)==('N+','E-'):x,y=F(3,5),F(1,2)
                    elif (mc,md)==('N-','N+'):x=y=F(3,5)
                p,q,r,det=restricted(R,x,y)
                assert p<0 and det>0,(L,M,mc,md,x,y,p,q,r,det)
                row.update(outcome='negative',x=str(x),y=str(y),diagonal=str(p),determinant=str(det),matrix=[[str(p),str(r)],[str(r),str(q)]])
                diags.append(p);dets.append(det)
            rows.append(row)
        expected={(9,9):(F(-5,18),F(1,81),58),(9,10):(F(-13,45),F(8,1125),63),(10,10):(F(-7,25),F(24,625),64)}[L,M]
        assert (max(diags),min(dets),len(dets))==expected,(L,M,max(diags),min(dets),len(dets),expected)
        summaries.append({'L':L,'M':M,'negative_witnesses':len(dets),'short':64-len(dets),'max_first_diagonal':str(max(diags)),'minimum_determinant':str(min(dets))})
    assert len(rows)==192
    result={'summaries':summaries,'records':rows}
    (OUT/'cycle_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'summaries':summaries,'short_records':[row for row in rows if row['outcome']=='short']},indent=2),flush=True)

if __name__=='__main__':run()
