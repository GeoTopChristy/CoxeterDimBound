"""Replay tree certificates with independent closed-form block inverses.

Does not import tree_audit.py. Uses no network or third-party library.
"""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter
import hashlib,json,time,platform

OUT=Path(__file__).absolute().parent

def rebuild(row):
    P,k,Q=row['P'],row['ell'],row['Q'];m=P+k+Q;a=P;b=P+k
    kinds=row['endtype'];N=m+1;leaves=[];edges={};diagonal={i:F(2) for i in range(N)}
    for i in range(m):edges[i,i+1]=F(-1)
    for endpoint,kind in zip((0,m),kinds):
        current=[]
        for _ in range({'A':0,'B':1,'D':2}[kind]):
            current.append(N);diagonal[N]=F(4 if kind=='B' else 2)
            edges[min(endpoint,N),max(endpoint,N)]=F(-2 if kind=='B' else -1);N+=1
        leaves.append(current)
    c,d=N,N+1
    for z,support,variant,opposite in ((c,row['support_c'],row['variant_c'],a),(d,row['support_d'],row['variant_d'],b)):
        heavy=any(label==4 for _,label in support)
        diagonal[z]=F(4 if heavy else 2)-F(1,2)*variant[1]
        for u,label in support:edges[min(u,z),max(u,z)]=F(-2 if heavy else -1)
        if variant[0]:edges[min(z,opposite),max(z,opposite)]=-F(variant[0])
    def entry(i,j):return diagonal[i] if i==j else edges.get((min(i,j),max(i,j)),F(0))
    # Three disconnected unmarked blocks. Trunks at B/D ends are ordered
    # from the decorated endpoint toward the deleted marked node.
    blocks=[(list(range(a)),kinds[0],leaves[0]),(list(range(a+1,b)),'A',[]),(list(range(m,b,-1)),kinds[1],leaves[1])]
    inverse={}
    for trunk,kind,extra in blocks:
        t=len(trunk)
        for i,u in enumerate(trunk):
            for j,v in enumerate(trunk):
                inverse[u,v]=(F((min(i,j)+1)*(t-max(i,j)),t+1) if kind=='A' else F(t-max(i,j)))
        if kind!='A':
            for u in extra:
                for j,v in enumerate(trunk):inverse[u,v]=inverse[v,u]=F(t-j,2)
                for v in extra:
                    inverse[u,v]=(F(t+1,4) if kind=='B' else F(t,4)+(F(1,2) if u==v else 0))
    unmarked=sorted({u for block,_,extra in blocks for u in block+extra})
    retained=[a,b,c,d]
    conn=[{u:entry(i,u) for u in unmarked if entry(i,u)} for i in retained]
    R=[[entry(i,j)-sum(x*inverse.get((u,v),F(0))*y for u,x in conn[ii].items() for v,y in conn[jj].items()) for jj,j in enumerate(retained)] for ii,i in enumerate(retained)]
    return R,edges,retained

def run():
    started=time.time();counts=Counter();bytype={};seen=set();chains=set();minimum=None
    for line in (OUT/'records.jsonl').open(encoding='utf-8'):
        row=json.loads(line);ident=tuple(json.dumps(row[key],sort_keys=True) for key in ('endtype','P','ell','Q','support_c','support_d','variant_c','variant_d'))
        assert ident not in seen;seen.add(ident)
        chains.add((row['endtype'],row['P'],row['ell'],row['Q']))
        R,edges,ret=rebuild(row);cert=row['certificate'];outcome=row['outcome']
        if outcome=='failed':
            if cert['condition']=='S+mark PD':
                z=cert['mark'];assert R[z][z]==F(cert['pivot']) and R[z][z]<=0
            else:
                i,j=cert['pair'];A,B,C=R[i][i],R[i][j],R[j][j]
                assert [[A,B],[B,C]]==[[F(v) for v in rr] for rr in cert['matrix']]
                assert A<0 or C<0 or A*C-B*B<0
        elif outcome=='negative':
            x,y=F(cert['x']),F(cert['y']);assert x>0 and y>0
            U=[[F(1),F(0),x,F(0)],[F(0),F(1),F(0),y]]
            M=[[sum(U[i][a]*R[a][b]*U[j][b] for a in range(4) for b in range(4)) for j in range(2)] for i in range(2)]
            det=M[0][0]*M[1][1]-M[0][1]*M[1][0]
            assert M==[[F(v) for v in rr] for rr in cert['matrix']]
            assert M[0][0]<0 and det>0 and det==F(cert['determinant'])
            minimum=det if minimum is None else min(minimum,det)
        elif outcome=='short':
            convert=lambda u:ret[2] if u=='c' else ret[3] if u=='d' else u
            path=[convert(u) for u in cert['path']];pair=[convert(u) for u in cert['pair']]
            assert path[0]==pair[0] and path[-1]==pair[1] and len(path)-1==cert['distance']<=7
            assert len(path)==len(set(path)) and not ((set(ret)-set(pair))&set(path))
            P,k,Q=row['P'],row['ell'],row['Q'];a,b=ret[:2];m=P+k+Q
            grow=set()
            for length,l,r in ((P,2,a-2),(k,a+2,b-2),(Q,b+2,m-2)):
                boundary=10 if (l,r)==(a+2,b-2) else 9
                if length==boundary:grow.update((i,i+1) for i in range(l,r))
            for u,v in zip(path,path[1:]):
                edge=(min(u,v),max(u,v));assert edges.get(edge) and edge not in grow
        else:raise AssertionError(outcome)
        counts[outcome]+=1;bytype.setdefault(row['endtype'],Counter())[outcome]+=1
    assert len(chains)==450 and len(seen)==79113
    assert counts=={'negative':39389,'short':11167,'failed':28557}
    assert minimum==F(1,32400)
    result={'independent_inverse_formulas':True,'all_certificates_replayed':True,'distinct_initial_chains':len(chains),'distinct_records':len(seen),'counts':dict(counts),'bytype':{k:dict(v) for k,v in bytype.items()},'minimum_determinant':str(minimum),'seconds':round(time.time()-started,2),'record_sha256':hashlib.sha256((OUT/'records.jsonl').read_bytes()).hexdigest(),'runtime':{'python':platform.python_version(),'platform':platform.platform(),'dependencies':'Python standard library only'}}
    (OUT/'replay_summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':run()
