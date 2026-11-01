"""Exact metric and weight data for every affine Coxeter family."""
from collections import Counter
from itertools import combinations
import json
import sympy as sp
from diagram_checks import ROOT,F,ALPHA,DELTA,graph,subgraph,total_weight

def finite(kind,r):
    if r==0: return graph(0,[])
    if kind in ('A','B'):
        return graph(r,[(i,i+1,4 if kind=='B' and i==0 else 3) for i in range(r-1)])
    if r==2: return graph(2,[])
    if r==3: return finite('A',3)
    return graph(r,[(i,i+1,3) for i in range(r-2)]+[(r-3,r-1,3)])

def loss(kind,r):
    if not r: return F(0)
    if r>=18: return {'A':F(41,3),'B':F(127,12),'D':F(23,4)}[kind]
    return ALPHA*r-total_weight(finite(kind,r))

def star(arms):
    edges=[]
    node=1
    for arm in arms:
        prev=0
        for _ in range(arm):
            edges.append((prev,node,3))
            prev=node
            node+=1
    return graph(node,edges)

def affine(kind,d):
    if kind=='A': return graph(d+1,[(i,i+1,3) for i in range(d)]+[(0,d,3)])
    if kind=='B': return graph(d+1,[(i,i+1,4 if i==d-2 else 3) for i in range(d-1)]+[(1,d,3)])
    if kind=='C': return graph(d+1,[(i,i+1,4 if i in (0,d-1) else 3) for i in range(d)])
    if kind=='D': return graph(d+1,[(i,i+1,3) for i in range(d-2)]+[(1,d-1,3),(d-3,d,3)])
    if kind=='E': return star({6:(2,2,2),7:(1,3,3),8:(1,2,5)}[d])
    if kind=='F': return graph(5,[(i,i+1,4 if i==1 else 3) for i in range(4)])
    if kind=='G': return graph(3,[(0,1,6),(1,2,3)])

def coordinates(kind,d):
    p=lambda i:[int(j<i) for j in range(d)]
    if kind=='C': return [p(i) for i in range(d+1)]
    if kind=='B': return [p(0),[2]+[0]*(d-1)]+[p(i) for i in range(2,d+1)]
    if kind=='D': return [p(0),[2]+[0]*(d-1)]+[p(i) for i in range(2,d-1)]+[p(d-1)[:-1]+[-1],p(d)]

def metric_data(kind,d):
    if kind=='A':
        lengths={(i,j):F((j-i)*(d+1-j+i),d+1) for i,j in combinations(range(d+1),2)}
    elif kind in ('B','C','D'):
        pts=coordinates(kind,d)
        lengths={(i,j):sum((x-y)**2 for x,y in zip(pts[i],pts[j])) for i,j in combinations(range(d+1),2)}
    else:
        G=affine(kind,d)
        M=sp.eye(d+1)*2
        for a,b,m in G[1]: M[a,b]=M[b,a]=-2*sp.cos(sp.pi/m)
        kernel=M.nullspace()
        assert len(kernel)==1
        a=kernel[0]
        assert all(x>0 for x in a) or all(x<0 for x in a)
        V=M[1:,1:].inv()
        lengths={}
        for i,j in combinations(range(d+1),2):
            val=V[j-1,j-1]/a[j]**2 if i==0 else V[i-1,i-1]/a[i]**2+V[j-1,j-1]/a[j]**2-2*V[i-1,j-1]/(a[i]*a[j])
            lengths[i,j]=sp.simplify(val)
    global_counts=Counter(lengths.values())
    S=[sum(v*v for v in Counter(t for (i,j),t in lengths.items() if z in (i,j)).values()) for z in range(d+1)]
    M=F(sum(v*v for v in global_counts.values()),(d+1)**2)
    return S,M,max(global_counts.values())

def deletion_node(kind,d,i):
    if kind in ('A','C'): return i
    if kind=='B': return 0 if i==0 else d if i==1 else i-1
    if kind=='D': return 0 if i==0 else d-1 if i==1 else d if i==d else i-1
    return i

def formula_loss(kind,d,i):
    if kind=='A': return loss('A',d)
    if kind=='C': return loss('B',i)+loss('B',d-i)
    if kind=='B': return loss('B',d) if i in (0,1) else loss('D',i)+loss('B',d-i)
    return loss('D',d) if i in (0,1,d-1,d) else loss('D',i)+loss('D',d-i)

def correction(kind,d):
    return 1 if kind=='A' and d==2 else 2 if kind=='B' and d>=9 else 3 if kind=='D' and d>=9 else 0

def run():
    finite_rows=[]
    for kind in ('A','B','D'):
        for r in range(2 if kind=='D' else 1,41):
            exact=ALPHA*r-total_weight(finite(kind,r))
            assert exact==loss(kind,r)
            assert exact>=DELTA
            finite_rows.append(dict(type=kind,order=r,loss=str(exact)))
    rows=[]
    for kind,minimum in (('A',2),('B',3),('C',2),('D',4)):
        for d in range(minimum,41):
            G=affine(kind,d)
            S,M,mult=metric_data(kind,d)
            losses=[]
            c=correction(kind,d)
            for i in range(d+1):
                node=deletion_node(kind,d,i)
                L=ALPHA*d-total_weight(subgraph(G,(node,)))
                assert L==formula_loss(kind,d,i),(kind,d,i,L,formula_loss(kind,d,i))
                residual=L+c-ALPHA*F(S[i],d)
                assert residual>=0,(kind,d,i,residual)
                losses.append(L)
                if d<=20:
                    rows.append(dict(type=kind,dimension=d,coordinate_vertex=i,diagram_node=node,
                                     loss=str(L),S=S[i],correction=c,residual=str(residual)))
            beta=sum(losses,F(0))/(d+1)-c
            assert beta>=ALPHA
            assert mult<=F(6*(d+1),5)
            if d>=9:
                assert beta>=4*ALPHA*M/d
                assert mult<=d+1
                T=F(d*(d+1)*(2*d+1),6)
                expected= F(d,2)-F(d%2,4) if kind=='A' else (T if kind=='C' else T-2*d+2 if kind=='B' else T-4*d+26)/(d+1)**2
                assert M==expected
            if d>=20:
                expected={'A':F(41,3),'B':F(49*d+673,3*(d+1)),
                          'C':F(127*d+1021,6*(d+1)),'D':F(23*d+557,2*(d+1))}[kind]
                assert beta+c==expected
    assert len(rows)==902
    exceptional=[]
    for kind,d in [('E',6),('E',7),('E',8),('F',4),('G',2)]:
        G=affine(kind,d)
        S,M,mult=metric_data(kind,d)
        losses=[ALPHA*d-total_weight(subgraph(G,(z,))) for z in range(d+1)]
        assert all(L>=ALPHA*F(s,d) for L,s in zip(losses,S))
        assert sum(losses,F(0))/(d+1)>=ALPHA
        assert mult<=F(6*(d+1),5)
        exceptional.append(dict(type=kind,dimension=d,S=S,loss_times_12=[int(L*12) for L in losses],maximum_multiplicity=mult))
    (ROOT/'tables/finite_losses.json').write_text(json.dumps(finite_rows,indent=2)+'\n')
    (ROOT/'tables/small_vertex_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    (ROOT/'tables/exceptional_metrics.json').write_text(json.dumps(exceptional,indent=2)+'\n')
    return dict(classical_vertex_checks=len(rows),exceptional_types=len(exceptional),formula_spot_checks_through_dimension=40)

if __name__=='__main__': print(run())
