"""Check the explicit coordinates and disjoint cases in Lemma 2.5.

This finite arithmetic check supplements the symbolic counting argument.
It uses literal standard-basis coordinates, including all five exceptions.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
import sympy as sp


def exceptional_vertices():
    def vec(n, *values):
        return sp.Matrix(list(map(sp.sympify, values)) + [0]*(n-len(values)))
    s=sp.sqrt
    return {
        'E6': [vec(6),vec(6,1),vec(6,1,s(3)),vec(6,0,0,1),
               vec(6,0,0,1,s(3)),vec(6,0,0,0,0,1),vec(6,0,0,0,0,1,s(3))],
        'E7': [vec(7),vec(7,s(3)),vec(7,0,s(2)),vec(7,0,s(2),2),
               vec(7,0,s(2),2,2*s(3)),vec(7,0,0,0,0,s(2)),
               vec(7,0,0,0,0,s(2),2),vec(7,0,0,0,0,s(2),2,2*s(3))],
        'E8': [vec(8),vec(8,2*s(5)),vec(8,0,s(15)),vec(8,0,s(15),3*s(5)),
               vec(8,0,0,0,2*s(3)),vec(8,0,0,0,2*s(3),3*s(2)),
               vec(8,0,0,0,2*s(3),3*s(2),s(30)),
               vec(8,0,0,0,2*s(3),3*s(2),s(30),2*s(15)),
               vec(8,0,0,0,2*s(3),3*s(2),s(30),2*s(15),6*s(5))],
        'F4': [vec(4),vec(4,s(3)),vec(4,s(3),1),vec(4,s(3),1,s(2)),
               vec(4,s(3),1,s(2),s(6))],
        'G2': [vec(2),vec(2,1),vec(2,1,s(3))],
    }


def expected_edges(name):
    if name.startswith('E'):
        arms={'E6':[2,2,2],'E7':[1,3,3],'E8':[1,2,5]}[name]
        edges={};i=1
        for arm in arms:
            parent=0
            for _ in range(arm):
                edges[(parent,i)]=3;parent=i;i+=1
        return edges
    if name=='F4':return {(0,1):3,(1,2):4,(2,3):3,(3,4):3}
    return {(0,1):6,(1,2):3}


def check_exceptional():
    from metric_checks import metric_data
    result={}
    for name,vs in exceptional_vertices().items():
        n=len(vs)-1
        basis=sp.Matrix.hstack(*vs[1:])
        inverse=basis.inv()
        gradients=sp.Matrix.vstack(-sp.ones(1,n)*inverse,inverse)
        gram=gradients*gradients.T
        edges=expected_edges(name)
        for i,j in combinations(range(n+1),2):
            m=edges.get((i,j),2)
            assert gram[i,j]<=0
            assert sp.simplify(gram[i,j]**2/(gram[i,i]*gram[j,j])-sp.cos(sp.pi/m)**2)==0,(name,i,j)
        counts=Counter(sp.simplify((v-u).dot(v-u)) for u,v in combinations(vs,2))
        assert max(counts.values())=={'E6':6,'E7':3,'E8':3,'F4':2,'G2':1}[name]
        _,moment,mx=metric_data(name[0],n)
        assert mx==max(counts.values())
        assert sp.Rational(sum(x*x for x in counts.values()),(n+1)**2)==sp.Rational(moment.numerator,moment.denominator)
        result[name]={'maximum':mx,'squared_length_multiplicities':{str(k):v for k,v in sorted(counts.items())},'facet_angles':'exactly verified'}
    return result


def classical_vertices(kind,n):
    p=lambda k:tuple(int(i<k) for i in range(n))
    if kind=='C':return [p(i) for i in range(n+1)]
    first=[p(0),tuple(2*x for x in p(1))]
    if kind=='B':return first+[p(i) for i in range(2,n+1)]
    return first+[p(i) for i in range(2,n-1)]+[tuple([1]*(n-1)+[-1]),p(n)]


def squared(u,v):return sum((x-y)**2 for x,y in zip(u,v))


def table_count(kind,n,k):
    if kind=='C':return n+1-k if 1<=k<=n else 0
    if kind=='B':
        if k in (1,4):return n-2
        if 2<=k<=n-1 and k!=4:return n+1-k
        return 2 if k==n else 0
    if k==1:return n-4
    if k==4:return n-1
    if 2<=k<=n-3 and k!=4:return n+1-k
    return 4 if k in (n-2,n) else 0


def check_classical():
    checked=0
    for kind,low in [('C',2),('B',3),('D',4)]:
        for n in range(low,41):
            vs=classical_vertices(kind,n)
            counts=Counter(squared(u,v) for u,v in combinations(vs,2))
            if n>=7 or kind=='C':
                assert all(counts[k]==table_count(kind,n,k) for k in range(1,n+2)),(kind,n)
            maximum=n if kind=='C' else n-1 if kind=='B' else {4:6,5:4,6:6}.get(n,n-1)
            assert max(counts.values())==maximum,(kind,n)
            assert maximum<=sp.Rational(6,5)*(n+1)
            if n>=9:
                base=sp.Rational(n*(n+1)*(2*n+1),6)
                expected={'C':base,'B':base-2*n+2,'D':base-4*n+26}[kind]
                assert sum(c*c for c in counts.values())==expected,(kind,n,'moment')
            for i,v in enumerate(vs):
                inc=Counter(squared(u,v) for j,u in enumerate(vs) if j!=i)
                S=sum(c*c for c in inc.values())
                if kind=='C':expect=n+2*min(i,n-i)
                elif kind=='B':
                    if i<2:expect=3 if n==3 else n+2
                    elif 2*i<=n:expect=n+6+2*min(i-2,n-i)
                    else:expect=n+2+2*min(i-2,n-i)
                elif i in (0,1,n-1,n):expect=n+6 if n==4 else n+2 if n==5 else n+4
                else:expect=n+4+2*min(i-2,n-2-i)+(8 if 2*i==n else 4 if abs(2*i-n)>=2 else 0)
                assert S==expect,(kind,n,i,S,expect)
            checked+=1
    # Derive the three moment numerators symbolically from the disjoint cases.
    n,k=sp.symbols('n k',integer=True,positive=True)
    c=sp.summation((n+1-k)**2,(k,1,n))
    b=2*(n-2)**2+sp.summation((n+1-k)**2,(k,2,n-1))-(n-3)**2+4
    d=(n-4)**2+(n-1)**2+sp.summation((n+1-k)**2,(k,2,n-3))-(n-3)**2+32
    assert sp.simplify(b-c+2*n-2)==0
    assert sp.simplify(d-c+4*n-26)==0
    return {'dimensions_checked_through':40,'simplices_checked':checked,'global_and_incident_counts':'passed','symbolic_moment_identities':'passed'}


def check_a_and_correction():
    for n in range(2,41):
        N=n+1
        counts=Counter(sp.Rational((j-i)*(N-j+i),N) for i,j in combinations(range(N),2))
        assert max(counts.values())==N
        moment=sp.Rational(sum(c*c for c in counts.values()),N*N)
        assert moment==sp.Rational(n,2)-(sp.Rational(1,4) if n%2 else 0)
    C=sp.Rational(73,12)
    for n in range(1,9):
        for correction in ([0,1] if n==2 else [0]):
            assert n*(C-sp.Rational(n-1,2))-correction>=C
    return 'passed'


def check_classical_angles():
    from metric_checks import affine
    count=0
    for kind,low in [('A',2),('B',3),('C',2),('D',4)]:
        for n in range(low,13):
            if kind=='A':
                N=n+1
                vs=[sp.Matrix([int(j<i)-sp.Rational(i,N) for j in range(N)]) for i in range(N)]
            else:vs=list(map(sp.Matrix,classical_vertices(kind,n)))
            basis=sp.Matrix.hstack(*vs[1:])
            h=sp.Matrix.vstack(-sp.ones(1,n),sp.eye(n))
            gram=h*(basis.T*basis).inv()*h.T
            diagram={(a,b):m for a,b,m in affine(kind,n)[1]}
            if kind=='B':mapping=[0,n]+list(range(1,n))
            elif kind=='D':mapping=[0,n-1]+list(range(1,n-1))+[n]
            else:mapping=list(range(n+1))
            for i,j in combinations(range(n+1),2):
                m=diagram.get(tuple(sorted((mapping[i],mapping[j]))),2)
                assert gram[i,j]<=0
                assert sp.simplify(gram[i,j]**2/(gram[i,i]*gram[j,j])-sp.cos(sp.pi/m)**2)==0,(kind,n,i,j,m)
            count+=1
    return {'simplices_checked':count,'dimensions_through':12,'facet_angles':'exactly verified'}


def main():
    result={'arithmetic':'exact','exceptional_coordinates':check_exceptional(),'classical':check_classical(),'A_moment_and_small_correction':check_a_and_correction(),'classical_angles':check_classical_angles()}
    path=Path(__file__).absolute().parents[1]/'tables/coordinate_proof_checks.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
