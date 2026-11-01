"""Low-facet witnesses and exact symbolic elimination."""
from itertools import combinations,product
import json
import networkx as nx
import sympy as s
from diagram_checks import ROOT

def faces(embedding):
    seen=set()
    out=[]
    for a,b in embedding.edges():
        if (a,b) not in seen:
            out.append(tuple(embedding.traverse_face(a,b,seen)))
    return out

def check_faces():
    rows=[]
    for G in nx.graph_atlas_g():
        n=len(G)
        if not 4<=n<=6 or min(dict(G.degree()).values())<3 or nx.node_connectivity(G)<3:
            continue
        planar,embedding=nx.check_planarity(G)
        if not planar: continue
        ff=faces(embedding)
        if any(len(f)>4 for f in ff): continue
        cofacial=lambda marks:any(set(marks)<=set(f) for f in ff)
        II=[]
        III=[]
        IV=[]
        for triple in combinations(G.nodes(),3):
            if cofacial(triple): continue
            if all(G.has_edge(a,b) for a,b in combinations(triple,2)):
                II.append(triple)
            for a in triple:
                b,c=[v for v in triple if v!=a]
                if G.has_edge(a,b) and G.has_edge(a,c) and not G.has_edge(b,c) and cofacial((b,c)):
                    III.append((a,b,c))
        for marks in combinations(G.nodes(),4):
            for X in combinations(marks,2):
                Y=tuple(v for v in marks if v not in X)
                if X>Y: continue
                if not cofacial(X) and not cofacial(Y) and all(G.has_edge(a,b) for a in X for b in Y):
                    IV.append((X,Y))
        degrees=sorted(dict(G.degree()).values())
        witness=[]
        if n==4:
            kind='I'; bound=3; witness=[list(G.nodes())]
        elif II:
            kind='II'; bound=2; witness=[II[0]]
        elif n==5:
            kind='two III'; bound=2
            for a,b in combinations(III,2):
                ee=lambda t:{tuple(sorted((t[0],v))) for v in t[1:]}
                if not ee(a)&ee(b): witness=[a,b]; break
            assert witness
        elif degrees==[3]*6:
            kind='exceptional bipyramid'; bound=0
        elif III:
            kind='III'; bound=1; witness=[III[0]]
        else:
            kind='three IV'; bound=1
            for triple in combinations(IV,3):
                edge_sets=[{tuple(sorted((a,b))) for a in X for b in Y} for X,Y in triple]
                if len(set.union(*edge_sets))==12 and sum(map(len,edge_sets))==12:
                    witness=list(triple); break
            assert witness
        assert bound>=7-n or kind=='exceptional bipyramid'
        rows.append(dict(facets=n,dual_edges=sorted(G.edges()),dual_faces=ff,dual_degrees=degrees,
                         configuration=kind,witness=witness,lower_bound=bound))
    assert len(rows)==9
    for rows_n,expected in [(4,1),(5,2),(6,6)]:
        assert sum(r['facets']==rows_n for r in rows)==expected
    (ROOT/'tables/three_faces.json').write_text(json.dumps(rows,indent=2)+'\n')
    # Exhaustive colourings for the two small rectangle bounds.
    for a,b,bound in [(2,3,2),(3,3,4)]:
        max_rect=0
        for data in product((0,1),repeat=a*b):
            ct=sum(data[i*b+k]==data[j*b+l] and data[i*b+l]==data[j*b+k]
                   and data[i*b+k]!=data[i*b+l]
                   for i,j in combinations(range(a),2) for k,l in combinations(range(b),2))
            max_rect=max(max_rect,ct)
        assert max_rect==bound
    return len(rows)

def check_algebra():
    n,t,j=s.symbols('n t j')
    eps=2*(5*n-1)/((n-1)*(n-4)*(9*n-1))
    kap=5*(n-1)**2*(n-2)/(4*n*(n-4)*(9*n-1))
    lam=2/(n-2)*(s.Rational(73,12)-s.Rational(23,4)/(n-1))
    H=9*n**3-46*n**2+11*n+2
    assert s.cancel(1-3*eps-H/((n-1)*(n-4)*(9*n-1)))==0
    polynomials={
        'H':H,
        'P_e':9*n**4-1486*n**3+4335*n**2-2994*n+256,
        'Q_e':5*n**4-457*n**3+165*n**2-11*n+10,
        'P_o':9*n**5-1513*n**4+8901*n**3-12393*n**2+1384*n+60,
        'Q_o':5*n**5-472*n**4+1968*n**3-986*n**2+91*n-30}
    for parity,u in [('e',6*(n-2)/(n-4)),('o',6*(n-1)/(n-3))]:
        U=(u-6*eps)/(1-3*eps)
        r=(kap+3*eps)/(1-3*eps)
        F0=7-U-lam*(3*U-6)
        F1=(1+3*lam)*r+12*lam-1
        denom=(n-2)*H*(1 if parity=='e' else n-3)
        assert s.cancel(F0-polynomials['P_'+parity]/denom)==0
        assert s.cancel(F0+F1-(2*n**2+67*n-138)*polynomials['Q_'+parity]/(8*n*(n-1)*denom))==0
        c=u/(n-2)
        choose=lambda z,k:s.prod(z-i for i in range(k))/s.factorial(k)
        coeff=c*(choose(j,3)+choose(n-j,3))-choose(j,2)-choose(n-j,2)
        expected=(2*j-n)**2*(n-1)/(2*(n-4)) if parity=='e' else n*(2*j-n-1)*(2*j-n+1)/(2*(n-3))
        assert s.factor(coeff-expected)==0,(parity,s.factor(coeff))
        local=choose(j,2)-c*choose(j,3)
        expected=(n-2*j)*j*(j-1)/(2*(n-4))
        if parity=='o': expected+=2*j*(j-1)*(j-2)/((n-4)*(n-3)*(n-2))
        assert s.factor(local-expected)==0
    shifts={key:s.Poly(s.expand(value.subs(n,163+t)),t).all_coeffs() for key,value in polynomials.items()}
    assert all(coefficient>0 for coefficients in shifts.values() for coefficient in coefficients)
    (ROOT/'tables/positive_shifts.json').write_text(json.dumps({k:[int(x) for x in v] for k,v in shifts.items()},indent=2)+'\n')
    return 'all polynomial identities and positive shifts passed'

if __name__=='__main__':
    print('three-faces',check_faces())
    print(check_algebra())
