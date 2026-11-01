"""Exact consistency checks for the normalized upper-half-space drawing of P^3.

The geometric argument proving equality of the two cusp lengths is in
Lemma 2.5. These checks verify the coordinate model used in Figure 2.
"""
from pathlib import Path
from datetime import datetime,timezone
import json
import sympy as s

ROOT=Path(__file__).absolute().parents[1]
x,y,z=s.symbols('x y z',real=True)
facets={'x0':x,'x1':1-x,'y0':y,'y1':1-y,
        'Sp':x*x+y*y+z*z-1,'Sq':(x-1)**2+(y-1)**2+z*z-1}
vertices={'p':(0,0,1),'q':(1,1,1),'u2':(1,0,0),'u3':(0,1,0)}
incidences={'p':{'x0','y0','Sp'},'q':{'x1','y1','Sq'},
            'u2':{'x1','y0','Sp','Sq'},'u3':{'x0','y1','Sp','Sq'}}
for name,point in vertices.items():
    values={k:s.simplify(f.subs(dict(zip((x,y,z),point)))) for k,f in facets.items()}
    assert all(v>=0 for v in values.values())
    assert {k for k,v in values.items() if v==0}==incidences[name]

edges=[('x0','y0',(0,0,2)),('x1','y1',(1,1,2)),
       ('x1','y0',(1,0,1)),('x0','y1',(0,1,1)),
       ('Sp','y0',(s.Rational(3,5),0,s.Rational(4,5))),
       ('Sp','x0',(0,s.Rational(3,5),s.Rational(4,5))),
       ('Sq','x1',(1,s.Rational(2,5),s.Rational(4,5))),
       ('Sq','y1',(s.Rational(2,5),1,s.Rational(4,5))),
       ('Sp','Sq',(s.Rational(1,2),s.Rational(1,2),s.sqrt(2)/2))]
for f,g,point in edges:
    subst=dict(zip((x,y,z),point))
    assert s.simplify(facets[f].subs(subst))==s.simplify(facets[g].subs(subst))==0
    assert all(s.simplify(k.subs(subst))>=0 for k in facets.values())
    nf=s.Matrix([s.diff(facets[f],v) for v in (x,y,z)])
    ng=s.Matrix([s.diff(facets[g],v) for v in (x,y,z)])
    assert s.simplify(nf.dot(ng).subs(subst))==0

# The two lower spherical triangles meet on x+y=1 and are orthogonal there.
assert s.expand(facets['Sq']-facets['Sp'])==2-2*x-2*y
grad_sp=s.Matrix([s.diff(facets['Sp'],v) for v in (x,y,z)])
grad_sq=s.Matrix([s.diff(facets['Sq'],v) for v in (x,y,z)])
assert s.expand(grad_sp.dot(grad_sq)-4*(facets['Sp']+1-x-y))==0

report={'status':'passed','utc':datetime.now(timezone.utc).isoformat(),
        'normalized_radius':1,'finite_vertices':['p','q'],'finite_boundary_ideal_vertices':['u2','u3'],
        'remaining_ideal_vertex':'infinity','vertex_facet_incidence_checks':4,
        'right_angle_edge_checks':9,'hemisphere_intersection_plane':'x+y=1',
        'horosphere_metric':'h^(-2)(dx^2+dy^2)','square_side_length':'1/h'}
(ROOT/'logs/square_link_geometry_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
