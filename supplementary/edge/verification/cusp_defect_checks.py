from pathlib import Path
from datetime import datetime,timezone
import sys,json
import sympy as s

OUT=Path(__file__).absolute().parents[1]
sys.path.insert(0,str(OUT/'verification'))
from global_checks import check_algebra

n,j=s.symbols('n j')
choose=lambda x,k:s.prod(x-i for i in range(k))/s.factorial(k)
for parity,u,c in [('even',6+12/(n-4),6/(n-4)),
                   ('odd',6+12/(n-3),6*(n-1)/((n-2)*(n-3)))]:
    assert s.cancel(u/(n-2)-c)==0
    pair=u/(n-2)*(choose(j,3)+choose(n-j,3))-choose(j,2)-choose(n-j,2)
    expected=(2*j-n)**2*(n-1)/(2*(n-4)) if parity=='even' else n*(2*j-n-1)*(2*j-n+1)/(2*(n-3))
    assert s.factor(pair-expected)==0
    kernel=choose(j,2)-u/(n-2)*choose(j,3)
    expected=(n-2*j)*j*(j-1)/(2*(n-4))
    if parity=='odd':expected+=2*j*(j-1)*(j-2)/((n-4)*(n-3)*(n-2))
    assert s.factor(kernel-expected)==0
t=s.symbols('t')
for dim in range(6,25):
    poly=s.Poly(s.expand((1+t+t*t)*sum(t**i for i in range(dim-2))),t)
    h=[poly.nth(i) for i in range(dim)]
    g=[h[0]]+[h[i]-h[i-1] for i in range(1,(dim-1)//2+1)]
    assert g[:3]==[1,1,1] and all(x==0 for x in g[3:]),dim
    u=6+s.Rational(12,dim-4 if dim%2==0 else dim-3)
    theta=sum(g[k]*(s.binomial(k,2)-u/(dim-2)*s.binomial(k,3)) for k in range(2,len(g)))
    assert theta==1,dim
result=check_algebra()
report={'status':'passed','utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Section 5 notation equivalence, coefficient identities, link example, and existing endpoint algebra; no full-suite rerun.',
        'checks':{'u_n_normalization_both_parities':True,'paired_coefficients_both_parities':True,
                  'cusp_kernel_both_parities':True,'product_link_example_dimensions':'6 through 24',
                  'existing_endpoint_algebra':result}}
(OUT/'logs/cusp44_derivation_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
