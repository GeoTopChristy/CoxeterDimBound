from fractions import Fraction as F
from pathlib import Path
import json, math

def parts(n,lo=1):
    if n==0:
        yield ()
    for x in range(lo,n+1):
        for tail in parts(n-x,x):
            yield (x,)+tail

def calc(ds):
    h=[1]
    for d in ds:
        out=[0]*(len(h)+d)
        for i,x in enumerate(h):
            for j in range(d+1): out[i+j]+=x
        h=out
    R=sum(ds);H=sum(h)
    g=[h[j]-(h[j-1] if j else 0) for j in range(R//2+1)]
    assert all(x>=0 for x in g)
    pi=[F((R+1-2*j)*x,H) for j,x in enumerate(g)]
    m=sum(j*p for j,p in enumerate(pi));N=sum(j*(j-1)*p for j,p in enumerate(pi))
    e2=sum(a*b for i,a in enumerate(ds) for b in ds[i+1:])
    assert sum(pi)==1
    assert sum((R+1-2*j)**2*p for j,p in enumerate(pi))==1+sum(d*(d+2) for d in ds)
    assert sum(j*(R+1-j)*p for j,p in enumerate(pi))==F(e2,2)
    assert m>=F(e2,2*(R+1)) and N>=m*m-m and m<=N+1 and m*m<=2*N+1
    n=R+1
    if n>=4:
        c=F(4,n-2) if n%2==0 else F(4*n,(n-1)**2)
        theta=sum(F(g[j])*(j-c*math.comb(j,2)) for j in range(len(g)))
        residual=6*theta/(c*H)-F(3,2)*m
        assert residual>=0
        assert n%2 or residual==0
    if n>=6:
        u=6+F(12,n-4 if n%2==0 else n-3)
        theta=sum(F(g[j])*(math.comb(j,2)-u*F(math.comb(j,3),n-2)) for j in range(len(g)))
        residual=theta/H-N/F(2*(n-4))
        assert residual>=0
        assert n%2 or residual==0
    return True
count=0
for R in range(15):
    for ds in parts(R):
        calc(ds);calc((0,)+ds+(0,));count+=1
result={'product_dimension_range':[0,14],'unordered_nonzero_dimension_lists':count,'zero_factor_checks':count,'all_generic_moment_identities_exact':True,'actual_cusp_substitution_R_equals_n_minus_1':True,'theta2_and_theta3_not_identified':True,'source_files_modified':False}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
