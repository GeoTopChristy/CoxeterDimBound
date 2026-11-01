from pathlib import Path
from fractions import Fraction as F
from math import comb, prod
from datetime import datetime, timezone
import json

OUT = Path(__file__).absolute().parents[1]

def multiply(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out

def partitions(total, smallest=1):
    if not total:
        yield ()
    for a in range(smallest, total+1):
        for tail in partitions(total-a, a):
            yield (a,) + tail

def choose(a, b):
    return comb(a, b) if a >= b else 0

for a in range(13):
    for b in range(13):
        p = multiply([1]*(a+1), [1]*(b+1))
        q = [0]*(a+b+1)
        for j in range(min(a,b)+1):
            for k in range(j, a+b-j+1):
                q[k] += 1
        assert p == q

count = 0
odd_strict = 0
examples = {}
for degree in range(5, 20):
    n = degree+1
    u = 6 + F(12, n-4 if n%2 == 0 else n-3)
    for dims in partitions(degree):
        h = [1]
        for a in dims:
            h = multiply(h, [1]*(a+1))
        g = [1] + [h[j]-h[j-1] for j in range(1, degree//2+1)]
        assert min(g) >= 0
        block_h = [0]*n
        for j, coefficient in enumerate(g):
            for k in range(j, n-j):
                block_h[k] += coefficient
        assert block_h == h
        vertex_count = prod(a+1 for a in dims)
        pi = [F((n-2*j)*gj, vertex_count) for j, gj in enumerate(g)]
        assert sum(pi) == 1
        m = sum(j*p for j,p in enumerate(pi))
        N = sum(j*(j-1)*p for j,p in enumerate(pi))
        e2 = sum(a*b for i,a in enumerate(dims) for b in dims[i+1:])
        assert sum(j*(n-j)*p for j,p in enumerate(pi)) == F(e2,2)
        assert m >= F(e2,2*n)
        assert N >= m*m-m
        assert m <= N+1 and m*m <= 2*N+1
        centered = sum(F(h[k], vertex_count)*(F(k)-F(degree,2))**2 for k in range(n))
        factor_centered = sum(F(a*(a+2),12) for a in dims)
        block_centered = sum(p*F((n-2*j)**2-1,12) for j,p in enumerate(pi))
        assert centered == factor_centered == block_centered
        theta = sum(gj*(choose(j,2)-u/F(n-2)*choose(j,3)) for j,gj in enumerate(g))
        weighted_theta = F(0)
        for j,p in enumerate(pi):
            kernel = (choose(j,2)-u/F(n-2)*choose(j,3))/F(n-2*j)
            expected = F(j*(j-1), 2*(n-4))
            if n%2:
                expected += F(2*j*(j-1)*(j-2), (n-2*j)*(n-4)*(n-3)*(n-2))
            assert kernel == expected
            weighted_theta += p*kernel
        assert weighted_theta == theta/vertex_count
        assert theta/vertex_count >= N/F(2*(n-4))
        if n%2 == 0:
            assert theta/vertex_count == N/F(2*(n-4))
        elif theta/vertex_count > N/F(2*(n-4)):
            odd_strict += 1
        # Check the unchanged application using the newly explicit finite sums.
        t = degree-max(dims)
        pair_mins = sum(min(a,b) for i,a in enumerate(dims) for b in dims[i+1:])
        assert pair_mins <= F(t*(t+1),2)
        assert e2 >= F(degree*t,2)
        A = F(12*n*(9*n-1)*(n-4),5*(n-1)**2)
        B = F(6*n*(5*n-1),5*(n-1)**2)
        assert F(3,5)*pair_mins <= A*theta/vertex_count+B
        if dims in [(2,3), (1,)*6]:
            examples[str(dims)] = {'n':n,'f0':vertex_count,'g':g,'pi':[str(x) for x in pi],
                                    'm':str(m),'N':str(N),'e2':e2,'Theta':str(theta)}
        count += 1

report = {'status':'passed','utc':datetime.now(timezone.utc).isoformat(),
          'scope':'Exact finite-sum identities and inequalities for the finite weighted sums in Section 5 and its existing local application; no full-suite rerun.',
          'two_factor_identity_cases':169,
          'all_factor_dimension_partitions':count,
          'ambient_dimensions':'6 through 20',
          'odd_cases_with_strict_Theta_inequality':odd_strict,
          'checks':['interval expansion and nonnegative g', 'normalization by f_0(L)',
                    'product and interval centered-square sums', 'factor-dimension weighted identity',
                    'weighted first and second sums', 'exact Theta coefficients for both parities',
                    'local bad-face comparison coefficients'],
          'examples':examples}
(OUT/'logs/lemma45_weighted_sum_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
