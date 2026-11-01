"""Exact envelope reduction and its 341 small endpoint records."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
from cusp_enumeration import multiply,first_num,independent_num,psi
OUT=Path(__file__).absolute().parent

def envelope(d):
    return min(F(5*d,8),F(19*d,40)+F(6,5),F(d,4)+F(24,5),F(13*d,160)+F(327,40),F(81,8))

def run():
    assert all(psi[d]<=envelope(d) for d in range(100))
    records=[]
    best=F(-1);witness=None
    for a,b,c,d in product(range(24),range(8),range(5),range(4)):
        eta=F(a)+F(16*b+26*c,5)+F(63*d,8)
        if eta>=24:continue
        dims=[8]*a+[16]*b+[20]*c+[24]*d
        h=[1];H=1;R=sum(dims)
        for dim in dims:h=multiply(h,dim);H*=dim+1
        num=first_num(h,R)
        assert num==independent_num(h,R,H)
        val=5*a+F(44*b+49*c,5)+F(81*d,8)-F(3*num,2*H)
        if val>best:best,witness=val,(a,b,c,d)
        records.append({'counts':[a,b,c,d],'eta':str(eta),'m':str(F(num,H)),'bound':str(val)})
    assert len(records)==341
    assert best==F(254166,24565) and witness==(0,3,0,0)
    assert F(83,8)-best==F(5567,196520)
    for d in (8,16,20,24):assert d*(d+2)<=90*(F(3*d,4)-envelope(d))
    result={'records':len(records),'max':str(best),'witness':witness,'gap':str(F(83,8)-best),'checks':records}
    (OUT/'envelope_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2),flush=True)

if __name__=='__main__':run()
