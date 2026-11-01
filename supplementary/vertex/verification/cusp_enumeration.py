"""Enumerate the reduced cusp partitions with exact rational comparison."""
from pathlib import Path
from fractions import Fraction as F
import json,time

OUT=Path(__file__).absolute().parent
data=json.loads((OUT/'budget_results.json').read_text())
psi={d:F(5*d,8) for d in range(6)}
for row in data['affine']:
    if row['d']>=6: psi[row['d']]=max(F(s) for s in row['scores'].values())

def multiply(h,d):
    out=[];acc=0
    for j in range(len(h)+d):
        if j<len(h): acc+=h[j]
        if j>d and j-d-1<len(h): acc-=h[j-d-1]
        out.append(acc)
    return out

def first_num(h,R):
    # Direct g_j=h_j-h_{j-1}, no probability terminology or floating point.
    acc=prev=0
    for j in range(R//2+1):
        gj=h[j]-prev;prev=h[j]
        assert gj>=0
        acc+=j*(R+1-2*j)*gj
    return acc

def independent_num(h,R,H):
    z=(R+1)*H-2*sum(abs(2*j-R)*a for j,a in enumerate(h))-(h[R//2] if R%2==0 else 0)
    assert z%2==0
    return z//2

def run():
    start=time.monotonic()
    targets={n:{'count':0,'num':-1,'den':1,'witness':None} for n in range(82,90)}
    visited=0
    def visit(parts,R,h,H,P,minpart):
        nonlocal visited
        visited+=1
        if R>=77:
            for small in range(5):
                n=R+small+1
                if n not in targets: continue
                hs=h if small==0 else multiply(h,small)
                Hs=H*(small+1)
                num=first_num(hs,R+small)
                ps=P+psi[small]
                A,B=ps.numerator,ps.denominator
                cost_num=2*A*Hs-3*num*B
                cost_den=2*B*Hs
                row=targets[n];row['count']+=1
                if cost_num*row['den']>row['num']*cost_den:
                    assert num==independent_num(hs,R+small,Hs)
                    row.update(num=cost_num,den=cost_den,witness=list(parts)+([small] if small else []),m=str(F(num,Hs)))
        for d in range(minpart,89-R):
            visit(parts+(d,),R+d,multiply(h,d),H*(d+1),P+psi[d],d)
        if visited%100000==0:
            print('Visited',visited,'elapsed',round(time.monotonic()-start,1),flush=True)
    visit((),0,[1],1,F(0),5)
    expected=['105177423/11358856','26064463/2839714','103537961/11358856','51508227/5679428','102692899/11358856','10474701/1169294','808290543/90870848','451668671/51114852']
    counts=[226774,253863,283993,317579,354899,396446,442587,493897]
    rows=[]
    for n,row in targets.items():
        val=F(row.pop('num'),row.pop('den'))
        assert val==F(expected[n-82]),(n,val)
        assert row['count']==counts[n-82],(n,row['count'])
        row.update(n=n,bound=str(val));rows.append(row)
    result={'visited_partial_partitions':visited,'enumerated_reduced_partitions':sum(r['count'] for r in rows),'seconds':time.monotonic()-start,'cusp_table':rows}
    (OUT/'cusp_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':run()
