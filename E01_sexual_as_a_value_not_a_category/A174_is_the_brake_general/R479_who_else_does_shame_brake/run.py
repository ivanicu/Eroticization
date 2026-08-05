import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #434b found shame SUBTRACTS from porn use's path to behaviour. Is that specific to porn
   use, or does shame brake every quantity on this page?

⚠ THE ARITHMETIC TRAP, named before running (ss9): b (shame->acted) is SHARED, so ANY focal
   with a>0 gets a negative indirect and any with a<0 gets a positive one. The SIGN is
   algebraically forced and carries almost no information. Therefore:
   (1) the panel MUST include a focal with a<0, or the round only demonstrates arithmetic --
       `-z(five-item)` has a = -0.0658 on shame, so it is included by design;
   (2) the QUANTITY OF INTEREST is not the sign but the SHARE |a*b / c|: how much of each
       total the brake removes. That is NOT forced -- it depends on how large each direct
       path is relative to its own route through shame.

Worlds
  A  one share       : the brake removes a similar fraction of every quantity -> shame acts
                       like a uniform tax on acting, and the page can say so once.
  B  it discriminates: the shares differ by more than their own resampling spread -> shame
                       brakes some routes to behaviour much harder than others, and WHICH
                       ones is the sentence about people.

Pre-registered: shares compared against a bootstrap spread; KILL = if every share's interval
overlaps every other's, world B is dead and A is reported.
CONTROL : the OLS identity c = c' + a*b must hold for EVERY focal (#434a).
CONTROL2: the panel contains >=1 focal with a<0 (else the round is the trap above).
FRONTIER.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show

GATE=Gate("R479 is the brake general")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
PH=pd.to_numeric(raw['pornhabit'],errors='coerce').values
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float); ACT=np.asarray(OUT['实践了多少'],dtype=float)
MM = M & np.isfinite(PH) & np.isfinite(AGE) & np.isfinite(SHAME) & np.isfinite(ACT)

FOCALS={'色情使用量':PH, '冷门程度 S':A, '常规也管用(−五题)':Bv, '广度型 c3⁻':C3, '类别数':ncat}
def other(nm): return [v for k,v in FOCALS.items() if k!=nm]

def paths(fname, fvec, idx):
    ctrl=[z(v,idx) for v in other(fname)]+[z(AGE,idx)]
    Xb=np.column_stack([np.ones(idx.sum())]+ctrl)
    p=z(fvec,idx); sh=z(SHAME,idx); ac=z(ACT,idx)
    a =float(np.linalg.lstsq(np.column_stack([Xb,p]),   sh,rcond=None)[0][-1])
    c =float(np.linalg.lstsq(np.column_stack([Xb,p]),   ac,rcond=None)[0][-1])
    bb=       np.linalg.lstsq(np.column_stack([Xb,p,sh]),ac,rcond=None)[0]
    return a, float(bb[-1]), c, float(bb[-2])

rng=np.random.default_rng(11); idxall=np.flatnonzero(MM)
rows=[]
for nm,v in FOCALS.items():
    a,b,c,cp=paths(nm,v,MM); ind=a*b
    nul=np.array([ (lambda t:t[0]*t[1])(paths(nm,perm_in(v,MM,seed=9500+i),MM)) for i in range(200)])
    sh_=[]
    for i in range(200):
        take=rng.choice(idxall,len(idxall),replace=True)
        mm=np.zeros(len(MM),bool); mm[np.unique(take)]=True
        t=paths(nm,v,mm)
        if abs(t[2])>1e-9: sh_.append(t[0]*t[1]/t[2])
    sh_=np.array(sh_); lo,hi=np.percentile(sh_,[2.5,97.5])
    rows.append(dict(focal=nm,a=a,b=b,c=c,cp=cp,indirect=ind,
                     ind_thr=float(np.percentile(np.abs(nul),95)),
                     share=ind/c if abs(c)>1e-9 else np.nan,
                     sh_lo=lo,sh_hi=hi, ident=abs(c-cp-ind)))
T=pd.DataFrame(rows)
show(T[['focal','a','b','c','cp','indirect','ident']], HERE/'results/paths.csv', n=8, label="三路径")
print()
show(T[['focal','share','sh_lo','sh_hi','indirect','ind_thr']],
     HERE/'results/shares.csv', n=8, label="份额")

GATE.asserted("CONTROL the OLS identity holds for every focal",
              bool((T.ident<1e-9).all()), f"max residual = {T.ident.max():.2e}", kind="control")
GATE.asserted("CONTROL2 the panel contains a focal with a<0 (else the round is arithmetic)",
              bool((T.a<0).any()), f"a<0 focals = {list(T[T.a<0].focal)}", kind="control")

# do the shares differ by more than their own spread?
iv=[(r.sh_lo,r.sh_hi) for _,r in T.iterrows()]
disjointpairs=sum(1 for i in range(len(iv)) for j in range(i+1,len(iv))
                  if iv[i][1]<iv[j][0] or iv[j][1]<iv[i][0])
npairs=len(iv)*(len(iv)-1)//2
print(f"\n份额区间两两**不相交**的对数 = **{disjointpairs} / {npairs}**")
GATE.asserted("KILL the brake discriminates between quantities",
              disjointpairs>0, f"disjoint share intervals = {disjointpairs}/{npairs}")
verdict = "DISCRIMINATES" if disjointpairs>0 else "UNIFORM_TAX"
print(f"判决 = {verdict}   (预注册预测:符号是代数必然;有信息的是份额)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),disjoint=disjointpairs,npairs=npairs,
               rows=T.to_dict('records')), open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
