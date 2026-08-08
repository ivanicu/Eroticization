import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #472c noted that MODERATION is a new evidence layer and has been used exactly once. The page
   has a second, already-established path -- `c1+` (engagement breadth) reaching the sense of
   being able to change. Is THAT path moderated by how sparsely a domain is entered?

Worlds
  A  not moderated -> "the moderator belongs to shame" becomes the stronger and simpler
     "the moderator belongs to shame and to no path at all".
  B  moderated -> each path has its own moderator, and the fork holds at the moderation layer
     too, which is a stronger statement than #472c.

PRE-REGISTERED PREDICTION: **A**. (Eleven predictions this session, five right.)
⚠ MANDATORY PRE-CHECK (#472's NEXT wrote it): `c1` is a block-level coordinate and the split is
   also block-level. If `c1`'s loadings track the blocks' pick rates, then splitting by pick
   rate would be moderating a quantity BY ITSELF. corr(loading, rate) is reported BEFORE
   anything else, and the round stops if it is high.
Design is the same shape as #469b/#516: build the person's `c1+` score from the high-rate
   blocks only and from the low-rate blocks only, put BOTH in one model, sweep equal k.
CONTROL : the two group scores must each correlate with the full-sample `c1+`.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R517 is the other path moderated")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
keep=list(range(NB)); Ra,Rb=prof_(A,keep),prof_(B,keep)
C=np.zeros((NB,NB))
for i in keep:
    for j in keep:
        g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        if g.sum()>200: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
Rr=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
Rr=np.where(np.isfinite(Ra)|np.isfinite(Rb),Rr,np.nan)
Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
rate=np.array([float(Mb.mean()) for Mb,_ in MB])
v1=V[:,0].copy()
def score_from(vec, blocks):
    sel=np.zeros(NB,bool); sel[list(blocks)]=True
    F=Fm&sel[:,None]; Z=np.where(F,Rr,0.0)
    nu=(vec[:,None]*Z).sum(0); de=(F*np.abs(vec)[:,None]).sum(0)
    return np.where(de>1e-9, nu/np.maximum(de,1e-9), np.nan)
raw0=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
FET0=pd.to_numeric(raw0['totalfetishcategory'],errors='coerce').values.astype(float)
s_full=score_from(v1,range(NB))
gg=m&np.isfinite(s_full)&np.isfinite(FET0)
if np.corrcoef(s_full[gg],FET0[gg])[0,1]<0: v1=-v1          # 朝「卷入更广」(#464c)
s_full=score_from(v1,range(NB))
r_pre=float(np.corrcoef(v1,rate)[0,1])
print(f"⚠ **强制前置**:corr(`c1` 块载荷, 块人群勾选率) = **{r_pre:+.4f}**")
GATE.asserted("PRE-CHECK the coordinate is not the same thing as the split",
              abs(r_pre)<0.6, f"corr(loading, rate) = {r_pre:+.4f}", kind="control")
_v1=v1.copy(); _rate=rate.copy(); _sf=s_full.copy()
_score=lambda blocks: score_from(_v1, blocks)

_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
BE=np.asarray(OUT['能不能改'],dtype=float); SH=np.asarray(OUT['羞耻'],dtype=float)
order=np.argsort(_rate)

def run(y,k,seed=0,nboot=300):
    H=_score(order[-k:]); L=_score(order[:k])
    mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L); n=int(mm.sum())
    if n<800: return None
    X=np.column_stack([np.ones(n), z(H,mm), z(L,mm), z(ncat,mm), z(AGE,mm)])
    b=np.linalg.lstsq(X,z(y,mm),rcond=None)[0]
    rg=np.random.default_rng(seed+k); idx=np.flatnonzero(mm); bs=[]
    for _ in range(nboot):
        take=rg.choice(idx,len(idx),replace=True)
        m2=np.zeros(len(mm),bool); m2[np.unique(take)]=True; kk=int(m2.sum())
        X2=np.column_stack([np.ones(kk), z(H,m2), z(L,m2), z(ncat,m2), z(AGE,m2)])
        bb=np.linalg.lstsq(X2,z(y,m2),rcond=None)[0]; bs.append(float(bb[2]-bb[1]))
    bs=np.array(bs); lo,hi=np.percentile(bs,[2.5,97.5])
    return dict(k=k,n=n,b_common=float(b[1]),b_uncommon=float(b[2]),
                diff=float(b[2]-b[1]),lo=lo,hi=hi,sig=bool((lo>0)==(hi>0)))

# CONTROL:两组分数各自与全量 c1⁺ 相关
H14=_score(order[-14:]); L14=_score(order[:14])
gh=M&np.isfinite(H14)&np.isfinite(_sf); gl=M&np.isfinite(L14)&np.isfinite(_sf)
rh=float(np.corrcoef(H14[gh],_sf[gh])[0,1]); rl=float(np.corrcoef(L14[gl],_sf[gl])[0,1])
GATE.asserted("CONTROL both group scores track the full-sample c1+",
              min(abs(rh),abs(rl))>0.3, f"corr with full c1+ = {rh:+.3f} / {rl:+.3f}", kind="control")
print(f"**对照**:两组分数与全量 `c1⁺` 的相关 = **{rh:+.3f} / {rl:+.3f}**")

rows=[]
for nm,y,sd in (('**能不能改(本轮)**',BE,3100),('羞耻(参照)',SH,3200)):
    for k in (8,10,12,14):
        r=run(y,k,seed=sd)
        if r: rows.append(dict(outcome=nm,**r))
T=pd.DataFrame(rows)
show(T[['outcome','k','n','b_common','b_uncommon','diff','lo','hi','sig']],
     HERE/'results/c1_moderation.csv', n=10, label="`c1⁺` 两组 × 两结局")
be=T[T['outcome'].str.contains('能不能改')]
pos=int((be['diff']>0).sum()); sig=int(be['sig'].sum())
print(f"\n`c1⁺` -> 能不能改:差为正 **{pos}/{len(be)}** · 不含零 **{sig}/{len(be)}** · "
      f"中位差 **{be['diff'].median():+.4f}** · k=14 分辨率 **±{(be[be['k']==14]['hi'].iloc[0]-be[be['k']==14]['lo'].iloc[0])/2:.4f}**")
moderated = sig>=len(be)//2
GATE.asserted("KILL the other path is moderated too (world B)", moderated,
              f"{pos}/{len(be)} positive, {sig}/{len(be)} significant")
verdict = "ALSO_MODERATED" if moderated else "MODERATOR_IS_SHAME_ONLY"
print(f"\n判决 = **{verdict}**  (预注册预测 = A / MODERATOR_IS_SHAME_ONLY)")
json.dump(dict(verdict=verdict,r_precheck=r_pre,rows=T.to_dict('records'),
               ctrl=[rh,rl],prediction="A"),open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
