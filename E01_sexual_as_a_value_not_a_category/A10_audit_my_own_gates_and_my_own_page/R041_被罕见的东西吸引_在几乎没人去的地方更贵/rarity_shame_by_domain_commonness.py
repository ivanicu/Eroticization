import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: `E02` (another session, same repo, opened on Ivan's counter-examples today) argues at the
   SOCIETY level that condemnation tracks rarity only where norms can actually suppress the
   behaviour -- and that sexual interest is not such a domain. That is D4 and not established.
   #468's NEXT: test the person-level corollary here.

⚠ MY OWN NEXT PRE-REGISTERED THE DIRECTION BACKWARDS, AND THAT IS FIXED BEFORE RUNNING.
   It said world B is "the slope is stronger where the domain is COMMON". Reasoned through:
   norm suppression operates where norms BITE, and where they bite the behaviour becomes RARE.
   So suppression predicts a TIGHTER rarity-shame coupling in UNCOMMON domains, not common
   ones. The signpost had the sign inverted; this is the fifth time this session that reading
   my own NEXT before acting caught an error in it (#452a #459a #462b #464b).

Worlds (directions now stated correctly)
  A  the slope does not depend on how common the domain is -> consistent with `E02`'s D4
     reading, as person-level evidence.
  B  the slope is STRONGER in UNCOMMON domains -> norm suppression leaves a person-level
     trace, and `E02`'s reading would have to narrow.
  C  stronger in COMMON domains -> neither account predicts this; it would need its own.

Design, built to avoid the circularity this question invites: rarity is computed WITHIN a
block (a person's mean surprisal over the options they picked there), and the SPLIT is BETWEEN
blocks by how often that block's options are picked at all. So the thing being split on is not
the thing being measured.
Both group scores enter the SAME model, so the comparison is within-person.
Spec curve over the cut point (#P16). Difference of slopes gets its own bootstrap.
CONTROL : each group's rarity score must correlate with the pooled position score -- if not,
   the per-group construction is not measuring rarity.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show

GATE=Gate("R513 does the coupling depend on commonness")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
# 每块:人的块内稀有度 v_b(勾到的选项的平均意外度)+ 该块的人群勾选率
rate=np.full(NB,np.nan); RAR=np.full((NB,NN),np.nan)
for b,(Mb,ppl) in enumerate(MB):
    rr=-np.log(np.clip(Mb.mean(0),1e-4,1.)); nb=Mb.sum(1)
    rate[b]=float(Mb.mean())
    v=np.where(nb>0,(Mb@rr)/np.maximum(nb,1),np.nan)
    RAR[b,ppl]=v
print(f"块数 **{NB}** · 人群勾选率 min **{np.nanmin(rate):.3f}** · 中位 **{np.nanmedian(rate):.3f}** · "
      f"max **{np.nanmax(rate):.3f}**")
_RAR=RAR.copy(); _rate=rate.copy(); _m=m.copy()

_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SH=np.asarray(OUT['羞耻'],dtype=float)

def group_scores(cut, need=3):
    hi=_rate>=cut; lo=~hi
    def agg(sel):
        sub=_RAR[sel]; cnt=np.isfinite(sub).sum(0)
        return np.where(cnt>=need, np.nanmean(sub,0), np.nan)
    return agg(hi), agg(lo), int(hi.sum()), int(lo.sum())

rows=[]
CUTS=np.quantile(_rate[np.isfinite(_rate)],[0.35,0.45,0.5,0.55,0.65])
for cut in CUTS:
    H,L,nh,nl=group_scores(cut)
    mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L)
    n=int(mm.sum())
    if n<1000: continue
    X=np.column_stack([np.ones(n), z(H,mm), z(L,mm), z(ncat,mm), z(AGE,mm)])
    b=np.linalg.lstsq(X,z(SH,mm),rcond=None)[0]
    rows.append(dict(cut=float(cut), n_hi_blocks=nh, n_lo_blocks=nl, n=n,
                     b_common=float(b[1]), b_uncommon=float(b[2]),
                     diff=float(b[2]-b[1])))
T=pd.DataFrame(rows)
show(T, HERE/'results/spec_curve.csv', n=6, label="切点规格曲线(常见 vs 罕见领域)")

mid=T.iloc[len(T)//2]
H,L,_,_=group_scores(mid['cut'])  # ⚠ #469a:`mid.diff` 取到 pandas 方法,本类第四次
mm = M & np.isfinite(AGE) & np.isfinite(H) & np.isfinite(L)
n=int(mm.sum())
rS=float(np.corrcoef(H[mm],A[mm])[0,1]); rS2=float(np.corrcoef(L[mm],A[mm])[0,1])
print(f"\n**对照**:两组稀有度分与合并位置分 `S` 的相关 = **{rS:+.3f}** / **{rS2:+.3f}**")
GATE.asserted("CONTROL both group scores really measure rarity",
              min(abs(rS),abs(rS2))>0.3, f"corr with S = {rS:+.3f} / {rS2:+.3f}", kind="control")

rg=np.random.default_rng(131); idx=np.flatnonzero(mm); bs=[]
for _ in range(400):
    take=rg.choice(idx,len(idx),replace=True)
    m2=np.zeros(len(mm),bool); m2[np.unique(take)]=True
    k=int(m2.sum())
    X=np.column_stack([np.ones(k), z(H,m2), z(L,m2), z(ncat,m2), z(AGE,m2)])
    bb=np.linalg.lstsq(X,z(SH,m2),rcond=None)[0]
    bs.append(float(bb[2]-bb[1]))
bs=np.array(bs); lo,hi=np.percentile(bs,[2.5,97.5])
print(f"中位切点(勾选率 {mid['cut']:.3f},{int(mid['n_hi_blocks'])} 常见块 / {int(mid['n_lo_blocks'])} 罕见块,n={n:,}):")
print(f"   常见领域斜率 **{mid['b_common']:+.4f}** · 罕见领域斜率 **{mid['b_uncommon']:+.4f}** · "
      f"差(罕见−常见)**{mid['diff']:+.4f}**")
print(f"   自助 95% **[{lo:+.4f}, {hi:+.4f}]** -> "
      f"{'**不含 0:两组不同**' if (lo>0)==(hi>0) else '**含 0:两组分不开**'}")
same = not ((lo>0)==(hi>0))
signs = T['diff'].gt(0).sum()
print(f"   规格曲线:{len(T)} 个切点里差为正(罕见更强)的有 **{signs}**")
GATE.asserted("KILL the slope depends on how common the domain is", not same,
              f"difference {mid['diff']:+.4f}, boot [{lo:+.4f},{hi:+.4f}]")
verdict = ("NO_DEPENDENCE" if same else
           ("STRONGER_IN_UNCOMMON" if mid['diff']>0 else "STRONGER_IN_COMMON"))
print(f"\n判决 = **{verdict}**")
json.dump(dict(verdict=verdict,n=n,rows=T.to_dict('records'),
               diff=float(mid['diff']),boot=[lo,hi],corr_S=[rS,rS2]),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
