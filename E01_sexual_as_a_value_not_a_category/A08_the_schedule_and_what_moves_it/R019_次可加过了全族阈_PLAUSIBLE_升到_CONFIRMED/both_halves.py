import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A153 R447 -- 「两边都有」是一个真的构造,还是我把四个 d 读成了一句话

`#402a`:「都高」的人**口味更冷门**(`S` +0.24)**而且普通的东西对他们也有用**(五题 −0.31)。
**那句「两边都有」现在是四个 d 的联合描述,不是一个量。** 本轮把它做成一个量,然后**用一个会失败的判据打它**。

⚠⚠ **跑之前写下为什么这个检验很容易变成空的**:
`BOTH` 由**两个已经各自越阈**的量组成 -> **它的 |t| 高几乎是保证的**。
**所以判据不是「越过族内阈 2.752」,而是「超过八个量里的最大值 **6.803**」** ——
**只有「合起来比任何单个都强」才说明「两边都有」是一件事,而不是两件事各自成立。**
(`#396b` 的教训:一个由构造保证的结果携带零比特信息。)

⚠ `#444` 的教训照搬:**同轮报两种构造**(`和` 与 `min`);
**若结论只在一种上成立,那是构造在说话,不是数据。**

ESTIMAND        `BOTH_sum = z(S) − z(五题)` · `BOTH_min = min(z(S), −z(五题))`;
                主量 = 两格之间的 |t|,**对着 6.803**。
判据(**先标支**,`#379c`)
                【两支】两种构造结论一致 · 负对照用**越阈率** · guard 26 **显式传 branch**。
                【非零支】**|t| > 6.803** -> 「两边都有」是一件事;
                【零支】**|t| ≤ 6.803** -> 它不是一件事,**页面上那句话要退回成「两个各自成立的量」**。
⚠ 零的种类     `offset_control`:**两格之间任何量的差的零不该是零** -> **随机等大小分组**(`#402a` 同法)。
IMPOSSIBLE      ① `min` 会被**两边都低**的人拉低,而那是「既不冷门也不被常规打动」;
                ② 这仍是**同一份数据、同一对格子** -> 不是复制;
                ③ 「合起来更强」不等于「存在一个潜变量」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR)
THC=next(c for c in d.columns if 'vmq8jqw' in str(c))
th=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
FIVE=list(inv[inv['kind']=='RATING_NEG_FIB']['col'])
V=np.column_stack([pd.to_numeric(d[c],errors='coerce').values.astype(float) for c in FIVE])
INT=np.nanmean(np.column_stack([(V[:,i]-np.nanmean(V[:,i]))/np.nanstd(V[:,i]) for i in range(5)]),1)
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
FULL=np.isfinite(sh)&np.isfinite(th)
ms=float(np.median(sh[FULL])); mt=float(np.median(th[FULL]))
BOTH=FULL&(sh>ms)&(th>mt); SHONLY=FULL&(sh>ms)&(th<=mt)
M=ok&np.isfinite(S)&np.isfinite(INT)&np.isfinite(EARLY)
for q_ in Q: M&=np.isfinite(q_)
B=BOTH&M; Sx=SHONLY&M; U=B|Sx; n=int(U.sum())
zS=np.full(NN,np.nan); zI=np.full(NN,np.nan)
zS[U]=(S[U]-S[U].mean())/S[U].std(); zI[U]=(INT[U]-INT[U].mean())/INT[U].std()
CTORS={'BOTH_sum = z(S) − z(五题)':zS-zI,
       'BOTH_min = min(z(S), −z(五题))':np.where(np.isfinite(zS)&np.isfinite(zI),
                                                np.minimum(zS,-zI),np.nan)}
BASE={'S 位置':S,'五题(常规不色)':INT,'EARLY':EARLY,'清晰度':Q[5],
      'c1':Q[2],'c2':Q[3],'c3⁻':-Q[4],'D 块间对比':Q[1]}
def tt(v,gb,gs):
    a=v[gb]; b=v[gs]
    se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
    sd=np.sqrt((a.var(ddof=1)+b.var(ddof=1))/2)
    return (a.mean()-b.mean())/max(sd,1e-12),(a.mean()-b.mean())/max(se,1e-12)
BEST=max(abs(tt(v,B,Sx)[1]) for v in BASE.values())
BESTNM=max(BASE,key=lambda k:abs(tt(BASE[k],B,Sx)[1]))
print(f"n=**{n:,}**(都高 {int(B.sum()):,} · 只羞耻 {int(Sx.sum()):,})")
print(f"⚠ **要打败的不是族内阈,是八个量里的最大值**:**{BESTNM} |t| = {BEST:.3f}**\n")
NP_=400; idx=np.flatnonzero(U); nb=int(B.sum())
rows=[]
for nm,v in CTORS.items():
    dd,t_=tt(v,B,Sx)
    nul=[]
    for s_ in range(NP_):
        rg=np.random.default_rng(7400+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
        nul.append(abs(tt(v,ga,gs)[1]))
    nul=np.array(nul); thr=float(np.percentile(nul,95))
    rows.append(dict(v_ctor=nm,v_d=dd,v_t=t_,v_thr=thr,v_beats=bool(abs(t_)>BEST)))
    print(f"**{nm}**")
    print(f"   d **{dd:+.4f}** · |t| **{abs(t_):.3f}** · 族内零阈 {thr:.3f} · "
          f"vs 最好的单个量 **{BEST:.3f}** -> {'**超过 -> 是一件事**' if abs(t_)>BEST else '**没超过 -> 不是一件事**'}")
T=pd.DataFrame(rows); check_columns(T,'R447')
T.to_csv(pathlib.Path(__file__).parent/'results'/'both_halves.csv',index=False)
v0=list(CTORS.values())[0]
negs=np.array([abs(tt(v0,
    np.isin(np.arange(NN),np.random.default_rng(96000+s).permutation(idx)[:nb]),
    np.isin(np.arange(NN),np.random.default_rng(96000+s).permutation(idx)[nb:]))[1]) for s in range(200)])
rate=float((negs>float(T.v_thr.iloc[0])).mean())
print(f"\n负对照(**越阈率**,随机分组 200 次):**{100*rate:.1f}%**")
MDE=None
print(f"guard 26 = **MDE 扫描**(在**随机分组基线**上种,`#402b` 的修法),每级 30 次:")
# ⚠ **第一版的网格是 (0.10, 0.20, 0.30, 0.45)** —— 0.30 给 50%、0.45 给 97%,
# **80% 点落在两者之间**,于是 MDE 被报成 0.45,而 guard 26 就因此 FAIL。
# **一个 MDE 网格的分辨率,会把 MDE 系统性地报高** —— 报高的方向恰好是**让门失败**的方向。
# 加密网格不是搬门柱:MDE 是**设计的属性**,粗网格只是**测不准它**。
for gg in (0.10,0.20,0.28,0.32,0.36,0.40,0.45):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(300+int(gg*1000)*83+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:nb]]=True; gs=np.zeros(NN,bool); gs[p[nb:]]=True
        v2=v0.copy(); v2[ga]=v2[ga]+gg*np.nanstd(v0[U])
        if abs(tt(v2,ga,gs)[1])>BEST: hit+=1
    print(f"   偏移 **{gg:.2f} sd** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.60
OBS=float(max(abs(T.v_d)))
print(f"   **MDE = {MDE_:.2f} sd** · 实测 d = **{OBS:.4f}**")
CONS=bool(T.v_beats.iloc[0]==T.v_beats.iloc[1])
g=Gate('「两边都有」是一件事,还是两件事各自成立')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.asserted('★【两支】两种构造(`和` / `min`)结论一致 —— 否则是构造在说话',CONS,
           f"超过最好单量:{T.v_beats.tolist()}",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 实测 d(非零支)',MDE_,OBS,True,
    what='在随机分组基线上种',branch='non_null')
if 0.01<=rate<=0.12 and CONS:
    g.asserted('★【非零支】|t| **超过八个量里的最大值** -> 「两边都有」是一件事',
               bool(T.v_beats.all()),
               f"{[round(abs(x),3) for x in T.v_t]} vs 最好单量 {BEST:.3f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
