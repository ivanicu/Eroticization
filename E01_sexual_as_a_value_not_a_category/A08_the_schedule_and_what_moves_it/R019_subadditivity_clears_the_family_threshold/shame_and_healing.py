import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A153 R445 -- 有没有一种人,既深深羞耻,又觉得这件事是疗愈的

页面把羞耻与「治疗性」当作**两个独立的结局**在处理(四条路 vs 两条路),
**但从没量过它们的联合分布。**
而这是一个**关于人的**问题,不是关于变量的:**那种人存在吗?有多少?**

两个活着的世界:
**A 一条轴** —— 相关强负,「两者都高」那一格**远低于**独立期望 -> 羞耻与疗愈是同一条轴的两端;
**B 两个维度** —— 相关弱,「都高」那一格**接近或高于**独立期望 -> **那一格里的人是真实存在的一群**。

⚠ **这一轮不控制任何东西。** 问的是**这些人存不存在**,不是「控制之后还存不存在」——
控制会把一群真实的人变成一个残差,而残差里没有人。

ESTIMAND        ① `corr(羞耻, 治疗性)`;② 按两者**中位数**分四格,报四格**人数占比**;
                主量 = **「都高」那一格的占比**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 用 **MDE 扫描**;`#392e` 先看两个结局各自。
                【非零支】「都高」占比**越过**打断配对的零 -> 世界 B(这群人比偶然更多);
                          **低于**零的下侧 -> 世界 A(被压制)。
⚠ 零的种类     `offset_control`:**「都高」占比的零绝不是零** ——
                两者的**边际分布**决定了即使完全独立也会有一群人落在那一格。
                零 = **边际不变、配对被打断**(`lib.nulls.perm_in` 打乱其中一个),取该格占比的分布。
IMPOSSIBLE      ① 中位数分割丢掉强度信息 -> **同轮报相关(连续)与四格(离散)两者**;
                ② 两题都自报,同一份问卷 -> 同源方差会**抬高**相关,即**偏向世界 A**,方向上保守;
                ③ 「存在这群人」不等于「这群人是一个类型」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
SHC=next(c for c in d.columns if 'ashamed' in str(c))
THC=next(c for c in d.columns if 'vmq8jqw' in str(c))
sh=pd.to_numeric(d[SHC],errors='coerce').values.astype(float)
th=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
print("⚠ **`#392e`:两个结局各自先看清楚**")
for nm,v in (('羞耻',sh),('治疗性',th)):
    g=np.isfinite(v)&np.isfinite(anc)
    print(f"   {nm}:取值 {np.unique(v[np.isfinite(v)]).tolist()} · "
          f"众数 **{float(pd.Series(v[np.isfinite(v)]).mode().iloc[0]):g}** · n={int(np.isfinite(v).sum()):,} · "
          f"与锚 `Totalsexacts` 相关 **{np.corrcoef(v[g],anc[g])[0,1]:+.4f}**")
M=np.isfinite(sh)&np.isfinite(th)
n=int(M.sum())
R=float(np.corrcoef(sh[M],th[M])[0,1])
print(f"\n① **`corr(羞耻, 治疗性)` = {R:+.4f}**(n={n:,})· "
      f"{'**强负 -> 偏向世界 A**' if R<-0.4 else '**弱 -> 偏向世界 B**'}")
ms=float(np.median(sh[M])); mt=float(np.median(th[M]))
def cells(a,b):
    hi_a=a>ms; hi_b=b>mt
    return dict(both_hi=float((hi_a&hi_b).mean()),sh_only=float((hi_a&~hi_b).mean()),
                th_only=float((~hi_a&hi_b).mean()),both_lo=float((~hi_a&~hi_b).mean()))
C=cells(sh[M],th[M])
print(f"\n② 四格(中位数分割:羞耻 >{ms:g} · 治疗性 >{mt:g}):")
print(f"   **都高 {100*C['both_hi']:.2f}%** · 只羞耻 {100*C['sh_only']:.2f}% · "
      f"只疗愈 {100*C['th_only']:.2f}% · 都低 {100*C['both_lo']:.2f}%")
print(f"   人数:**都高 {int(C['both_hi']*n):,} 人**")
NP_=1000
nul=np.array([cells(sh[M],perm_in(th,M,4400+s)[M])['both_hi'] for s in range(NP_)])
LO=float(np.percentile(nul,2.5)); HI=float(np.percentile(nul,97.5))
print(f"\n⚠ offset 零(**边际分布不变、配对被打断**;"
      f"**这个零绝不是零 —— 即使完全独立也会有一群人落在那一格**):")
print(f"   **{100*nul.mean():.2f}% ± {100*nul.std():.2f}%** · 95% 区间 **[{100*LO:.2f}%, {100*HI:.2f}%]**")
sd_=(C['both_hi']-nul.mean())/max(nul.std(),1e-12)
print(f"   实测 **{100*C['both_hi']:.2f}%** -> **{sd_:+.2f} sd** · "
      f"{'**高于零 -> 世界 B**' if C['both_hi']>HI else ('**低于零 -> 世界 A(被压制)**' if C['both_hi']<LO else '**落在零里 -> 与独立无异**')}")
negs=np.array([cells(sh[M],perm_in(th,M,90000+s)[M])['both_hi'] for s in range(200)])
rate=float(((negs<LO)|(negs>HI)).mean())
print(f"\n负对照(**越界率**,打断配对 200 次):**{100*rate:.1f}%**(合格 1–12%)")
print(f"\nguard 26 = **MDE 扫描**,每级 30 次(在配对上种一个真实的压制):")
MDE=None
for gg in (0.02,0.04,0.06,0.10):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(700+int(gg*100)*73+s_)
        t2=th.copy(); j=np.flatnonzero(M&(sh>ms)&(th>mt))
        k=int(gg*len(j)); 
        if k: t2[rg.choice(j,k,replace=False)]=mt-1        # 把一部分「都高」的人推下去
        c2=cells(sh[M],t2[M])['both_hi']
        if c2<LO or c2>HI: hit+=1
    print(f"   压制 **{gg:.0%}** 的「都高」-> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.15
print(f"   **MDE = {MDE_:.0%}**(相对压制)")
pd.DataFrame([dict(v_corr=R,v_both_hi=C['both_hi'],v_n_both=int(C['both_hi']*n),
                   v_null_mean=float(nul.mean()),v_lo=LO,v_hi=HI,v_sd=sd_,v_n=n,
                   v_sh_only=C['sh_only'],v_th_only=C['th_only'],v_both_lo=C['both_lo'])]).to_csv(
    pathlib.Path(__file__).parent/'results'/'joint.csv',index=False)
g=Gate('有没有一种人,既深深羞耻,又觉得这件事是疗愈的')
g.asserted('★【两支】负对照:**越界率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描(相对压制)vs 实测偏离零的相对幅度',
    MDE_,max(abs(C['both_hi']-nul.mean())/max(nul.mean(),1e-9),1e-9),True,what='MDE 扫描 80% 检出')
g.asserted('★【两支】offset 零非退化(边际决定了独立时也有人落在那一格)',nul.std()>0,
           f"{100*nul.mean():.2f}% ± {100*nul.std():.2f}%",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】「都高」占比**高于**打断配对的零 -> 世界 B(这群人比偶然更多)',
               C['both_hi']>HI,
               f"{100*C['both_hi']:.2f}% vs 零区间 [{100*LO:.2f}%, {100*HI:.2f}%] · {sd_:+.2f} sd")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
