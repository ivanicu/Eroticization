import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A65 R287 -- 清 `#211`:三个受限档的零,有多少是截断,有多少是换了人

**类型:CLOSURE**(如实标注)。`#241b` 的待办表第一条。

`#211b` 报:Δ 在 cap ∈ {14,16,18,∞} 上是 +0.0055 / +0.0012 / −0.0033 / −0.0328,
「三个受限档全部不可分辨」。**但每个 cap 都改变了 KEEP(≥8 个起始年龄),
所以四行是四批不同的人** —— 这正是 `#239a` 的形状,和 `#210` 一模一样。

ESTIMAND        每个 cap 补上第三格:**未截断的 Δ,在该 cap 的那批人上**。
                截断部分 = Δ_capped − Δ_uncapped|同批人;样本部分 = Δ_uncapped|同批人 − Δ_∞。
KILL            **若每一档的「同一批人未截断 Δ」都仍明显为负(接近 −0.0328)->
                受限档的零确实由截断造成,`#211` 成立;
                若它们本身就已接近零 -> 那些零大部分是换了人造成的,`#211` 的作用域要改写。**
POSITIVE CTRL   ① ∞ 档的两个读数必须**逐位相同**(同一批人,同一个量);
                ② 四个已发表的 Δ 必须复现。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      「这批人更早熟」与「这批人 rho 估得更准」本轮不分(`#285` 同款登记)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_of(Vm):
    D=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=8)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
rho_inf=rho_of(V0)
rng=np.random.default_rng(20260804)
rows=[]
print(f"{'cap':>5}{'KEEP':>8}{'Δ(截断)':>12}{'Δ(同批人未截断)':>18}{'Δ_∞ 全样本':>13}"
      f"{'截断部分':>10}{'样本部分':>10}")
for cap in (14,16,18,np.inf):
    Vc=np.where(V0<=cap,V0,np.nan) if np.isfinite(cap) else V0
    r=rho_of(Vc); m=np.isfinite(r)
    dc=float(np.nanmean(r[m])); du=float(np.nanmean(rho_inf[m&np.isfinite(rho_inf)]))
    d_inf=float(np.nanmean(rho_inf[np.isfinite(rho_inf)]))
    sd=float(np.std([np.nanmean(r[np.flatnonzero(m)[i]]) for i in
        (rng.choice(int(m.sum()),int(m.sum()),True) for _ in range(200))]))
    sdu=float(np.std([np.nanmean(rho_inf[np.flatnonzero(m&np.isfinite(rho_inf))[i]]) for i in
        (rng.choice(int((m&np.isfinite(rho_inf)).sum()),int((m&np.isfinite(rho_inf)).sum()),True)
         for _ in range(200))]))
    rows.append(dict(cap=(14 if cap==14 else 16 if cap==16 else 18 if cap==18 else 999),
                     keep=int(m.sum()),d_capped=dc,d_capped_sd=sd,d_uncapped_common=du,
                     d_uncapped_sd=sdu,d_full=d_inf))
    print(f"{('≤%g'%cap if np.isfinite(cap) else '∞'):>5}{int(m.sum()):>8,}"
          f"{dc:>+10.4f}±{sd:.4f}{du:>+14.4f}±{sdu:.4f}{d_inf:>+13.4f}"
          f"{dc-du:>+10.4f}{du-d_inf:>+10.4f}")
T=pd.DataFrame(rows); check_columns(T,'R287')
T.to_csv(pathlib.Path(__file__).parent/'results'/'cap_decomposition.csv',index=False)
pub={14:0.0055,16:0.0012,18:-0.0033,999:-0.0328}
rep=all(abs(float(r.d_capped)-pub[int(r.cap)])<0.004 for _,r in T.iterrows())
print(f"\n复现 `#211b` 的四个已发表 Δ:{'✅' if rep else '❌'} "
      + ' · '.join(f"{('≤%d'%int(r.cap)) if r.cap<999 else '∞'} {r.d_capped:+.4f} vs {pub[int(r.cap)]:+.4f}"
                   for _,r in T.iterrows()))
inf=T[T.cap==999].iloc[0]
print(f"∞ 档两读数逐位相同:{'✅' if abs(inf.d_capped-inf.d_uncapped_common)<1e-12 else '❌'}"
      f"({inf.d_capped:+.6f} vs {inf.d_uncapped_common:+.6f})")
nul=[]
for _ in range(20):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    nul.append(float(np.nanmean(rho_of(Vp))))
print(f"置换零 {np.mean(nul):+.4f} ± {np.std(nul):.4f}")

g=Gate('清 `#211`:截断 vs 换了人')
g.asserted('正对照①:∞ 档的两个读数必须逐位相同',
           abs(inf.d_capped-inf.d_uncapped_common)<1e-12,
           f"{inf.d_capped:+.6f} vs {inf.d_uncapped_common:+.6f}")
g.asserted('正对照②:复现 `#211b` 的四个已发表 Δ',rep,
           ' · '.join(f"{r.d_capped:+.4f}" for _,r in T.iterrows())+" vs +0.0055/+0.0012/−0.0033/−0.0328")
g.negative_control('人内跨人置换',abs(float(np.mean(nul))),abs(float(inf.d_capped)),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
for _,r in T[T.cap<999].iterrows():
    g.control_kept_the_sample(f"cap ≤{int(r.cap)} 的截断控制",
        before=float(inf.d_capped),after=float(r.d_capped),
        n_before=int(inf.keep),n_after=int(r.keep),
        before_common=float(r.d_uncapped_common),after_common=float(r.d_capped),n_common=int(r.keep))
g.asserted('★ 注册的 kill:每档「同批人未截断 Δ」仍明显为负 -> `#211` 成立',
           all(float(r.d_uncapped_common)<-2*float(r.d_uncapped_sd) for _,r in T[T.cap<999].iterrows()),
           ' · '.join(f"≤{int(r.cap)}: {r.d_uncapped_common:+.4f}±{r.d_uncapped_sd:.4f}"
                      for _,r in T[T.cap<999].iterrows()))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
