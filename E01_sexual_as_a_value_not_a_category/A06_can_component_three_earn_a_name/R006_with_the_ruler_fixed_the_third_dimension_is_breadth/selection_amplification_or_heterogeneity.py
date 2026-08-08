import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A66 R288 -- 那个 4.3 倍,是真实异质还是选择放大

`#242b`:14 岁前就报了 ≥8 个起始年龄的那 2,806 人,未截断的 **Δ = −0.1407**,全样本的 **4.3 倍**。
但那批人是**按「早就有很多兴趣」选出来的**,而 Δ 是 (起始残差 × 稀有度) 的人内相关。

WORLDS          ① **真实异质**:早熟的人身上这个关系确实更强 -> 一条关于发展的事实
                ② **选择放大**:「早期类别数 ≥8」与 Δ 的分子相关,选择本身把它抬起来
ESTIMAND        按**总类别数**(与「早/晚」无关的条件)匹配一个同样大小、同样类别数分布、
                且**不在**早熟组里的子集,算它的未截断 Δ。
KILL            **若匹配子集的 Δ 也接近 −0.14 -> 是「报得多的人 Δ 更强」,与早熟无关,`#242b` 改写;
                若仍停在 −0.03 附近 -> 早熟这一条是真的,`#242b` 升为声明。**
POSITIVE CTRL   两端:① 直接按「早期类别数」切,必须重现 −0.1407;
                ② 随机切同样大小的子集,必须落在 −0.0328 附近。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      「早熟」在这里被操作化为**14 岁前报了 ≥8 个类别**,它同时携带
                「记得早」与「真的早」两种成因(`#114`),本轮不分。
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
rho=rho_of(V0); base=np.isfinite(rho)
NCAT=np.isfinite(V0).sum(1); NEARLY=(np.isfinite(V0)&(V0<=14)).sum(1)
EARLY=base&(NEARLY>=8)
rng=np.random.default_rng(20260804)
def D_of(mask):
    m=mask&base; v=rho[m]
    sd=float(np.std([np.mean(v[i]) for i in (rng.choice(len(v),len(v),True) for _ in range(300))]))
    return float(np.mean(v)),sd,int(m.sum())
d_all=D_of(base); d_early=D_of(EARLY)
print(f"全样本 Δ = {d_all[0]:+.4f} ± {d_all[1]:.4f}(n={d_all[2]:,})")
print(f"早熟组(14 岁前 ≥8 类别)Δ = **{d_early[0]:+.4f}** ± {d_early[1]:.4f}(n={d_early[2]:,})"
      f"  [`#242b` 报 −0.1407]")
print(f"⚠ 两组的总类别数中位:全样本 {int(np.median(NCAT[base]))} · 早熟组 {int(np.median(NCAT[EARLY]))}"
      f"  <- 混杂在这里")

# 按总类别数匹配的、不在早熟组里的子集
pool=base&(~EARLY)
sel=[]
for c in np.unique(NCAT[EARLY]):
    need=int((NCAT[EARLY]==c).sum()); cand=np.flatnonzero(pool&(NCAT==c))
    if len(cand)==0: continue
    sel.append(rng.choice(cand,min(need,len(cand)),replace=False))
MATCH=np.zeros(N,bool); MATCH[np.concatenate(sel)]=True
d_match=D_of(MATCH)
print(f"\n**按总类别数匹配、且不在早熟组的子集:Δ = {d_match[0]:+.4f} ± {d_match[1]:.4f}"
      f"(n={d_match[2]:,};类别数中位 {int(np.median(NCAT[MATCH]))})**")
# 正对照②:随机同样大小
rd=[]
for s in range(5):
    r_=np.zeros(N,bool); r_[np.random.default_rng(700+s).choice(np.flatnonzero(base),d_early[2],replace=False)]=True
    rd.append(D_of(r_)[0])
print(f"正对照②(随机同样大小 n={d_early[2]:,}):Δ = {np.mean(rd):+.4f} ± {np.std(rd):.4f}"
      f"  [必须落在 {d_all[0]:+.4f} 附近]")
# 正对照①:换一个早熟阈值,必须仍在 −0.14 附近
E15=base&((np.isfinite(V0)&(V0<=15)).sum(1)>=8); d15=D_of(E15)
print(f"正对照①(阈值换成 15 岁 ≥8):Δ = {d15[0]:+.4f} ± {d15[1]:.4f}(n={d15[2]:,})")
# ⚠ 未注册但必须报的混杂:两组的人内起始年龄离散度不同 -> 范围受限
sd_on=lambda m: float(np.nanmean(np.nanstd(np.where(np.isfinite(V0)&m[:,None],V0,np.nan),axis=1)))
sd_rar=lambda m: float(np.nanmean(np.nanstd(np.where(np.isfinite(V0)&m[:,None],rar0[None,:],np.nan),axis=1)))
print(f"\n⚠ 范围受限体检:人内起始年龄 sd —— 早熟组 {sd_on(EARLY):.3f} · 匹配组 {sd_on(MATCH):.3f} · "
      f"全样本 {sd_on(base):.3f}")
print(f"                人内稀有度 sd —— 早熟组 {sd_rar(EARLY):.3f} · 匹配组 {sd_rar(MATCH):.3f} · "
      f"全样本 {sd_rar(base):.3f}")
print(f"  (若早熟组的离散度**更小**却给出**更强**的相关,范围受限解释不掉它 —— 受限只会衰减)")

nul=[]
for _ in range(20):
    Vp=V0.copy()
    for j in range(Mc):
        idx=np.flatnonzero(np.isfinite(Vp[:,j])); Vp[idx,j]=Vp[rng.permutation(idx),j]
    r2=rho_of(Vp); nul.append(float(np.nanmean(r2[np.isfinite(r2)&EARLY])))
print(f"置换零(早熟组上){np.mean(nul):+.4f} ± {np.std(nul):.4f}")

T=pd.DataFrame([dict(arm='全样本',d=d_all[0],sd=d_all[1],n=d_all[2],med_ncat=int(np.median(NCAT[base]))),
                dict(arm='早熟组',d=d_early[0],sd=d_early[1],n=d_early[2],med_ncat=int(np.median(NCAT[EARLY]))),
                dict(arm='类别数匹配_非早熟',d=d_match[0],sd=d_match[1],n=d_match[2],
                     med_ncat=int(np.median(NCAT[MATCH])))])
check_columns(T,'R288'); T.to_csv(pathlib.Path(__file__).parent/'results'/'selection_test.csv',index=False)

g=Gate('4.3 倍:真实异质还是选择放大')
g.asserted('正对照①:换一个早熟阈值(15 岁),Δ 必须仍在 −0.14 附近',
           d15[0]<-0.10, f"15 岁阈值 {d15[0]:+.4f}(n={d15[2]:,})vs 14 岁 {d_early[0]:+.4f}")
g.asserted('正对照②:随机同样大小的子集必须落在全样本值附近',
           abs(np.mean(rd)-d_all[0])<0.01, f"随机 {np.mean(rd):+.4f} vs 全样本 {d_all[0]:+.4f}")
g.negative_control('置换零(早熟组)',abs(float(np.mean(nul))),abs(d_early[0]),
                   null_spread=float(np.std(nul)),null_kind='题内跨人置换起始年龄 —— 只打掉配对')
g.offset_control('★ 类别数匹配的非早熟子集 vs 早熟组',d_match[0],d_early[0],
                 float(np.hypot(d_match[1],d_early[1])),
                 null_kind='早熟组的 Δ —— 不是零假设,是「若这只是报得多,匹配子集该落在哪」')
g.asserted('⚠ 未注册但必须报:范围受限体检(受限只会衰减,不会造出更强的相关)',
           sd_on(EARLY)<sd_on(MATCH), f"人内起始年龄 sd 早熟 {sd_on(EARLY):.3f} < 匹配 {sd_on(MATCH):.3f}")
g.asserted('★ 注册的 kill:匹配子集仍停在全样本值附近 -> 早熟是真的;接近 −0.14 -> 只是报得多',
           abs(d_match[0]-d_all[0])<abs(d_match[0]-d_early[0]),
           f"匹配 {d_match[0]:+.4f};全样本 {d_all[0]:+.4f};早熟 {d_early[0]:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
