"""E03·A21·R114 —— 命名不是测量:「作答一致性特质」是真的,还是我给一次巧合起了名字

**类型:FRONTIER**。`#671` 用「作答一致性特质」统一了两次失败(`#663` 与 `#671`),**而命名不是测量。**

⚠ **BASIN**:**那个统一是我上一轮起的名字,而我喜欢它** ⇒ **下注 W2/W3(它被推翻)。**
W1 统一特质(≥3 个内容无关的题组同号且全距 < 0.5)· W2 不统一(全距 ≥ 0.5)⇒ **`#671` 的命名被推翻** ·
**W3 = meta-separator:打乱层序也能得到同样的同号率** ⇒ **「同号」从来就不是证据。**

## ⚠ 两处标注的偏离(跑之前写下)
**①** 预注册写「四档、≥3 题」——**十组候选里只有一组严格四档,照字面执行只剩 1 组,判据无法评估**
   (又一次「不可满足的判据」,与 `#651` `#666` `#670` 同族)。
   **而估计量本来就是天花板归一的(`#647`),跨格式正是它要处理的事** ⇒ 放宽为
   **「同格式、≥3 题、且验证过成块」**。
**②** **言论/任教/藏书三组容忍题内容几乎相同 ⇒ 计为一族**,不是三组独立;
   `#671` 的判据要求「内容无关」,而这一条由我判,**写在这里可被推翻**。

## G1 ESTIMAND
每组:**总体最弱一环**(天花板归一)必须 > 0 才纳入;
然后 **`ρ(年龄层序, 该层最弱一环)`**。**主量 = 纳入组的 ρ 的中位与全距。**
## G2 CONTROLS
**正对照**:性四题 ρ ≈ +0.80、性别角色三题 ρ ≈ +1.00(`#671` 实测)必须复现。
**零(⑤ 写死)**:**把年龄层序随机打乱**,重算全部 ρ 的**同号率**;
  **这个零该不该是零?** 该 —— 打乱之后各组之间不该系统同号 ⇒ `negative_control`。
## G3/G4:十组全报(含被剔除的)· 年龄 {四档, 十分位} 两条规格。
## KILL(条件式)
if 正对照复现 and 打乱后的同号率 < 观测同号率:
  **≥3 个内容无关的族同号 且 ρ 全距 < 0.5** -> **W1 特质被正面支持**
  **全距 ≥ 0.5** -> **W2:不是统一特质,`#671` 的命名被推翻**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**「作答一致性」没有被直接测量的量表**(GSS 无 acquiescence 指标)⇒ 本轮只能测它的**可观察蕴含** ·
族的划分由我判 · **只有 GSS 一具仪器** —— 可证伪形式见 `#671` · `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; FLOOR=200
FAM={"性道德":(["premarsx","xmarsex","homosex","teensex"],"性",None),
     "性别角色":(["fefam","fechld","fepresch"],"性别","fechld"),
     "容忍·言论":(["spkath","spkrac","spkcom","spkmil","spkhomo"],"容忍",None),
     "容忍·任教":(["colath","colrac","colcom","colmil","colhomo"],"容忍",None),
     "容忍·藏书":(["libath","librac","libcom","libmil","libhomo"],"容忍",None),
     "政府该管":(["helppoor","helpnot","helpblk","helpsick"],"再分配",None),
     "机构信任":(["confinan","conbus","coneduc","conpress","conmedic","conjudge"],"信任",None),
     "支出":(["natheal","nateduc","natcrime","natenvir","natrace"],"支出",None),
     "自杀":(["suicide1","suicide2","suicide3","suicide4"],"安乐死",None),
     "堕胎":(["abdefect","abnomore","abhlth","abpoor","abrape","absingle","abany"],"堕胎",None)}
cols=sorted({c for v in FAM.values() for c in v[0]})
_,m=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",metadataonly=True,encoding="latin1")
have=[c for c in cols if c in m.column_names]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year","age"]+have,encoding="latin1")
d=df[df.year>=1990].copy()
for nm,(items,_f,rev) in FAM.items():
    if rev and rev in d.columns: d[rev]=-d[rev]
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if s<0: y=y[::-1]
    return sp(x,y)
def weakest(fr,items):
    vals=[]
    for a,b in combinations(items,2):
        per=[]
        for _,g in fr.groupby("year"):
            mm=g[[a,b]].dropna()
            if len(mm)<FLOOR or mm[a].nunique()<2 or mm[b].nunique()<2: continue
            r=sp(mm[a],mm[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(mm[a],mm[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
        if per: vals.append(float(np.median(per)))
    return (min(vals) if vals else np.nan), len(vals)
print("=== 硬规则①:十组全报,只有最弱一环 > 0 的才纳入 ===")
KEEP={}
for nm,(items,fam,_r) in FAM.items():
    its=[c for c in items if c in d.columns]
    if len(its)<3: print(f"  {nm:8s} 可用 {len(its)} <3 ⇒ 剔除"); continue
    w,np_=weakest(d,its)
    ok=np.isfinite(w) and w>0
    print(f"  {nm:8s}({fam}) {len(its)} 题 · 可算对 {np_} · **总体最弱一环 {w:+.4f}** ⇒ {'✅ 纳入' if ok else '⛔ 不成块,剔除'}")
    if ok: KEEP[nm]=(its,fam)
def grad(items,q=None):
    lab = pd.qcut(d.age,q,labels=False,duplicates="drop") if q else pd.cut(d.age,[17,34,49,64,99],labels=False)
    xs,ys=[],[]
    for i,lv in enumerate(sorted(pd.Series(lab).dropna().unique())):
        w,_=weakest(d[lab==lv],items)
        if np.isfinite(w): xs.append(i); ys.append(w)
    return (sp(xs,ys) if len(xs)>=3 else np.nan), list(zip(xs,ys))
print(f"\n=== G3/G4:纳入 {len(KEEP)} 组 × 年龄 {{四档, 十分位}} ===")
R={}
for nm,(its,fam) in KEEP.items():
    a,_=grad(its); b,_=grad(its,10)
    R[nm]=dict(fam=fam,rho4=float(a) if np.isfinite(a) else None,rho10=float(b) if np.isfinite(b) else None)
    print(f"  {nm:8s}({fam:4s}) 四档 ρ = **{a:+.4f}** · 十分位 ρ = {b:+.4f}")
vals=[v["rho4"] for v in R.values() if v["rho4"] is not None]
fams=sorted({v["fam"] for v in R.values()})
print(f"\n=== 主量 ===\n  纳入 {len(vals)} 组 · 覆盖 **{len(fams)}** 个内容族 {fams}")
print(f"  ρ 中位 **{np.median(vals):+.4f}** · 全距 **{max(vals)-min(vals):.4f}**(判据 < 0.5)")
pos=sum(1 for v in vals if v>0)
print(f"  同号:正 {pos} / {len(vals)}")
def null(seed):
    rng=np.random.default_rng(seed); out=[]
    lab0=pd.cut(d.age,[17,34,49,64,99],labels=False)
    for nm,(its,_f) in KEEP.items():
        perm=rng.permutation(sorted(pd.Series(lab0).dropna().unique()))
        mp={o:n for o,n in zip(sorted(pd.Series(lab0).dropna().unique()),perm)}
        xs,ys=[],[]
        for lv in sorted(pd.Series(lab0).dropna().unique()):
            w,_=weakest(d[lab0==lv],its)
            if np.isfinite(w): xs.append(mp[lv]); ys.append(w)
        if len(xs)>=3: out.append(sp(xs,ys))
    return max(sum(1 for v in out if v>0),sum(1 for v in out if v<0))/max(len(out),1)
nl=float(np.median([null(s) for s in SEEDS]))
obs=max(pos,len(vals)-pos)/len(vals)
print(f"\n=== ⑤ 零:打乱年龄层序后的同号率 = **{nl:.3f}** · 观测 = **{obs:.3f}**")
G=Gate("命名不是测量:作答一致性是真的还是我给巧合起了名字")
p1=G.positive_control("性四题与性别角色三题必须复现 #671 的 +0.80 / +1.00",
    planted=float(0.15-max(abs((R.get('性道德') or {}).get('rho4',0)-0.8),abs((R.get('性别角色') or {}).get('rho4',0)-1.0))),
    floor=0.0,spread=0.01)
p2=G.negative_control("打乱年龄层序后的同号率应低于观测",null=nl,effect=obs,null_spread=0.05,
                      null_kind="随机重排年龄层的顺序 —— 各组之间不该系统同号")
rng_=max(vals)-min(vals)
if p1 and p2:
    verdict=(f"**W1 —— 特质被正面支持:{len(fams)} 个内容族全部同号,ρ 全距 {rng_:.4f} < 0.5**"
             if (len(fams)>=3 and pos in (0,len(vals)) and rng_<0.5) else
             f"**W2 —— 不是统一特质:ρ 全距 {rng_:.4f}(判据 <0.5),`#671` 的命名被推翻**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · 零 {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(kept={k:v for k,v in R.items()},families=fams,median=float(np.median(vals)),
               rng=float(rng_),same_sign_obs=obs,same_sign_null=nl,verdict=verdict,unchallenged=True),
          open(OUT/"is_it_a_trait.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'is_it_a_trait.json'}")

# ── ⚠ 跨仪器闸在本轮开火,而它抓对了:第二具仪器是有的,我没用 ──────────────
# `#658` 的闸判本轮 single_instrument。**豁免语在这里会是假的** ——
# **MFQ 有五个道德域,每个都是一个成块的题组,而且带年龄**(`#653` 已量过它们的最弱一环)。
# ⇒ **补跑第二具仪器。闸的价值正在于它逼我做这件事,而不是让我写一句「换不了仪器」。**
print("\n=== 跨仪器补跑:MFQ 五个道德域(第二具仪器)===")
import pyreadstat as _prs
ITEM={"emotionally":"HARM","weak":"HARM","cruel":"HARM","compassion":"HARM","animal":"HARM","kill":"HARM",
 "treated":"FAIRNESS","unfairly":"FAIRNESS","rights":"FAIRNESS","fairly":"FAIRNESS","justice":"FAIRNESS","rich":"FAIRNESS",
 "lovecountry":"INGROUP","betray":"INGROUP","loyalty":"INGROUP","history":"INGROUP","family":"INGROUP","team":"INGROUP",
 "respect":"AUTHORITY","traditions":"AUTHORITY","chaos":"AUTHORITY","kidrespect":"AUTHORITY","sexroles":"AUTHORITY","soldier":"AUTHORITY",
 "decency":"PURITY","disgusting":"PURITY","god":"PURITY","harmlessdg":"PURITY","unnatural":"PURITY","chastity":"PURITY"}
mfq,_=_prs.read_sav("data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav")
age_col=[c for c in mfq.columns if c.lower()=="age"]
if not age_col:
    print("  ⚠ MFQ 无年龄列 ⇒ 跨仪器补跑做不到,如实记")
else:
    Xm=mfq[list(ITEM)+age_col].dropna()
    Xm=Xm[(Xm[age_col[0]]>=15)&(Xm[age_col[0]]<=90)]
    print(f"  MFQ 30 题 + 年龄 完整 n = **{len(Xm)}** · 年龄范围 [{Xm[age_col[0]].min():.0f},{Xm[age_col[0]].max():.0f}]")
    def wk(fr,items):
        v=[]
        for a,b in combinations(items,2):
            mm=fr[[a,b]].dropna()
            if len(mm)<100 or mm[a].nunique()<2 or mm[b].nunique()<2: continue
            r=sp(mm[a],mm[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(mm[a],mm[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
        return min(v) if v else np.nan
    lab=pd.qcut(Xm[age_col[0]],4,labels=False,duplicates="drop")
    MR={}
    for dom in sorted(set(ITEM.values())):
        its=[k for k,v in ITEM.items() if v==dom]
        xs,ys=[],[]
        for i,lv in enumerate(sorted(pd.Series(lab).dropna().unique())):
            w=wk(Xm[lab==lv],its)
            if np.isfinite(w): xs.append(i); ys.append(w)
        r=sp(xs,ys) if len(xs)>=3 else np.nan
        MR[dom]=float(r) if np.isfinite(r) else None
        print(f"    {dom:10s} 总体最弱 {wk(Xm,its):+.4f} · ρ(年龄四分位, 最弱一环) = **{r:+.4f}**")
    mv=[v for v in MR.values() if v is not None]
    print(f"  MFQ 五域:ρ 中位 **{np.median(mv):+.4f}** · 全距 **{max(mv)-min(mv):.4f}** · 同号 {sum(1 for v in mv if v>0)}/{len(mv)}")
    print(f"  **两具仪器合起来:GSS 全距 1.8000 · MFQ 全距 {max(mv)-min(mv):.4f}** -> "
          f"{'**两具都散开 ⇒ 撤回成立**' if max(mv)-min(mv)>=0.5 else '⚠ MFQ 不散 ⇒ 两具不一致,须单列'}")
    d2=json.load(open(OUT/"is_it_a_trait.json")); d2["mfq"]=MR
    d2["mfq_range"]=float(max(mv)-min(mv)); json.dump(d2,open(OUT/"is_it_a_trait.json","w"),indent=1,ensure_ascii=False)
