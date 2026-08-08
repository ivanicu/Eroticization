"""E03·A23·R118 —— 「更成体系」到底是水平位移,还是真的更一致

**类型:FRONTIER。A22 关弧**(社会侧结构性做不到 + 人侧交叉,决定已安全)**,A23 开弧。**
⚠ A22 的目录名是 `for_which_societies`,而工作在 `#675` 转到了人侧 —— **弧名漂移,不改名**
(738 个账本地址指着它,L81 annotate never rewrite),在这里记下。

`#675` 的 NEXT:Δ 是「更成体系」,但四题是**一起变宽容 / 一起变严厉 / 还是只是彼此更一致**?

## 硬规则①(已跑,在算任何相关之前)
低教育端 `educ<=10` n=1,268 · 高教育端 `educ>=16` n=3,274(总 11,576)。
**均值位移巨大**:`homosex` 1.633→2.698(**+1.064**)· `premarsx` 2.533→3.088(+0.555)·
`teensex` +0.206 · `xmarsex` +0.184。
⚠ **最强混淆量出来了:性四题最弱一对的天花板 0.4658 → 0.5362(+0.0705)**,
而 `#675` 观测的 Δ 是 +0.1349 ⇒ **天花板的移动占观测的 52%。**
(对照:性别角色三题的天花板反而**下降** −0.1295。)

## G1 ESTIMAND —— 三条曲线分开报,不许只报一条
**① 水平** `L(e)` = 四题标准化后的均值(共同水平)· **② 生**`Rw(e)` = **未除天花板**的最弱一环 ·
**③ 归一** `W(e)` = 除过天花板的最弱一环(`#675` 用的那个)· 外加 **④ 天花板本身** `C(e)`。
**主量仍是 Δ,但判据落在「③ 与 ④ 的关系」上。**
## G2 CONTROLS
**正对照**:全样本归一最弱一环复现 `#675` 的 **+0.4732**。
**边际匹配规格(决定性的那一刀)**:在每个窗口内把每题换成**窗内经验 CDF 中位秩** ——
强行把各窗口的边际拉齐 ⇒ **天花板按构造几乎不动**;若 Δ 在这个规格下仍在,
**这条梯度就不可能是边际偏斜。**
**安慰剂**:性别角色三题走同一条流水线(它的天花板是**反向**移动的,所以它同时是一个反向对照)。
## KILL(条件式)
if 正对照复现 and 天花板曲线 C 的 Δ 已单独报告:
  边际匹配规格下 Δ 仍超零 -> **真的更一致** ·
  边际匹配下 Δ 掉到地板 -> **是边际偏斜/水平位移,`#675` 那句话要改写** ·
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
四档序数 + 窗内并列 ⇒ **边际只能拉近,不可能完全拉齐**;因果:横断面无干预。
**跨仪器:MFQ 上已在 `#675` 测过不复制**,本轮是同一具仪器上的分解,**不新增跨仪器主张**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEX=["premarsx","xmarsex","homosex","teensex"]; FEM=["fefam","fepol","fepresch"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
      usecols=["year","educ"]+SEX+FEM, apply_value_formats=False, encoding="latin1")
j=df.dropna(subset=SEX+FEM+["educ"])
educ=j["educ"].to_numpy(float); grid=np.arange(8,19.01,1.0); BW=2.5
def wc(x,y,w):
    mx=np.average(x,weights=w);my=np.average(y,weights=w);cx=x-mx;cy=y-my
    s=np.sqrt(np.average(cx*cx,weights=w)*np.average(cy*cy,weights=w))
    return np.average(cx*cy,weights=w)/s if s>1e-12 else np.nan
def ceil_w(x,y,sign):
    a=np.sort(x);b=np.sort(y); b=b if sign>0 else b[::-1]
    return float(np.corrcoef(a,b)[0,1])
def links(X,W,norm):
    """X: n x k 已排秩。返回最弱一环(norm=True 除天花板)与该组的最小天花板。"""
    v=[];c=[]
    for i,k in combinations(range(X.shape[1]),2):
        r=wc(X[:,i],X[:,k],W)
        if not np.isfinite(r) or abs(r)<1e-12: continue
        idx=W>np.quantile(W,0.5)           # 天花板用窗内主体样本估,避免尾部权重噪声
        cc=ceil_w(X[idx,i],X[idx,k],1 if r>0 else -1)
        if not np.isfinite(cc) or abs(cc)<1e-9: continue
        v.append(r/abs(cc) if norm else r); c.append(abs(cc))
    return (min(v) if v else np.nan),(min(c) if c else np.nan)
def curves(items,marg_match=False):
    R=np.column_stack([pd.Series(j[c]).rank().to_numpy(float) for c in items])
    raw=[];nor=[];cei=[];lev=[]
    Z=np.column_stack([(j[c]-j[c].mean())/j[c].std() for c in items])
    for g in grid:
        W=np.exp(-0.5*((educ-g)/BW)**2)
        if W.sum()<200: raw.append(np.nan);nor.append(np.nan);cei.append(np.nan);lev.append(np.nan); continue
        X=R
        if marg_match:   # 窗内经验 CDF 中位秩 -> 强行拉齐各窗口边际
            idx=W>np.quantile(W,0.5); X=np.column_stack([
                pd.Series(np.where(idx,j[c],np.nan)).rank(pct=True).fillna(0.5).to_numpy() for c in items])
        r,_=links(X,W,False); n_,c_=links(X,W,True)
        raw.append(r);nor.append(n_);cei.append(c_);lev.append(float(np.average(Z.mean(1),weights=W)))
    return map(np.array,(raw,nor,cei,lev))
def D(y):
    m=np.isfinite(y)
    return float(np.mean(y[m][-3:])-np.mean(y[m][:3])) if m.sum()>=4 else np.nan
res={}
for nm,items in [("性四题",SEX),("性别角色三题",FEM)]:
    raw,nor,cei,lev=curves(items)
    print(f"\n=== {nm} · n = {len(j):,} ===")
    print(f"  ① 水平 L   {[f'{v:+.3f}' for v in lev]}   Δ = **{D(lev):+.4f}**")
    print(f"  ② 生 Rw    {[f'{v:+.3f}' for v in raw]}   Δ = **{D(raw):+.4f}**")
    print(f"  ④ 天花板 C {[f'{v:+.3f}' for v in cei]}   Δ = **{D(cei):+.4f}**")
    print(f"  ③ 归一 W   {[f'{v:+.3f}' for v in nor]}   Δ = **{D(nor):+.4f}**  <- #675 用的那个")
    mraw,mnor,mcei,_=curves(items,marg_match=True)
    print(f"  ★ 边际匹配后的归一 {[f'{v:+.3f}' if np.isfinite(v) else 'na' for v in mnor]}   Δ = **{D(mnor):+.4f}**")
    rng=np.random.default_rng(20260806); nul=[]
    for _ in range(300):
        e0=educ.copy(); np.random.default_rng().shuffle(e0)
        globals()['educ']=rng.permutation(e0)
        _,mn2,_,_=curves(items,marg_match=True); nul.append(abs(D(mn2)))
    globals()['educ']=j["educ"].to_numpy(float)
    nul=np.array(nul); q=float(np.nanquantile(nul,0.95)); p=float(np.nanmean(nul>=abs(D(mnor))))
    print(f"     零 95% 分位 {q:.4f} · p = **{p:.4f}**  {'✅ 边际匹配后仍在' if p<0.05 else '⛔ 边际匹配后掉到地板'}")
    res[nm]=dict(level=D(lev),raw=D(raw),ceiling=D(cei),norm=D(nor),marg=D(mnor),q95=q,p=p)
G=Gate("「更成体系」是水平位移还是真的更一致")
p1=G.positive_control("全样本归一最弱一环复现 #675 的 +0.4732(>0.20)",
                      planted=float(np.nanmean(list(curves(SEX))[1])),floor=0.20,spread=0.01)
s=res["性四题"]
if p1:
    verdict=(f"**真的更一致:边际匹配后 Δ = {s['marg']:+.4f},p = {s['p']:.4f}**" if s['p']<0.05
             else f"**改写 `#675`:边际匹配后 Δ = {s['marg']:+.4f} 掉到地板(p = {s['p']:.4f})—— 那条梯度里有边际偏斜**")
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(res=res,verdict=verdict,unchallenged=True),open(OUT/"level_vs_coherence.json","w"),indent=1,ensure_ascii=False)
