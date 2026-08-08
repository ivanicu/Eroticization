"""E03·A22·R116 —— 「一个社会对所有孩子一视同仁」,对哪些社会成立

**类型:FRONTIER**。**A21 关弧,A22 开弧。**
`#673` 的 NEXT:回到对象,**只许问一条声明**。而 ① 的清点结果是:
**从未分层的站得住声明有 8 条,而「社会」这个单位一条都没问过。**
⇒ 选页面头号的社会侧声明:**换对象 +0.845**(体罚四对象两两中位,`#640`)。

## 硬规则①(已跑)
`SCCS63`「社区规模」**8 档有序**,与体罚四件套联合 **n = 139**(3 档 47/46/46 · 4 档 35/35/34/35)。
⚠ **`SCCS61`「定居性」的码不单调**(码 5「不定期迁移」夹在 4 与 6 之间)⇒ **不可当有序分层,弃用。**
安慰剂用的疼爱四件套与 `SCCS63` 联合 **n = 152**。

⚠ **BASIN**:W2(随规模变)是有故事的那个 ⇒ **下注 W1/W3。**
W1 普遍(梯度含零)—— **这也是结果,不是失败** · W2 随社会规模变 ·
**W3 = 最强混淆:梯度是民族志详略的伪影** ——
  **社会越大、民族志越详、编码者越有把握 ⇒ 四个对象之间自然更一致。这是编码伪影,不是心理。**

## G1 ESTIMAND
每个社区规模层内:**体罚四对象的最弱一环**(天花板归一后最小的那一对)。
**主量 = `ρ(规模层序, 最弱一环)`。**
## G2 CONTROLS
**正对照**:全样本最弱一环必须复现 `#640` 的 **+0.7924**(六对中的最小者,容差 0.03)。
**安慰剂 = W3 的控制**:**同一分层作用在「疼爱」四对象上的同一个 ρ**。
  **这个零该不该是零?** 该 —— 若梯度是关于**体罚**的,疼爱上不该有同样的梯度;
  **若两者相同 ⇒ 是民族志详略,不是体罚** ⇒ `negative_control`。
**零(补)**:**打乱层序**后的 ρ 分布,报观测的经验 p。
## G3/G4:3 档与 4 档两条规格 · 主与安慰剂各报全部层值。
## KILL(条件式)
if 正对照复现 and |安慰剂 ρ| < 0.5×|主 ρ| and 观测 ρ 超出打乱零的 95% 分位:
  -> **W2 随规模变**;安慰剂不小于一半 -> **W3 民族志详略**;
  观测未超零分位 -> **W1 对所有社会一样(这也是结果)**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**每层仅 46 个社会**(3 档)⇒ 功率有限,**只可能测到大梯度** ·
**「民族志详略」没有被直接测量的量** ⇒ 只能用疼爱四件套作代理 ·
**跨仪器:换不了仪器,只此一具**(`#664` 穷举五库;`#667` 已证 Lang 只有一个对象)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; FLOOR=30
B="data/external/dplace/repo/datasets/SCCS/"
D=pd.read_csv(B+"data.csv"); W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
CORP=["SCCS453","SCCS454","SCCS455","SCCS456"]; AFF=["SCCS469","SCCS470","SCCS471","SCCS472"]; SZ="SCCS63"
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if s<0: y=y[::-1]
    return sp(x,y)
def weakest(fr,items):
    v=[]
    for a,b in combinations(items,2):
        m=fr[[a,b]].dropna()
        if len(m)<FLOOR or m[a].nunique()<2 or m[b].nunique()<2: continue
        r=sp(m[a],m[b])
        if not np.isfinite(r) or r==0: continue
        c=rmax(m[a],m[b],1 if r>0 else -1)
        if np.isfinite(c) and abs(c)>1e-9: v.append(r/abs(c))
    return (min(v) if v else np.nan), len(v)
ovc,_=weakest(W,CORP); ova,_=weakest(W,AFF)
print(f"=== 正对照 ===\n  体罚四件套 全样本最弱一环 = **{ovc:+.4f}**(`#640` 六对最小者 +0.7924,归一后应略高)")
print(f"  疼爱四件套 全样本最弱一环 = **{ova:+.4f}**(安慰剂须成块)")
def grad(items,k):
    j=W[[SZ]+items].dropna()
    lab=pd.qcut(j[SZ].rank(method='first'),k,labels=False,duplicates='drop')
    xs,ys,ns=[],[],[]
    for i,lv in enumerate(sorted(pd.Series(lab).dropna().unique())):
        sub=j[lab==lv]; w,npair=weakest(sub,items)
        if np.isfinite(w): xs.append(i); ys.append(w); ns.append(len(sub))
    return (sp(xs,ys) if len(xs)>=3 else np.nan), list(zip(xs,ys,ns))
print("\n=== G3/G4:3 档与 4 档 ===")
res={}
for k in (3,4):
    a,pa=grad(CORP,k); b,pb=grad(AFF,k)
    res[k]=dict(main=float(a) if np.isfinite(a) else None,placebo=float(b) if np.isfinite(b) else None,
                main_pts=[[int(x),float(y),int(n)] for x,y,n in pa],
                placebo_pts=[[int(x),float(y),int(n)] for x,y,n in pb])
    print(f"  {k} 档 · 体罚各层最弱一环 {[f'{y:+.3f}(n={n})' for _,y,n in pa]}")
    print(f"        主 ρ = **{a:+.4f}** · 安慰剂(疼爱) ρ = **{b:+.4f}** {[f'{y:+.3f}' for _,y,_ in pb]}")
def null(seed,k=3,B_=4000):
    rng=np.random.default_rng(seed); j=W[[SZ]+CORP].dropna()
    lab=pd.qcut(j[SZ].rank(method='first'),k,labels=False,duplicates='drop')
    ys=[]
    for lv in sorted(pd.Series(lab).dropna().unique()):
        w,_=weakest(j[lab==lv],CORP)
        ys.append(w)
    ys=[y for y in ys if np.isfinite(y)]
    out=[]
    for _ in range(B_):
        out.append(abs(sp(list(range(len(ys))),list(rng.permutation(ys)))))
    return np.array(out)
nd=np.concatenate([null(s) for s in SEEDS]); q95=float(np.quantile(nd,0.95))
obs=abs(res[3]["main"]); p=float((nd>=obs).mean())
print(f"\n=== 零:打乱层序 ===\n  零 95% 分位 **{q95:.4f}** · 观测 |ρ| **{obs:.4f}** · 经验 p = **{p:.4f}**")
G=Gate("「一个社会对所有孩子一视同仁」对哪些社会成立")
p1=G.positive_control("体罚四件套全样本最弱一环必须复现 #640(>0.75)",planted=float(ovc),floor=0.75,spread=0.01)
p2=G.negative_control("安慰剂:疼爱四件套上不该有同样的梯度(若有 ⇒ 是民族志详略)",
                      null=abs(res[3]["placebo"]),effect=obs,null_spread=0.05,
                      null_kind="同一分层作用在另一组四对象上 —— 若梯度相同,那是民族志详略而非体罚")
if p1 and p2:
    verdict=(f"**W2 —— 随社会规模变:ρ = {res[3]['main']:+.4f},超出打乱零的 95% 分位 {q95:.4f}(p={p:.4f})**"
             if obs>q95 else
             f"**W1 —— 对所有社会一样:|ρ| {obs:.4f} 未超零的 95% 分位 {q95:.4f}(p={p:.4f})。这也是结果。**")
elif p1: verdict=f"**W3 —— 民族志详略:安慰剂(疼爱) ρ = {res[3]['placebo']:+.4f},不小于主效应的一半**"
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(overall_corp=float(ovc),overall_aff=float(ova),specs=res,null_q95=q95,obs=obs,p=p,
               verdict=verdict,unchallenged=True),open(OUT/"which_societies.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'which_societies.json'}")
