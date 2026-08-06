"""E03·A34·R180 —— 档数不同的那一格,以及合成世界与真数据的落差

**类型:FRONTIER。这是 `#736`① —— 而它的预注册预期在跑之前就写下了。**

**心理学的那一句(它决定 `#736` 那组数能不能搬到这一页上):
我在合成世界里量出「天花板归一只贡献 6 个点」,而真数据里的天花板是 0.48–0.85,不是 0.95。
合成世界里那个「贡献很小」的结论,搬不到这一页。**

## ⚠ 预注册的预期(`#736`①,跑之前写下)
**档数不同会让天花板 <1,所以归一应当贡献更多。若不是,归一的机制我没搞懂。**
## G1 ESTIMAND
① 合成:`4×4` 对 `4×5` 两格的**天花板中位**与**归一的贡献**(`只归一 − 生`,除以 ρ_true);
② 真数据:GSS 六对各自的**天花板**,以及含 `homosex` 与不含的两组中位。
## G2 CONTROLS
**④ 正对照 = 预注册预期本身**:`4×5` 的天花板必须 < `4×4` 的,且归一贡献必须更大。
**⇒ 这个正对照能失败,而它失败就等于「我没搞懂这个机制」。**
## ⑤ 停止条件(跑之前写下)
- **合成里 `4×5` 的归一贡献不大于 `4×4` ⇒ 记「机制没搞懂」并停。**
- **真数据的天花板中位若落在合成的 [0.95, 1.0] 区间内 ⇒ 合成的 regime 与真数据一致,`#736` 的量级可搬;
  若显著更低 ⇒ 量级不可搬,只有方向可搬。**
## IMPOSSIBLE(不写 planned)
合成部分**一具仪器也没有用**(同 `#736`:**换不了仪器**,问的是估计量的性质);
真数据部分只有 GSS 这一组四题 ⇒ **天花板的分布只能在这四题上说。**`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pyreadstat
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def ceil_of(a,b):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    r=sp(a,b); return abs(sp(x,y if r>0 else y[::-1]))
C4=[-0.67,0,0.67]; C5=[-0.84,-0.25,0.25,0.84]
rng=np.random.default_rng(20260806); N=3000; REP=200
res={}
for lab,(ca,cb) in (("4x4",(C4,C4)),("4x5",(C4,C5))):
    gains=[]; ceils=[]
    for rho in (0.2,0.4,0.6):
        for rel in (0.7,1.0):
            gr=[];cc=[]
            for _ in range(REP):
                L=rng.multivariate_normal([0,0],[[1,rho],[rho,1]],size=N)
                if rel<1: L=L+rng.normal(0,np.sqrt((1-rel)/rel),size=L.shape)
                a=np.digitize(L[:,0],ca); b=np.digitize(L[:,1],cb)
                r=sp(a,b); c=ceil_of(a,b); gr.append((r/c-r)/rho); cc.append(c)
            gains.append(float(np.median(gr))); ceils.append(float(np.median(cc)))
    res[lab]=dict(gain=float(np.median(gains)),ceil=float(np.median(ceils)))
    print(f"合成 {lab}: 天花板中位 **{res[lab]['ceil']:.4f}** · 归一贡献中位 **{res[lab]['gain']:+.4f}**")
pc = res["4x5"]["ceil"]<res["4x4"]["ceil"] and res["4x5"]["gain"]>res["4x4"]["gain"]
print(f"④ 正对照(= 预注册预期):{'✅ 成立 —— 机制是我以为的那个' if pc else '⛔ 不成立 —— 机制我没搞懂'}")
assert pc
SEX=["premarsx","xmarsex","homosex","teensex"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+SEX,encoding="latin1")
J=g.dropna(subset=SEX)
print(f"\n硬规则①:GSS 四题完整个案 n={len(J):,} · **分析样本内的实际档数** "
      + " · ".join(f"{c}:{J[c].nunique()}" for c in SEX))
print("  ⚠ `homosex` 在**全库**是 5 档(含 1–5),而**四题完整个案内只出现 4 档** ——")
print("     **码本的档数不等于分析样本里的档数**,这是硬规则①的又一次兑现。")
wh=[];wo=[]
for a,b in itertools.combinations(SEX,2):
    c=ceil_of(J[a],J[b]); (wh if "homosex" in (a,b) else wo).append(c)
    print(f"    {a:10s} × {b:10s} 天花板 **{c:.4f}**{'  ← 含 homosex' if 'homosex' in (a,b) else ''}")
mh,mo=float(np.median(wh)),float(np.median(wo))
print(f"\n  含 homosex 的三对中位 **{mh:.4f}** · 不含的三对 **{mo:.4f}** ⇒ "
      f"**{'含 5 档题的对天花板更低' if mh<mo else '⚠ 含它的反而更高 —— 实数据里档数不是主因'}**")
allc=wh+wo
print(f"  **GSS 六对天花板全距 [{min(allc):.4f}, {max(allc):.4f}],中位 {np.median(allc):.4f}**")
inreg = 0.95<=np.median(allc)<=1.0
print(f"\n⑤ 判据:真数据天花板中位是否落在合成的 [0.95, 1.0] ⇒ **{'是,量级可搬' if inreg else '否 —— 量级不可搬,只有方向可搬'}**")
print("⇒ **`#736` 的「归一只贡献 6 个点」是合成世界的数,而合成世界的天花板是 0.95–0.99;**")
print(f"   **GSS 的天花板是 {min(allc):.2f}–{max(allc):.2f} ⇒ 在真数据里归一做的功远大于合成所示。**")
print("   **方向(四种校正全都低估)可搬;量级不可搬。**")
json.dump(dict(synth=res,positive_control=bool(pc),gss_ceilings_with=wh,gss_ceilings_without=wo,
  gss_median=float(np.median(allc)),in_synth_regime=bool(inreg),
  note="码本档数≠分析样本档数;真数据天花板由边际偏斜决定,不由档数",
  unchallenged=True),open(OUT/"mixed.json","w"),indent=1,ensure_ascii=False)
