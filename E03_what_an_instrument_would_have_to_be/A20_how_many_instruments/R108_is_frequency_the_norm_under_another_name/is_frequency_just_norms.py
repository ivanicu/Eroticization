"""E03·A20·R108 —— 那把「稀有度」尺,是不是「规范」换了个说法

**类型:FRONTIER**。`#665` 的 NEXT 要求穷举第二份频率编码,而穷举顺带打开了一个更该问的问题。

## ① 穷举(已跑,写在这里)
SCCS 全库 45 个含频率词的变量里,**对象是性实践的 7 个,全部出自 `broude1976cross`**:
`SCCS160` 婚内 · `SCCS166/167` 婚前 男/女 · `SCCS170/171` 婚外 男/女 · `SCCS174` 强奸 · `SCCS177` 同性恋。
**⇒ 没有第二个项目编过性实践频率。** 正对照:穷举重新找到 `SCCS166`/`SCCS167` ✅

**⇒ W2:「整条关系的跨项目复制」做不到。**
**可证伪形式(写成列名条件)**:*若 SCCS 或任一跨文化库出现一个 `source ≠ broude1976cross`
且标题含 frequency/prevalence 且对象为某性实践的变量,与 `SCCS282` 共同 n ≥ 60,这一条即被推翻。*

## ⑤ 读码(已跑):定义污染被排除
七个变量的码文**一个规范词都没有**(`Universal / Moderate / Occasional / Uncommon`)——
不像 `#665` 剔掉的 Frayser `SCCS961` 码 7(*strongly disapproved **and rare***)。

## ⚠ 但码文干净 ≠ 编码者没从规范里推。这一层可以直接测,而它是本轮的主体

⚠ **BASIN**:`#665` 刚记下「我的先验比数据悲观」⇒ 本轮**不反射性赌悲观**;
   **而真正不受欢迎的是 W3 —— 它会把「稀有度→谴责」从一条关系降成一条定义。下注 W3。**

**W3 若为真**:Broude 的频率只是他自己的态度码换个说法 ⇒
  **控制住 Broude 的态度之后,频率对 Murdock 那一侧的预测力应当消失。**
**W-indep 若为真**:频率携带**规范之外**的信息,而 Murdock 的编码者也读到了同样的东西 ⇒ 偏相关仍为正。
**Meta-separator**:W3 成立会让 `#665` 的 +0.6725 有一部分是同义反复,**而那正是本轮要冒的险。**

## G1 ESTIMAND(先于方法)
**主量 = 偏相关 `ρ(SCCS167 频率, SCCS282 Murdock 规范 | SCCS165 Broude 态度)`**(Spearman 秩偏相关)。
## G2 CONTROLS
**正对照**:`ρ(SCCS165, SCCS282 | SCCS167)` 必须仍为正 —— **规范应当在频率之外仍预测规范**;
  它若也塌,说明偏相关这个工具在这里没分辨力。
**安慰剂**:把控制变量换成与性无关的社会变量,主量应基本不动。
  **这个零该不该是零?** 不该 —— 换一个无关的控制变量,偏相关本就该回到原相关,
  **它是一条系统性基线,不是应当趋零的干扰** ⇒ **用 `offset_control` 并命名零的种类。**
## G3/G4:三对(频率×Murdock、态度×Murdock、频率×态度)的原始与偏相关全报;女/男两条规格。
## KILL(条件式)
if 正对照仍为正 and 安慰剂偏相关 ≈ 原相关:
  主量的块 bootstrap 区间**含零** -> **W3:频率是规范换了个说法,`#665` 必须缩小**
  不含零且为正 -> **W-indep:频率携带规范之外的信息**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**偏相关只能排除「线性/单调可解释的那部分」**,不能证明编码者主观上没参考规范 ·
四个项目读同一批民族志 · **跨仪器:换不了仪器,只此一具**(`#664` 穷举 D-PLACE 五库已验证)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]
B="data/external/dplace/repo/datasets/SCCS/"
D=pd.read_csv(B+"data.csv"); SOC=pd.read_csv(B+"societies.csv")
W=D.pivot_table(index="soc_id",columns="var_id",values="code",aggfunc="first")
BLK={r.id:(int(np.floor(r.Lat/10)),int(np.floor(r.Long/10))) for r in SOC.dropna(subset=["Lat","Long"]).itertuples()}
MUR=(-W["SCCS282"].where(~W["SCCS282"].isin([1.0,5.0]))).rename("mur")   # 定向:越大越谴责

def prep(freq,att):
    m=pd.concat([W[freq].rename("f"),W[att].rename("a"),MUR],axis=1).dropna()
    return m
def pcorr(m,x,y,z):
    rx=m[x].rank(); ry=m[y].rank(); rz=m[z].rank()
    ex=rx-np.poly1d(np.polyfit(rz,rx,1))(rz); ey=ry-np.poly1d(np.polyfit(rz,ry,1))(rz)
    return float(spearmanr(ex,ey).statistic)
def boot(m,fn,B_=600):
    bl=sorted({BLK.get(x,("na","na")) for x in m.index}); by={b:[x for x in m.index if BLK.get(x,("na","na"))==b] for b in bl}
    out=[]
    for sd in SEEDS:
        rng=np.random.default_rng(sd)
        for _ in range(B_//len(SEEDS)):
            socs=[x for i in rng.integers(0,len(bl),len(bl)) for x in by[bl[i]]]
            try: v=fn(m.loc[socs])
            except Exception: continue
            if np.isfinite(v): out.append(v)
    return np.quantile(out,[.025,.975]) if out else (np.nan,np.nan)

print("=== G3/G4:女(SCCS167/165)与男(SCCS166/SCCS164?)两条规格 ===")
SPECS=[("女","SCCS167","SCCS165")]
if "SCCS164" in W.columns: SPECS.append(("男","SCCS166","SCCS164"))
res={}
for tag,freq,att in SPECS:
    m=prep(freq,att)
    r_fm=float(spearmanr(m.f,m["mur"]).statistic); r_am=float(spearmanr(m.a,m["mur"]).statistic)
    r_fa=float(spearmanr(m.f,m.a).statistic)
    p_fm=pcorr(m,"f","mur","a"); p_am=pcorr(m,"a","mur","f")
    lo,hi=boot(m,lambda d:pcorr(d,"f","mur","a"))
    lo2,hi2=boot(m,lambda d:pcorr(d,"a","mur","f"))
    res[tag]=dict(n=int(len(m)),r_freq_mur=r_fm,r_att_mur=r_am,r_freq_att=r_fa,
                  partial_freq=p_fm,ci_freq=[float(lo),float(hi)],
                  partial_att=p_am,ci_att=[float(lo2),float(hi2)])
    print(f"\n  【{tag}】n={len(m)}")
    print(f"    原始 频率×Murdock **{r_fm:+.4f}** · 态度×Murdock {r_am:+.4f} · 频率×态度 {r_fa:+.4f}")
    print(f"    **偏相关 频率×Murdock | 态度 = {p_fm:+.4f}**  95%CI [{lo:+.4f}, {hi:+.4f}]")
    print(f"    正对照 态度×Murdock | 频率 = {p_am:+.4f}  95%CI [{lo2:+.4f}, {hi2:+.4f}]")

t="女"; R=res[t]; m=prep("SCCS167","SCCS165")
alt=[c for c in ["SCCS31","SCCS63","SCCS1","SCCS64"] if c in W.columns]
off=np.nan
for c in alt:
    mm=pd.concat([W["SCCS167"].rename("f"),W[c].rename("a"),MUR],axis=1).dropna()
    if len(mm)>=50:
        off=pcorr(mm,"f","mur","a")
        print(f"\n=== offset:把控制变量换成与性无关的 {c}(n={len(mm)})===")
        print(f"  原始 {spearmanr(mm.f,mm['mur']).statistic:+.4f} -> 偏相关 **{off:+.4f}**(应基本不动)")
        break
G=Gate("那把稀有度尺,是不是规范换了个说法")
p1=G.positive_control("态度×Murdock | 频率 必须仍为正",planted=R["partial_att"],floor=0.0,spread=0.03)
p2=G.offset_control("换一个与性无关的控制变量,偏相关应回到原相关",
                    effect=R["r_freq_mur"], offset=float(off), spread=0.05,
                    null_kind="与性道德无关的社会结构变量作控制 —— 系统性基线,不该趋零")
if p1 and p2:
    lo,hi=R["ci_freq"]
    verdict=("**W3 —— 频率是规范换了个说法:控制住 Broude 自己的态度后,偏相关 "
             f"{R['partial_freq']:+.4f},区间含零。`#665` 必须缩小。**" if lo<0<hi else
             f"**W-indep —— 频率携带规范之外的信息:偏相关 {R['partial_freq']:+.4f},区间不含零**")
else: verdict=f"UNVERIFIED —— 控制未齐(正 {p1} · offset {p2})"
print(f"\n{verdict}"); print(G)
json.dump(dict(specs=res,offset=float(off) if np.isfinite(off) else None,verdict=verdict,
               enumeration="7 个性实践频率变量全部 broude1976cross;码文 0 个规范词",
               falsifier="若出现 source≠broude1976cross 的性实践频率变量且与 SCCS282 共同 n≥60,W2 即被推翻",
               unchallenged=True),
          open(OUT/"is_frequency_just_norms.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'is_frequency_just_norms.json'}")

# ── ⚠ 标注的更正:我把控制工具的方向用反了 ────────────────────────────────
# 预注册写的是「换一个与性无关的控制变量,主量应**基本不动**」—— 那是一个**不变性检查**,
# 而 `offset_control` 检验的是「效应必须**超出**系统性基线」。**方向反了,是选工具的错,不是判据的错。**
# (与 `#639` 的单位错配、`#654` 的转录漏 `|·|` 同族:预注册的**文字**是对的,**代码**实现错了。)
# ⇒ 正确的量:**无关控制造成的变化**必须远小于**真控制造成的变化**。两者在上面同一次运行里已经算出。
print("\n=== 标注的更正:改用正确的不变性检查 ===")
chg_real = abs(R["r_freq_mur"]-R["partial_freq"])
chg_null = abs(R["r_freq_mur"]-off)
print(f"  真控制(Broude 态度)造成的变化 = **{chg_real:.4f}**")
print(f"  无关控制({c})造成的变化      = **{chg_null:.4f}**  ({chg_null/chg_real*100:.1f}% of 真控制)")
G2=Gate("那把稀有度尺,是不是规范换了个说法(更正后的控制)")
q1=G2.positive_control("态度×Murdock | 频率 必须仍为正",planted=R["partial_att"],floor=0.0,spread=0.03)
q2=G2.negative_control("不变性:换一个与性无关的控制变量,主量的变化应远小于真控制造成的变化",
                       null=chg_null, effect=chg_real, null_spread=0.02,
                       null_kind="与性道德无关的社会结构变量作控制 —— 它不该改变任何东西")
lo,hi=R["ci_freq"]
if q1 and q2:
    v2=("**W3 —— 频率是规范换了个说法,`#665` 必须缩小**" if lo<0<hi else
        f"**W-indep —— 频率携带规范之外的信息:偏相关 {R['partial_freq']:+.4f},95% CI [{lo:+.4f},{hi:+.4f}] 不含零**")
else: v2="UNVERIFIED"
print(f"\n{v2}"); print(G2)
d=json.load(open(OUT/"is_frequency_just_norms.json"))
d.update(dict(verdict_corrected=v2, change_real=float(chg_real), change_null=float(chg_null),
              correction="预注册文字是不变性检查,而我用了 offset_control(方向相反)。选工具的错,不是判据的错。"))
json.dump(d,open(OUT/"is_frequency_just_norms.json","w"),indent=1,ensure_ascii=False)
