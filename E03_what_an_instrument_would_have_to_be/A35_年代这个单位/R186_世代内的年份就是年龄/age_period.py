"""E03·A35·R186 —— 世代内部,「年份」就是「年龄」

**类型:FRONTIER。而它攻的是我自己三轮前刚立的那条。**

**心理学的那一句(本轮要判的):`#740` 说融合发生在世代内部,不是老一代退出。
但在一个出生世代内部,**年份与年龄是同一个变量** ——
所以那句话也可能是「人老了之后才把这两件事连起来」。哪一个?**

## ⚠ 元分离器先写下,因为它可能让 `#740` 的措辞整个作废
在固定出生世代里,**`age = year − cohort` 是恒等式** ⇒ **年份与年龄完全共线。**
⇒ **`#740` 的「世代内斜率」按构造是「年代 + 年龄」的合体,不是纯年代。**
**这不是一个可以被更好的统计量修好的问题,是识别问题(APC),必须如实登记。**
⚠ **而 `#742`① 那条元层 NEXT 本轮**推后一次**,理由:BASIN 规则要求设计一个
「正结果我会不高兴」的步,而本轮正是它 —— 它能推翻页上最新那一行。已记,不是漏。**

## 但两刀合起来能排除两个纯解释(这是本轮能做的,也是它的上限)
- `#740`(**固定世代,年龄随年份变**):斜率为正 ⇒ **不是纯队列替换**(已判,世代内占 113%);
- **本轮(固定年龄段,世代随年份变)**:斜率若仍为正 ⇒ **不是纯年龄**;
- **两者都为正 ⇒ 唯一同时出现在两刀里的是「年代」。** ⚠ **这是排除,不是识别。**

## 硬规则①(脚本内先印)
`age` 的 n 与范围;各年龄段逐年 n;以及世代内 `educ` 的漂移(它是次级混淆)。

## G1 ESTIMAND
`corr(premarsx, homosex)` 对年份的斜率,**在固定年龄段内**(25–39 · 40–54 · 55–69),按 n 加权合并。
## W1–W4(双边)
| 世界 | 年龄段内合并斜率 | 读法 |
|---|---|---|
| **W1 年代** | **≥ 世代内斜率的 60%** | 两刀都为正 ⇒ 排除纯年龄与纯队列,剩下年代 |
| **W2 其实是年龄** | **≤ 25%** | **`#740` 的措辞要改成「人老了才连起来」** |
| **W3 之间** | 25–60% | 报份额 |
| **W4 判不了** | 年龄段内斜率落在自己的零里 | 这具仪器分不开 |

⚠ **W2 的正结果直接削页上最新那一行 —— 这正是本轮的目的。**

## G2 CONTROLS
**④ 正对照**:不分层时必须复现 `#740` 的世代内 **+0.00672/年**(容差 0.0005)。
**零** = `negative_control`,**零的种类 = 在每个年龄段内部打乱受访者的年份标签 ——
保住年龄构成、每年 n 与作答分布,只毁掉「同一年龄段里谁属于哪一年」。**
**PLACEBO(必须,且它防的是本轮自己的陷阱)**:分层会减少每格 n ⇒ 斜率可能因噪声缩小。
⇒ **按性别做一次同样粒度的分层**(与年龄段同样的组数与 n 损失),**若它也把斜率压掉同样多,
那压掉的是 n,不是年龄。**
**次级控制**:在(年龄段 × 教育三分位)内再算一次 —— **`#676` 已证教育本身产生同一个融合。**
## G3:3 个年龄段 × {不分教育, 分教育} + 性别安慰剂 全报。G4:年龄段切法两种(3 档 / 5 档)。
## ⑤ 停止条件(**双边**,跑之前写死)
- **不分层复现不到 0.0005 ⇒ UNVERIFIED 并停。**
- **年龄段内斜率落在零里 ⇒ W4;≥60% ⇒ W1;≤25% ⇒ W2;之间 ⇒ W3。**
- **安慰剂(性别分层)把斜率压到与年龄分层同样低 ⇒ 判「压掉的是 n」,本轮的年龄结论作废。**
## IMPOSSIBLE(不写 planned)
**APC 不可识别**:年龄、年代、队列三者线性相关,**本轮只能排除两个纯解释,不能估出年代的系数**;
仍是**重复横断面**;**换不了仪器**。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
A,Bc="premarsx","homosex"
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","cohort","age","educ","sex",A,Bc],encoding="latin1")
J=g.dropna(subset=[A,Bc,"cohort"]).copy()
J["gen"]=pd.cut(J.cohort,[1880,1928,1946,1965,1981,1997,2010],
                labels=["前1929","1929–45","婴儿潮46–64","X 65–80","千禧81–96","Z 97+"])
GENS=["1929–45","婴儿潮46–64","X 65–80","千禧81–96"]; FN,FY=150,5
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def pooled_slope(frame,key,levels):
    """在 key 的每个层内算 corr~year 斜率,按 n 加权合并。"""
    num=[];wt=[]
    for lv in levels:
        fr=frame[frame[key]==lv][["year",A,Bc]].dropna(); pts=[]
        for y,sub in fr.groupby("year"):
            if len(sub)<FN or sub[A].nunique()<2 or sub[Bc].nunique()<2: continue
            pts.append((float(y),sp(sub[A],sub[Bc]),len(sub)))
        if len(pts)>=FY:
            num.append(float(np.polyfit([p[0] for p in pts],[p[1] for p in pts],1)[0])); wt.append(sum(p[2] for p in pts))
    return (float(np.average(num,weights=wt)) if num else np.nan), len(num)
print("=== 硬规则① ===")
print(f"  两题+cohort 非缺失 n={len(J):,} · age 非缺失 {J.age.notna().sum():,} · "
      f"educ 非缺失 {J.educ.notna().sum():,} · sex 非缺失 {J.sex.notna().sum():,}")
Ja=J.dropna(subset=["age"]).copy()
Ja["ageb"]=pd.cut(Ja.age,[24,39,54,69],labels=["25–39","40–54","55–69"])
print(f"  年龄段 n:"+" · ".join(f"{l} {int((Ja.ageb==l).sum()):,}" for l in ["25–39","40–54","55–69"]))
for l in ["25–39","40–54","55–69"]:
    yy=Ja[Ja.ageb==l].groupby("year").size()
    print(f"    {l}:{int((yy>=FN).sum())} 个 n≥{FN} 的调查年(共 {yy.size})")
print("\n=== ④ 正对照:不分层的世代内斜率必须复现 `#740` 的 +0.00672 ===")
s_gen,ng=pooled_slope(J,"gen",GENS)
print(f"  世代内合并 **{s_gen:+.5f}**({ng} 个世代)· 账本 +0.00672 · 差 {abs(s_gen-0.00672):.5f} "
      f"{'✅' if abs(s_gen-0.00672)<=0.0005 else '⛔ ⑤ 触发'}")
if abs(s_gen-0.00672)>0.0005:
    print("⛔ 停"); sys.exit(0)
print("\n=== 本轮的刀:固定年龄段,世代随年份变 ===")
s_age,na=pooled_slope(Ja,"ageb",["25–39","40–54","55–69"])
print(f"  年龄段内合并 **{s_age:+.5f}**({na} 段)⇒ **占世代内斜率的 {s_age/s_gen:.0%}**")
for l in ["25–39","40–54","55–69"]:
    s,_=pooled_slope(Ja,"ageb",[l]); print(f"    {l}: {s:+.5f}")
print("\n=== PLACEBO:按性别分层(同样的组数量级与 n 损失)===")
Js=J.dropna(subset=["sex"]).copy()
s_sex,ns=pooled_slope(Js,"sex",sorted(Js.sex.unique()))
print(f"  性别内合并 **{s_sex:+.5f}**({ns} 组)⇒ 占世代内的 {s_sex/s_gen:.0%}")
print("  ⇒ **若年龄分层把斜率压掉而性别分层没有,压掉的才是年龄;两者一起压掉,压掉的是 n。**")
print("\n=== 次级控制:年龄段 × 教育三分位(`#676` 已证教育本身产生同一个融合)===")
Je=Ja.dropna(subset=["educ"]).copy()
Je["edt"]=pd.qcut(Je.educ,3,labels=["低","中","高"],duplicates="drop")
# ⚠ `pd.cut`/`qcut` 对区间外的行给 NaN;先把两个分层键的 NaN 行去掉,再拼字符串,
#   否则 `astype(str)` 会造出 "nan" 层,而 `unique()` 里混着 float 会让 sorted() 直接崩。
Je=Je.dropna(subset=["ageb","edt"]).copy()
Je["cell"]=Je.ageb.astype(str)+"|"+Je.edt.astype(str)
s_ae,nae=pooled_slope(Je,"cell",sorted(Je.cell.unique()))
print(f"  (年龄段 × 教育)内合并 **{s_ae:+.5f}**({nae} 格)⇒ 占世代内的 {s_ae/s_gen:.0%}")
rng=np.random.default_rng(20260806); nul=[]
for _ in range(200):
    P=Ja.copy(); P["year"]=P.groupby("ageb",observed=True)["year"].transform(lambda s: rng.permutation(s.to_numpy()))
    v,_=pooled_slope(P,"ageb",["25–39","40–54","55–69"])
    if np.isfinite(v): nul.append(v)
q=np.quantile(nul,[0.025,0.975])
print(f"\n=== 零(年龄段内打乱年份,B={len(nul)})===")
print(f"  95% 区间 [{q[0]:+.5f}, {q[1]:+.5f}] · 实测 {s_age:+.5f} ⇒ "
      f"{'✅ 在零之外' if not (q[0]<=s_age<=q[1]) else '⚠ 落在零里 ⇒ W4'}")
G=Gate("世代内的年份就是年龄")
p1=G.positive_control("不分层的世代内斜率必须复现 #740(容差 0.0005)",
    planted=float(0.0005-abs(s_gen-0.00672)),floor=0.0,spread=0.00002)
p2=G.negative_control("年龄段内打乱年份后斜率应回到零",null=float(max(abs(q[0]),abs(q[1]))),
    effect=abs(s_age),null_spread=0.00002,
    null_kind="在每个年龄段内部打乱受访者的年份标签 —— 保住年龄构成、每年 n 与作答分布,只毁掉「同一年龄段里谁属于哪一年」")
r=s_age/s_gen
if not p1: v="**UNVERIFIED:旧值不可复现**"
elif q[0]<=s_age<=q[1]: v=f"**W4:年龄段内斜率 {s_age:+.5f} 落在自己的零里 ⇒ 这具仪器分不开**"
elif (s_sex/s_gen)<=0.60 and r<=0.60: v=f"**判不了:年龄分层压到 {r:.0%},而性别安慰剂也压到 {s_sex/s_gen:.0%} ⇒ 压掉的可能是 n**"
elif r>=0.60: v=f"**W1:年龄段内仍有世代内斜率的 {r:.0%} ⇒ 排除纯年龄,两刀合起来剩下年代**"
elif r<=0.25: v=f"**W2:年龄段内只剩 {r:.0%} ⇒ `#740` 的措辞要改成「人老了才连起来」**"
else: v=f"**W3:年龄段内剩 {r:.0%} ⇒ 报份额,不报判决**"
print(f"\n{v}"); print(G)
json.dump(dict(slope_within_cohort=s_gen,slope_within_age=s_age,share=float(r),
  slope_placebo_sex=s_sex,slope_age_x_educ=s_ae,null_ci=[float(q[0]),float(q[1])],
  verdict=v,identification="APC 不可识别:本轮只排除纯年龄与纯队列,不估年代系数",
  unchallenged=True),open(OUT/"ap.json","w"),indent=1,ensure_ascii=False)
