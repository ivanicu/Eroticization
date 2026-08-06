"""E03·A19·R102 —— 那 0.57 是不是 2021 年那次换问法造成的

**类型:FRONTIER**。`#659` 的安慰剂失败(+0.5697 = 主效应 167%),指认了一个**年份层共同成分**,
而它有一个具体嫌疑人:**2021 年 GSS 改成网络作答**。

⚠ **`#111c`:这是这一问题的第一次修正后重试。**
⚠ **BASIN**:我连着几轮都在发现「我自己的仪器是问题」——**本轮下注 W1(就是 2021,可修)**,
   即**下注一个方便的结果**,这是对最近那条故事的反向下注。

## 硬规则①(已跑)
GSS **有** `mode` 列(`Interview done in-person or over the phone`),取值 1–4。
**2021 年 mode 1(面访)= 0**,mode 4 = 3521;而 **2022 与 2024 两种方式同年并存**。
⇒ **不必删波去猜:可以在同一年内直接量。**

W1 **就是 2021** —— 剔波后安慰剂掉到主量一半以下 ⇒ 主量可读。
W2 **不是 2021** —— 安慰剂仍高 ⇒ **年代单位在这具仪器上判不了**,写进「做不到」并停。
W3 **方式效应可直接量且很大** —— 同年内 mode 1 vs mode 4 的谴责份额差,**这是一句关于人的话**:
   换一种问法,人对性的道德判断会移动多少。

## G1 ESTIMAND(先于方法)
**① 预注册的那个**:剔除 2021+2022 两波后,重算 ⑴ 主量(发生率差分 × 谴责差分)⑵ 安慰剂
   (谴责差分 × 死刑差分)⑶ 正对照(premarsx 差分 × homosex 差分)。
**② 标注的追加**(`mode` 列存在使它成为可能,而预注册⑤只写了「找不到就用年份哑变量」):
   **同年内 mode 1 与 mode 4 的谴责份额差**,在 2022 与 2024 各算一次 —— **直接量,不是推断。**

## KILL(条件式,预注册)
if 正对照 > 0:
  安慰剂 |r| < 0.5 × |主量| -> **W1**,然后才读主量
  否则 -> **W2:年代单位在这具仪器上判不了,停**
else UNVERIFIED
⚠ **最强混淆先写死**:剔两波只剩 **14 个差分**,**MDE 从 0.651 升到约 0.70** ——
  **所以这一轮只可能推翻,不可能确立**;主量若变得不显著,**那是功率丢的,不是效应没的**,必须这样写。

## IMPOSSIBLE(不写 planned)
14 个差分 ⇒ 功率极低是结构性的 · 无干预 · 自报 ·
**跨仪器:换不了仪器,只此一具**(`#659` 已验证:只有 GSS 同时具备两侧且波数 ≥10)· `[unchallenged]`
"""
import os, sys, pathlib, json, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from scipy.stats import pearsonr, norm
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
P="data/external/gss/GSS_stata/gss7224_r3a.dta"
df,_=pyreadstat.read_dta(P, usecols=["year","evstray","xmarsex","premarsx","homosex","cappun","mode"], encoding="latin1")

def share(col,val,d=None,floor=200):
    d=df if d is None else d
    g=d[d[col].isin([1,2,3,4])].groupby("year")
    s=g.apply(lambda x:(x[col]==val).mean()); n=g.size()
    return s[n>=floor]
def prevalence(d=None,floor=200):
    d=df if d is None else d
    g=d[d.evstray.isin([1,2])].groupby("year")
    s=g.apply(lambda x:(x.evstray==1).mean()); n=g.size()
    return s[n>=floor]
def corr(a,b):
    y=sorted(set(a.index)&set(b.index))
    if len(y)<5: return np.nan,0
    return float(pearsonr(np.diff(a[y].values),np.diff(b[y].values)).statistic), len(y)-1

print("=== ① 预注册:剔除 2021 + 2022 ===")
res={}
for tag,d in [("全部年份",df),("剔 2021+2022",df[~df.year.isin([2021,2022])])]:
    prev=prevalence(d); cond=share("xmarsex",1,d); cap=share("cappun",1,d)
    pm,hm=share("premarsx",1,d),share("homosex",1,d)
    main,n1=corr(prev,cond); plac,n2=corr(cond,cap); pos,n3=corr(pm,hm)
    mde=float(np.tanh((norm.ppf(0.975)+norm.ppf(0.80))/np.sqrt(max(n1-3,1))))
    res[tag]=dict(main=main,n_main=n1,placebo=plac,n_plac=n2,positive=pos,n_pos=n3,mde=mde)
    print(f"  {tag:14s} 主量 {main:+.4f}(n={n1}) · 安慰剂 {plac:+.4f}(n={n2}) · 正对照 {pos:+.4f}(n={n3}) · MDE {mde:.4f}")
A=res["剔 2021+2022"]
ratio=abs(A["placebo"])/max(abs(A["main"]),1e-9)
print(f"  ⇒ 剔波后 安慰剂/主量 = **{ratio:.2f}**(判据 < 0.50)")

print("\n=== ② 标注的追加:同年内直接量方式效应(不是推断)===")
mode_rows=[]
for y in [2022,2024]:
    d=df[(df.year==y)&df.xmarsex.isin([1,2,3,4])]
    a=d[d["mode"]==1]; b=d[d["mode"]==4]
    if len(a)<200 or len(b)<200: print(f"  {y}: n 不足({len(a)}/{len(b)})"); continue
    sa=(a.xmarsex==1).mean(); sb=(b.xmarsex==1).mean()
    se=np.sqrt(sa*(1-sa)/len(a)+sb*(1-sb)/len(b))
    mode_rows.append(dict(year=int(y),n_face=int(len(a)),n_web=int(len(b)),face=float(sa),web=float(sb),
                          diff=float(sa-sb),ci=[float(sa-sb-1.96*se),float(sa-sb+1.96*se)]))
    print(f"  {y}  面访 n={len(a):5d} 谴责 {sa:.4f} · 网络 n={len(b):5d} 谴责 {sb:.4f} · "
          f"**差 {sa-sb:+.4f}** 95%CI [{sa-sb-1.96*se:+.4f}, {sa-sb+1.96*se:+.4f}]")
# 安慰剂:同一对比作用在与性无关的道德题上
for y in [2022,2024]:
    d=df[(df.year==y)&df.cappun.isin([1,2])]
    a=d[d["mode"]==1]; b=d[d["mode"]==4]
    if len(a)<200 or len(b)<200: continue
    sa=(a.cappun==1).mean(); sb=(b.cappun==1).mean()
    print(f"  {y}  安慰剂(死刑)面访 {sa:.4f} · 网络 {sb:.4f} · 差 {sa-sb:+.4f}")

G=Gate("那 0.57 是不是 2021 年那次换问法造成的")
p1=G.positive_control("剔波后 premarsx×homosex 仍须为正",planted=A["positive"],floor=0.0,spread=0.05)
p2=G.negative_control("剔波后的安慰剂(谴责×死刑)必须掉下来",null=abs(A["placebo"]),effect=abs(A["main"]),
                      null_spread=0.05,null_kind="与性无关的道德态度序列,同样做一阶差分")
if p1 and p2:
    verdict=(f"**W1 —— 就是 2021:剔波后安慰剂 {A['placebo']:+.4f} 掉到主量 {A['main']:+.4f} 的一半以下,"
             f"主量可读。⚠ 而 MDE 升到 {A['mde']:.3f},只可能推翻,不可能确立。**")
elif p1:
    verdict=(f"**W2 —— 不是 2021:剔波后安慰剂仍是主量的 {ratio:.2f} 倍。"
             f"⇒ 年代单位在这具仪器上判不了,写进「做不到」并停。**")
else: verdict="UNVERIFIED —— 正对照失败"
print(f"\n{verdict}"); print(G)
json.dump(dict(preregistered=res,ratio=ratio,mode_effect=mode_rows,verdict=verdict,unchallenged=True),
          open(OUT/"mode_or_structure.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'mode_or_structure.json'}")
