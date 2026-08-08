"""#773 第二臂 —— 「性态度题在虔诚层饱和」这条能力边界,换一具仪器还成不成立?

本轮的判词是 UNVERIFIED,而唯一能报的内容是一条**能力边界**:
**GSS 的四档性态度题在最虔诚层饱和到 61–88%,该层的关系强度读数不可读。**
⚠ 一条能力边界若只在一具仪器上量过,它可能是那具仪器的毛病,**而不是这类题的毛病** ——
所以它必须换仪器。NSFG 2011–13 女性卷有 `sxok18`/`sxok16`(未婚 18/16 岁发生性行为可不可以)
与 `attndnow`(当前礼拜出席),**不同问卷、不同编码团队、不同人群。**

⚠ **它没有 `obey` ⇒ 不能复现交互**,只能复现**饱和**这一件事。如实标为**部分复现**。
G1 估计量:各礼拜出席层内,态度题落在**最严端点**的占比。
预注册:**最虔诚层端点占比 ≥60% ⇒ 这类题在虔诚层饱和,不是 GSS 的毛病;
<60% ⇒ 边界收窄为 GSS 特有,页面上要写明。**
"""
import pandas as pd, numpy as np, re, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
NS=ROOT/"data/external/nsfg"
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NS/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    m=pat.search(line)
    if m: LAY[m.group(2).lower()]=(int(m.group(1))-1,int(m.group(3)),m.group(4))
need=[c for c in ("sxok18","sxok16","attndnow","ager") if c in LAY]
print("=== 硬规则①:变量名不是测量 ===")
for c in need: print(f"  {c:10s} 列位 {LAY[c][0]+1} 宽 {LAY[c][1]}  {LAY[c][2][:58]}")
buf={n:[] for n in need}
for line in open(NS/"2011_2013_FemRespData.dat",errors="replace"):
    for n in need:
        s,w,_=LAY[n]; v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame(buf)
for c in ("sxok18","sxok16"): X[c]=X[c].where(X[c].between(1,5))
# ⚠⚠ 三个真错,全被本臂的分层控制拦下,记在这里:
#   ① `attndnow` 有**码 7(n=1,378)**而我第一版截到 1–6 —— **丢掉了最大的一档**;
#   ② **方向是反的**:mean(sxok18) 从码 1 的 3.204 单调降到码 7 的 2.216
#      ⇒ **码 1 = 出席最多**,而我把「最虔诚层」当成了码 6;
#   ③ `sxok18` 的最严档(码 5)**全样本只有 48 人(0.8%)** —— 这是 agree/disagree 量表,
#      严端本来就罕见,**与 GSS 的「总是错」不是同一种题**。
X["attndnow"]=X["attndnow"].where(X["attndnow"].between(1,7))
print(f"\n原始行数 {len(X)} · attndnow 取值分布 {X.attndnow.value_counts().sort_index().to_dict()}")
print(f"  ⚠ sxok18 取值 {sorted(X.sxok18.dropna().unique().astype(int))} —— 1=strongly agree … 5=strongly disagree")
print(f"     ⇒ **最严端点是 5(强烈不同意「未婚 18 岁发生性关系没问题」)**")
res={}
for c in ("sxok18","sxok16"):
    sub=X[[c,"attndnow"]].dropna()
    hi=sub[c].max()
    rows=[]
    for k in sorted(sub.attndnow.unique()):
        g=sub[sub.attndnow==k]
        if len(g)<200: continue
        rows.append(dict(level=int(k),n=len(g),top=float((g[c]==hi).mean()),sd=float(g[c].std(ddof=1))))
    res[c]=rows
    print(f"\n=== {c} · 按礼拜出席分层(端点 = 最严的一档)===")
    for r in rows: print(f"  出席码 {r['level']}  n={r['n']:5d}  最严端点 {r['top']*100:5.1f}%  sd {r['sd']:.3f}")
mx=max(r["top"] for c in res for r in res[c])
# ⚠ 出席最多的是**码 1**(由 mean 单调性确定,不是由码序假设)
most=min(res["sxok18"],key=lambda r:r["level"])
G=Gate("#773 第二臂 · NSFG 饱和")
G.asserted("① 分层必须真的分开了(最虔诚层与最不虔诚层的端点占比须不同)",
           bool(abs(most["top"]-max(res["sxok18"],key=lambda r:r["level"])["top"])>0.05),
           f"出席最多层(码1) {most['top']*100:.1f}% vs 出席最少层(码7) {max(res['sxok18'],key=lambda r:r['level'])['top']*100:.1f}%",
           kind="control")
G.asserted("② 预注册:最虔诚层端点占比 ≥60% ⇒ 这类题在虔诚层饱和",
           bool(most["top"]>=0.60), f"最高出席层端点占比 {most['top']*100:.1f}%(阈值 60%)", kind="kill")
print(); print(G)
print("\n"+"="*70)
# ⚠⚠ 而修完之后,真正的问题不是「没重现」,是**这一臂测的不是同一个量**:
#    GSS 的题是「总是错 … 完全不错」,严端是虔诚者的众数;
#    NSFG 的是 agree/disagree,严端全样本 0.8%。**端点占比在两种题上不是同一个统计量。**
if most["top"]>=0.60:
    v=f"**重现:出席最多层端点占比 {most['top']*100:.1f}% ≥ 60%**"
else:
    v=(f"**UNVERIFIED —— 仪器不匹配,不是「没重现」**:出席最多层端点占比只有 {most['top']*100:.1f}%,"
       f"但 `sxok18` 的严端全样本仅 0.8%,**它是 agree/disagree 量表,严端对所有人都罕见**;"
       f"GSS 的「总是错」则是虔诚者的众数。**端点占比在两种题上不是同一个统计量,这一臂答不了那个问题。**\n"
       f"  ⇒ **`#772` 的能力边界维持为 GSS 特有,页面上已写明。**")
print(v)
print("⚠ 本臂**不能**复现交互 —— NSFG 没有 `obey`;这是**部分复现**。")
json.dump(dict(res=res,max_top=mx,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"nsfg_ceiling.json","w"),ensure_ascii=False,indent=1)
