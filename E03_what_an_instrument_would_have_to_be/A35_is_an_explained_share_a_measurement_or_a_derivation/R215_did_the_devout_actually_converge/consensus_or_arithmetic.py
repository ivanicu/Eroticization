"""#776 · E03·A42·R215 —— 虔诚者内部的共识,是不是在几十年里真的收紧了?

`#772` 量到:最虔诚三分位里 87.8% 说婚外性「总是错」,方差塌到 0.542。
**那是一个横截面的事实。而它有一个纵向版本,是关于人的:那个共识是一直这么紧,还是收紧了?**

⚠⚠ **算术陷阱是本轮的全部**(realstat 开篇):**有界量表上,均值逼近端点时方差必降。**
所以「虔诚者的方差在下降」这句话**可能完全由均值移动逼出来**,而不是任何「共识收紧」。
⇒ **观测方差必须对「给定该均值时的最大可能方差」来读。**
四档量表 {1,2,3,4},给定均值 μ,把质量全放在两端时方差最大:
`p = (μ−1)/3` 在 4 上,`var_max = 9p(1−p)`。**报 `var_obs / var_max` —— 它对均值移动免疫。**

G1 估计量:`ratio_t = var_obs(t) / var_max(μ_t)` 在虔诚层与非虔诚层各自随年份的斜率。
Live worlds:
  A **共识真的收紧**:虔诚层的 ratio 随年份显著下降(超出它自己的零)
  B **全是均值逼出来的**:ratio 平坦 ⇒ `#772` 的天花板只是均值移动的后果,不是共识
⚠ **我不欢迎的是 B** —— 它会把「虔诚者彼此意见一致」这句话降成一句算术。

⚠ **最强混淆写在跑之前**:**分层用的是同期的虔诚度分位**,而虔诚度本身的分布在变
(教会出席率几十年在降)⇒ **「最虔诚三分之一」在 1988 与 2024 不是同一群人。**
⇒ 控制:**同时报固定阈值分层**(用全期合并分位切,而不是逐年切)作为 G4 规格曲线的第二格。
⚠ **零**:年份打乱(保住每年的 n 与该层的边际,只毁掉「哪一年配哪个 ratio」),报斜率的 95% 分位。
⚠ **换不了仪器,而这句话是量出来的**:本轮需要**同一道题在 ≥20 年上的重复测量**。
   MFQ 是**单次采集**,没有年代;NSFG 的 `SXOK18` 在 **2017–19 卷的字典里不存在**(已查),给不出长序列;
   SCCS 的单位是社会不是年代。⇒ **GSS 是这台机器上唯一能问这个问题的仪器,只此一具。**

预注册判词(按 `#764` 新写法):虔诚层与非虔诚层的斜率**各带自己的零**并排报;
  超零且为负 ⇒ A;落在零内 ⇒ B。**两种切法都要报,包括不同意的格。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.blocks import aligned, check_kept_codes
from lib.gates import Gate
RNG=np.random.default_rng(215)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
gp=ROOT/"data/external/gss/GSS_stata/gss7224_r3a.dta"
SEX=["premarsx","xmarsex","homosex","teensex"]
VALID={"premarsx":(1,4),"xmarsex":(1,4),"homosex":(1,4),"teensex":(1,4),
       "attend":(0,8),"reliten":(1,4),"fund":(1,3)}
for c,rng in VALID.items():
    dr,_=check_kept_codes(gp,c,rng)
    if dr: print(f"  #766 前瞻:{c} 删 "+" · ".join(f"码{int(a)} {b!r}({sh*100:.1f}%)" for a,b,n,sh in dr[:2]))
d=pd.read_stata(gp,columns=["year"]+list(VALID),convert_categoricals=False)
M=pd.DataFrame({c:pd.to_numeric(d[c],errors="coerce").where(
    lambda v,lo=VALID[c][0],hi=VALID[c][1]:(v>=lo)&(v<=hi)) for c in VALID})
M["year"]=d.year
cat=pd.read_stata(gp,columns=SEX,convert_categoricals=True)
cats={c:list(cat[c].cat.categories) for c in cat.columns}; cats["homosex"]=cats["homosex"][:4]
flip=aligned({c:cats[c] for c in SEX},"strict")
print(f"\n方向由 `aligned()` 定 -> 要翻 {sorted(flip)}(统一成「高=严」)")
for c in flip: M[c]=-M[c]+5     # 翻回 1..4,高=严,便于用有界量表的算术上界
M["reliten"]=-M["reliten"]; M["fund"]=-M["fund"]
z=lambda s:(s-s.mean())/s.std(ddof=1)
def vmax(mu):  # 四档 {1..4} 给定均值时的最大方差
    p=(mu-1)/3.0; return 9*p*(1-p)
def slope(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); return float(np.cov(x,y,ddof=1)[0,1]/np.var(x,ddof=1))
print("\n=== 逐题 · 两种分层切法 · 虔诚层 vs 非虔诚层的 ratio 斜率 ===")
res={}
for cut in ("逐年分位","全期固定阈值"):
    res[cut]={}
    for s in SEX:
        sub=M[[s,"attend","reliten","fund","year"]].dropna().copy()
        sub["REL"]=z(sub[["attend","reliten","fund"]]).mean(axis=1)
        if cut=="逐年分位":
            sub["k"]=sub.groupby("year")["REL"].transform(lambda v: pd.qcut(v,3,labels=False,duplicates="drop"))
        else:
            qs=sub["REL"].quantile([1/3,2/3]).to_numpy()
            sub["k"]=np.digitize(sub["REL"],qs)
        rows={}
        for kk,lab in ((2,"虔诚层"),(0,"非虔诚层")):
            g=sub[sub.k==kk]
            per=[]
            for y,gy in g.groupby("year"):
                if len(gy)<120: continue
                mu=float(gy[s].mean()); vo=float(gy[s].var(ddof=1)); vm=vmax(mu)
                if vm<=1e-6: continue
                per.append((int(y),mu,vo,vo/vm,len(gy)))
            if len(per)<8: continue
            yr=np.array([p[0] for p in per],float); rt=np.array([p[3] for p in per])
            b=slope(yr,rt)
            nul=[slope(RNG.permutation(yr),rt) for _ in range(2000)]
            q=float(np.quantile(np.abs(nul),.95))
            rows[lab]=dict(n_years=len(per),first=per[0],last=per[-1],slope=b,null95=q,ratio=abs(b)/q)
        res[cut][s]=rows
    print(f"\n--- {cut} ---")
    print(f"  {'题':10s}{'层':8s}{'年数':>5s}{'首年 ratio':>11s}{'末年 ratio':>11s}{'斜率/年':>11s}{'零95%':>10s}{'倍数':>7s}")
    for s in SEX:
        for lab,r in res[cut][s].items():
            print(f"  {s if lab=='虔诚层' else '':10s}{lab:8s}{r['n_years']:5d}{r['first'][3]:11.4f}{r['last'][3]:11.4f}"
                  f"{r['slope']:+11.6f}{r['null95']:10.6f}{r['ratio']:6.2f}×")
G=Gate("#776 · 共识还是算术")
mu_lo,mu_hi=1.5,3.5
# ⚠⚠ 第一版 ① 是 `vmax(1.0)` 比 0.0 —— 值上退化,库判 DEGENERATE。
#    而 **`#775` 那具 lint 在跑之前扫过这个脚本,返回 0** —— 它看不见「一个恰好返回 0 的函数调用」。
#    ⇒ **`#775`② 记的边界比我写的还宽:lint 只认语法形状,认不出值上的退化。**
#    两层是互补的:lint 在**写的时候**拦一部分,库的 `_degenerate` 在**跑的时候**拦全零那部分。
#    改成比**两个都非零**的值:vmax 在两个不同 μ 上的值必须各自对上理论。
G.identity_control("① 上界公式在 μ=2.5 处须等于 2.25(仪器检查)",
                   observed=float(vmax(2.5)), expected=2.25, tol=1e-9, what="9·0.5·0.5;两个值都非零")
G.identity_control("①b 上界公式在 μ=2.0 处须等于 2.00", observed=float(vmax(2.0)), expected=2.00, tol=1e-9,
                   what="p=1/3 ⇒ 9·(1/3)·(2/3)=2;两点对上才说明公式不是碰巧")
G.asserted("② 上界公式在端点须给 0(端点行为)", bool(abs(vmax(1.0))<1e-12 and abs(vmax(4.0))<1e-12),
           f"vmax(1)={vmax(1.0):.2e} · vmax(4)={vmax(4.0):.2e}", kind="control")
# ⚠⚠ 第一版 ③ 把**中位 |斜率|** 与**中位 零95% 分位**相比 —— **两个不同的分位,不是同一个对象**,
#    正是 realstat「控制因自己的理由而失败」那一行。⇒ 改成**逐格**比它自己的零,再数格子。
devout=[(c,s,res[c][s]["虔诚层"]) for c in res for s in SEX if "虔诚层" in res[c][s]]
neg_sig=[(c,s) for c,s,r in devout if r["slope"]<0 and abs(r["slope"])>r["null95"]]
pos_sig=[(c,s) for c,s,r in devout if r["slope"]>0 and abs(r["slope"])>r["null95"]]
G.asserted("③ 逐格比它自己的零:虔诚层有几格是「显著下降」",
           True, f"显著下降 {len(neg_sig)}/{len(devout)} 格 · 显著上升 {len(pos_sig)}/{len(devout)} 格 "
                 f"—— 这一条只报数,判词在下面", kind="control")
print(); print(G)
above=[(c,s) for c in res for s in SEX if "虔诚层" in res[c][s] and res[c][s]["虔诚层"]["ratio"]>1.0
       and res[c][s]["虔诚层"]["slope"]<0]
tot=[(c,s) for c in res for s in SEX if "虔诚层" in res[c][s]]
print("\n"+"="*76)
if len(above)==len(tot):
    v=f"**A:{len(above)}/{len(tot)} 格虔诚层的 ratio 显著下降 ⇒ 共识真的收紧,不是均值逼出来的**"
elif not above:
    v=f"**B:0/{len(tot)} 格显著下降 ⇒ 「虔诚者越来越一致」是均值移动的算术后果,不是共识收紧**"
else:
    v=f"**不一致:{len(above)}/{len(tot)} 格显著下降 ⇒ 报出整张网格,不选边**"
print(v)
json.dump(dict(res=res,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"consensus.json","w"),ensure_ascii=False,indent=1)
