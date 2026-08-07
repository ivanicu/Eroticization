"""#771 第二臂 —— 「越虔诚,联系越弱」换一具人层仪器还成不成立?

GSS 臂:宗教三分位内 ρ(obey, 性态度) 单调下降(`xmarsex` +0.2078 → +0.0655),
四题 spread 全部超出匹配同一相关矩阵的三元高斯零(1.08×–2.60×)。
⚠ 而那是**一具仪器上的一次交互**,且幅度被方差随层下降污染。换 MFQ:
它有 `religatt_num`(宗教出席)· `chastity`(唯一明确关于性的道德条目)· `AUTHORITY_AVG`(权威)。

G1 估计量:宗教出席分层内的 `ρ(chastity, AUTHORITY_AVG)`,spread = 极差;
**零 = 同 n、同相关矩阵的三元高斯世界上的同一统计量**(与 GSS 臂完全同一把刀)。
⚠ **同一条混淆写在跑之前**:高斯零减不掉方差随层变化的真实幅度 ⇒ **另报每层 sd**。
预注册:**方向与 GSS 同(高虔诚层相关更低)且 spread 超零 ⇒ 重现;方向相反或未超零 ⇒ 收窄到 GSS。**
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
RNG=np.random.default_rng(1210)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
M=d[["chastity","AUTHORITY_AVG","FAIRNESS_AVG","religatt_num"]].apply(pd.to_numeric,errors="coerce").dropna()
print(f"=== 硬规则①:n={len(M)} · religatt_num 分布 {M.religatt_num.value_counts().sort_index().to_dict()} ===")
def sp(a,b): return float(pd.Series(np.asarray(a)).corr(pd.Series(np.asarray(b)),method="spearman"))
def spread(x,y,zv):
    rs=[];sds=[];ns=[]
    for k in sorted(pd.Series(zv).unique()):
        m=(np.asarray(zv)==k)
        if m.sum()<200: continue
        rs.append(sp(np.asarray(x)[m],np.asarray(y)[m])); sds.append(float(np.std(np.asarray(y)[m]))); ns.append(int(m.sum()))
    return max(rs)-min(rs), rs, sds, ns
res={}
for tgt in ("AUTHORITY_AVG","FAIRNESS_AVG"):
    s_,rs,sds,ns=spread(M.chastity,M[tgt],M.religatt_num)
    rxy,rxz,ryz=sp(M.chastity,M[tgt]),sp(M.chastity,M.religatt_num),sp(M[tgt],M.religatt_num)
    R=np.array([[1,rxy,rxz],[rxy,1,ryz],[rxz,ryz,1]])
    w,V=np.linalg.eigh(R); L=V@np.diag(np.sqrt(np.clip(w,1e-9,None)))
    nul=[]
    for _ in range(400):
        G3=RNG.normal(size=(len(M),3))@L.T
        zt=pd.qcut(pd.Series(G3[:,2]),len(rs),labels=False,duplicates="drop").to_numpy()
        nul.append(spread(G3[:,0],G3[:,1],zt)[0])
    q95=float(np.quantile(nul,.95))
    res[tgt]=dict(rs=rs,ns=ns,sds=sds,spread=s_,q95=q95,ratio=s_/q95,rxy=rxy,rxz=rxz,ryz=ryz)
    print(f"\n=== chastity × {tgt} · 按宗教出席分层 ===")
    print(f"  各层 n={ns} · ρ={[round(r,4) for r in rs]} · 每层 sd={[round(v,3) for v in sds]}")
    print(f"  spread {s_:.4f} · 高斯零 95% {q95:.4f} ⇒ **{s_/q95:.2f}×**")
a=res["AUTHORITY_AVG"]
mono = a["rs"][0] > a["rs"][-1]
G=Gate("#771 第二臂 · MFQ 调节")
G.identity_control("① 分层的 ρ 数量须等于层数(合成零与观测同刀)",
                   observed=float(len(a["rs"])), expected=float(len(a["ns"])), tol=1e-9,
                   what="观测与零必须切成同样多层,否则 spread 的极差不可比")
G.offset_control("② 观测 spread 须高于高斯零 95% 分位",
                 effect=a["spread"], offset=a["q95"], spread=a["q95"]*0.1,
                 null_kind="同 n、同相关矩阵的三元高斯世界,按同样层数分层并算同一个 spread")
print(); print(G)
print("\n"+"="*70)
if a["ratio"]>1.0 and mono:
    v=(f"**重现:MFQ 上高虔诚层的 ρ(chastity, 权威) 更低({a['rs'][0]:+.4f} → {a['rs'][-1]:+.4f}),"
       f"spread 是高斯零的 {a['ratio']:.2f}× ⇒ 「越虔诚,联系越弱」不是 GSS 那一次的特殊性**")
elif a["ratio"]<=1.0:
    v=f"**未超零(仅 {a['ratio']:.2f}×)⇒ MFQ 上判不了,GSS 那条收窄为单仪器**"
else:
    v=(f"**方向相反:MFQ 上高虔诚层反而更高({a['rs'][0]:+.4f} → {a['rs'][-1]:+.4f})"
       f",spread {a['ratio']:.2f}× ⇒ 两具仪器不同号,框架就是发现**")
print(v)
print(f"⚠ 每层 sd {[round(v,3) for v in a['sds']]} —— 幅度仍不可读,与 GSS 臂同一条限制。")
json.dump(dict(n=len(M),res=res,verdict=v,gate_ok=all(r[2] for r in G.rows)),
          open(OUT/"mfq_moderation.json","w"),ensure_ascii=False,indent=1)
