"""#762 第二臂 —— 同一个问题换一具人层仪器:宗教能不能解释掉 chastity×权威?

GSS 臂:偏掉三个宗教量后 `obey`×四题性保留 **66.9%**,再加政治 62.4%。
⚠ 闸(`#658`)要求同一个问题在 ≥2 具仪器上问过 —— 而 MFQ 自带宗教出席变量,**换得动就不许写豁免**。

G1 估计量:`ρ(chastity, AUTHORITY_AVG | 控制量)`,控制量 = 空 · 政治 · 宗教出席 · 两者。
⚠ 预测(从 GSS 迁移过来,可证伪):**宗教解释的份额应当大于政治**,而两者都不该解释掉大半。
⚠ 而 MFQ 是自选样本、政治左偏 ⇒ **两具仪器不一致时,以概率样本(GSS)为准**,与 `#760` 同一条预注册。
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
from lib.gates import Gate
RNG=np.random.default_rng(1202)
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
print("=== 硬规则①:变量名不是测量 ===")
CAND=["religatt_num","Religion_attend_num","religion_attend","religion_current"]
for c in CAND:
    if c in d.columns:
        v=pd.to_numeric(d[c],errors="coerce")
        print(f"  {c:22s} n={int(v.notna().sum()):5d} 取值 [{v.min():.0f},{v.max():.0f}] sd={v.std():.3f}")
    else: print(f"  {c:22s} 不存在")
REL=[c for c in ("religatt_num",) if c in d.columns]
assert REL, "MFQ 没有可用的宗教出席量 —— 那才该写豁免"
M=d[["chastity","AUTHORITY_AVG","FAIRNESS_AVG","politics"]+REL].apply(pd.to_numeric,errors="coerce").dropna()
print(f"\n完整个案 n={len(M)}")
def resid(y,X):
    X=np.c_[np.ones(len(X)),X]; return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def prho(a,b,c=None):
    r=lambda v: pd.Series(v).rank().to_numpy(float)
    if c is None: return float(np.corrcoef(r(a),r(b))[0,1])
    C=np.asarray(c,float); C=C.reshape(-1,1) if C.ndim==1 else C
    rc=np.column_stack([r(C[:,j]) for j in range(C.shape[1])])
    return float(np.corrcoef(resid(r(a),rc),resid(r(b),rc))[0,1])
print(f"  ⚠ 共线性预检:ρ(chastity, 宗教出席) = {prho(M.chastity,M[REL[0]]):+.4f} · "
      f"ρ(chastity, politics) = {prho(M.chastity,M.politics):+.4f} · "
      f"ρ(宗教出席, politics) = {prho(M[REL[0]],M.politics):+.4f}")
SPECS={"(空,偏前)":[], "政治 politics":["politics"], "宗教出席":REL, "宗教 + 政治":REL+["politics"]}
out={}
for tgt in ("AUTHORITY_AVG","FAIRNESS_AVG"):
    raw=prho(M.chastity.to_numpy(),M[tgt].to_numpy())
    floor=3*1.65/np.sqrt(len(M))
    print(f"\n=== chastity × {tgt}(偏前 {raw:+.4f};比值分母下限 {floor:.4f})===")
    row={}
    for nm,cols in SPECS.items():
        v=raw if not cols else prho(M.chastity.to_numpy(),M[tgt].to_numpy(),M[cols].to_numpy())
        keep=(v/raw) if abs(raw)>=floor else None
        row[nm]=dict(val=v,keep=keep)
        print(f"  {nm:16s} {v:+.4f}  保留 "+(f"{keep*100:6.1f}%" if keep is not None else "不可读(分母近零)"))
    out[tgt]=row
raw=prho(M.chastity.to_numpy(),M.AUTHORITY_AVG.to_numpy())
pc=prho(M.chastity.to_numpy(),M.AUTHORITY_AVG.to_numpy(),RNG.normal(0,1e-9,(len(M),1)))
nul=[prho(M.chastity.to_numpy(),M.AUTHORITY_AVG.to_numpy(),
          M[REL].to_numpy()[RNG.permutation(len(M))]) for _ in range(300)]
G=Gate("#762 第二臂 · MFQ")
G.identity_control("① 常数控制须回到偏前", observed=pc, expected=raw, tol=0.005, what="仪器活着吗")
G.identity_control("② 打乱宗教出席须回到偏前", observed=float(np.median(nul)), expected=raw, tol=0.005,
                   what="打乱只毁掉「谁的宗教配谁的道德」")
G.offset_control("③ 偏掉真宗教出席后须显著低于偏前", effect=out["AUTHORITY_AVG"]["宗教出席"]["val"],
                 offset=raw, spread=float(np.std(nul)),
                 null_kind="同一批人、同一对量,唯一差别是控制量是真的宗教出席还是打乱后的宗教出席")
print(); print(G)
gk=out["AUTHORITY_AVG"]["宗教出席"]["keep"]; pk=out["AUTHORITY_AVG"]["政治 politics"]["keep"]
print("\n"+"="*66)
print(f"**MFQ 臂:宗教解释 {(1-gk)*100:.1f}% · 政治解释 {(1-pk)*100:.1f}%;"
      f"两者全放后仍保留 {out['AUTHORITY_AVG']['宗教 + 政治']['keep']*100:.1f}%**")
print(f"  ⚠ 对照 GSS 臂:宗教 33.1% · 政治 9.4% · 全放后保留 62.4% —— "
      f"{'两具仪器同向' if (1-gk)>(1-pk) else '**两具仪器不同向 —— 框架就是发现**'}")
json.dump(dict(n=len(M),specs=out,gate_ok=all(r[2] for r in G.rows),
               gss_ref=dict(religion=0.331,politics=0.094,残余=0.624)),
          open(OUT/"mfq_religion_arm.json","w"),ensure_ascii=False,indent=1)
