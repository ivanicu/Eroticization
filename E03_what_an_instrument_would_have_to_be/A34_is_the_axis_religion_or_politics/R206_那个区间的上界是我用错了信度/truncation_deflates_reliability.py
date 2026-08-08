"""#767 第二臂 —— 「截断样本会压低量表的表观信度」在第二具仪器上还成不成立?

第一臂在 GSS 上量到:同一把宗教尺子,**截断样本(删掉 never-attenders)α = 0.4804,
完整样本 α = 0.7287** —— 缺陷不只是移动了份额,它把量表的内部一致性砍掉了三分之一。
⚠ 而那是**一次意外**,不是设计出来的对照。**一次意外不是一条规律。** 换仪器重测。

G1 估计量:同一量表在**完整样本**与在**按某一题截断后的样本**上的 α 与 ω 之差。
仪器:MFQ 六题纯洁量表(`disgusting`·`decency`·`god`·`harmlessdg`·`chastity`·`unnatural`)。
⚠ 截断的做法要与 GSS 那次**同构**:GSS 是删掉**某一题的最低档**(attend=never)。
   这里同样删掉**某一题最低档**的人,并**在多个题上各做一次**(G3 全网格,不是一格)。

预注册(判词按 `#764` 新写法:只比已测量的量):
  若**每一格**的截断 α 都低于完整样本 α ⇒ 规律在第二具仪器上重现
  若**符号不一**(有的升有的降)⇒ 那是 GSS 那一次的特殊性,**第一臂的推广要收回**
  ⚠ 无论哪种,都报**整张网格**,包括不同意的格
"""
import pandas as pd, numpy as np, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[3]))
ROOT=pathlib.Path(__file__).resolve().parents[3]
OUT=pathlib.Path(__file__).resolve().parent/"results"; OUT.mkdir(exist_ok=True)
d=pd.read_spss(ROOT/"data/external/dataverse/mfq/GrahamHaidtNosek.2009.JPSP.Study_3.sav",convert_categoricals=False)
PUR=["disgusting","decency","god","harmlessdg","chastity","unnatural"]
PUR=[c for c in PUR if c in d.columns]
D=d[PUR].apply(pd.to_numeric,errors="coerce").dropna()
print(f"=== 硬规则①:六题纯洁量表 n={len(D)} ===")
for c in PUR:
    print(f"  {c:12s} 取值 {sorted(D[c].unique().astype(int))} · 最低档人数 {int((D[c]==D[c].min()).sum())}"
          f" ({(D[c]==D[c].min()).mean()*100:.1f}%)")
z=lambda s:(s-s.mean())/s.std(ddof=1)
def alpha(df):
    k=df.shape[1]; return float(k/(k-1)*(1-df.var(ddof=1).sum()/df.sum(axis=1).var(ddof=1)))
def omega(df):
    """一般 k 题:单因子载荷由主成分近似;k=3 时退化为闭式。报的是 ω_total。"""
    R=df.corr(method="spearman").to_numpy()
    w,v=np.linalg.eigh(R); lam=v[:,-1]*np.sqrt(w[-1])
    if lam.sum()<0: lam=-lam
    u=1-lam**2
    if (u<=0).any(): return None
    return float(lam.sum()**2/(lam.sum()**2+u.sum()))
Z=z(D[PUR])
a_full,w_full=alpha(Z),omega(Z)
print(f"\n=== 完整样本:α = {a_full:.4f} · ω = {w_full:.4f} ===")
print(f"\n=== G3 全网格:逐题删掉最低档的人(与 GSS 那次同构)===")
print(f"  {'删掉谁的最低档':16s}{'n':>7s}{'α':>9s}{'Δα':>9s}{'ω':>9s}{'Δω':>9s}")
rows={}
for c in PUR:
    S=D[D[c]>D[c].min()]
    if len(S)<300: continue
    Zs=z(S[PUR]); a,w=alpha(Zs),omega(Zs)
    rows[c]=dict(n=len(S),alpha=a,d_alpha=a-a_full,omega=w,d_omega=(w-w_full) if w else None)
    print(f"  {c:16s}{len(S):7d}{a:9.4f}{a-a_full:+9.4f}"+(f"{w:9.4f}{w-w_full:+9.4f}" if w else f"{'n/a':>9s}{'n/a':>9s}"))
neg=sum(1 for r in rows.values() if r["d_alpha"]<0)
print(f"\n  α 下降的格数:{neg}/{len(rows)}")
print("\n"+"="*70)
if neg==len(rows):
    v=(f"**规律重现:MFQ 上 {len(rows)}/{len(rows)} 格截断后 α 都下降(幅度 "
       f"{min(r['d_alpha'] for r in rows.values()):+.4f} … {max(r['d_alpha'] for r in rows.values()):+.4f})"
       f" ⇒ 「截断压低表观信度」不是 GSS 那一次的特殊性**")
elif neg==0:
    v=f"**反向:MFQ 上 0/{len(rows)} 格下降 ⇒ 第一臂的推广收回,那是 GSS 那一次的特殊性**"
else:
    v=(f"**符号不一:{neg}/{len(rows)} 格下降,其余上升 ⇒ 「截断压低信度」不是普遍规律,"
       f"要看截断落在量表的哪一端 —— 整张网格已在上面全列**")
print(v)
json.dump(dict(n=len(D),alpha_full=a_full,omega_full=w_full,grid=rows,n_down=neg,verdict=v),
          open(OUT/"truncation.json","w"),ensure_ascii=False,indent=1)
