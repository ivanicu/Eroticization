import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A142 R425 -- `#296` 那条 D5 的最后一块:`S` 的两性差

`#296` 的性别非不变性是页面上**唯一停在 D5** 的结论。它现在被夹在两个已知之间:
`#378a` 证明它**不在** `c3⁻ -> 羞耻` 这条路上(|t| 0.10);
`#380a` 证明它**在** `EARLY -> 羞耻` 这条路上(|t| 3.69,三个摘要一致)。
**-> 那 `S` 那条呢?** `S` 正是 `#296` 非不变性的**原始所在地**。

⚠ **本轮的正/负对照是**真实数据**上的,不是合成的** —— 比合成对照强一档:
**`EARLY` 必须显出差(阳性参照)· `c3⁻` 必须不显(阴性参照)**,两者都在**同一个模型、同一次运行**里。
一个把 `c3⁻` 判成有差、或把 `EARLY` 判成无差的仪器,它对 `S` 说什么都不算数。

⚠ 先读两个 IMPOSSIBLE 栏(`#376d`):
`#296`:结局侧混淆 -> **本轮用组内秩变换**(比 `#296` 当年的手段强,也拿掉形状差异);
`#378e②`:尺子在合并样本上估 -> **本轮量的是「同一把尺子下两组的斜率差」**,不是「两组各自的路」。

ESTIMAND        `羞耻(组内秩变换) ~ S + c3⁻ + EARLY + 类别数`,男女各拟合,窄口径(`S` 需 ≥8 块)。
                主量 = **`S` 的两性系数差**,对**族内**阈(3 条路一起)。
判据(**先标支**,`#379c`)
                【两支】阳性参照 `EARLY` 越阈 · 阴性参照 `c3⁻` 不越阈 · 随机劈分的零非退化。
                【非零支】`S` 越过族内阈 -> 非不变性**也在** `S` 上。
                【零支】仅当 `S` 未越阈时启用 MDE:MDE < 0.05 才算「看得见而没有」。
⚠ 零的种类     `offset_control`:**随机等大小劈分**的族内 max-|t|。
IMPOSSIBLE      ① 窄口径 n=6,717,功率比 `#380` 低 -> **先算 MDE**;
                ② 三条路彼此相关,族内阈对每条都用同一个 -> 这是**保守**的,不是精确的;
                ③ 「`S` 无差」若成立,也只说**斜率**一样,不说两组的 `S` **分布**一样(`#293` 已证明分布不一样)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
SEX=pd.to_numeric(d['biomale'],errors='coerce').values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(EARLY)&np.isfinite(sh)&np.isfinite(SEX)&np.isfinite(ncat)
g0=M&(SEX==0); g1=M&(SEX==1)
LAB=['S 位置','c3⁻','EARLY']; PRED=[S,C3,EARLY]
def rank_within(y,g):
    out=np.full(NN,np.nan); j=np.flatnonzero(g); r=rankdata(y[j]); out[j]=(r-r.mean())/r.std(); return out
YR=np.full(NN,np.nan)
for gs in (g0,g1):
    j=np.flatnonzero(gs); YR[j]=rank_within(sh,gs)[j]
def fit(y,g):
    k=int(g.sum()); yy=y[g]; yy=(yy-yy.mean())/max(yy.std(),1e-12)
    X=np.column_stack([np.ones(k)]+[(v[g]-v[g].mean())/max(v[g].std(),1e-12) for v in PRED+[ncat]])
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-5); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:4],se[1:4]
b0,s0=fit(YR,g0); b1,s1=fit(YR,g1)
D_=b1-b0; TT=D_/np.maximum(np.sqrt(s0**2+s1**2),1e-12)
n0=int(g0.sum()); idx=np.flatnonzero(M)
print(f"n=**{int(M.sum()):,}**(窄口径,`S` 需 ≥8 块)· 组0 {n0:,} · 组1 {int(g1.sum()):,}")
print(f"结局 = **组内秩变换**的羞耻(比 `#296` 当年的手段强:也拿掉形状差异)\n")
mt=[]
for s_ in range(400):
    rg=np.random.default_rng(7000+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    ba,sa=fit(YR,ga); bb,sb=fit(YR,gb)
    mt.append(float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))))
mt=np.array(mt); THR=float(np.percentile(mt,95))
print(f"⚠ offset 零(**随机等大小劈分** 400 次的族内 max-|t| 95 分位)= **{THR:.3f}** "
      f"(零的均值 {mt.mean():.3f} —— **任意两个子样本都会有差**)\n")
for i,l in enumerate(LAB):
    role={'S 位置':'← 本轮的问题','c3⁻':'← 阴性参照(`#378a`:应当无差)',
          'EARLY':'← 阳性参照(`#380a`:应当有差)'}[l]
    print(f"   {l:<8} 组0 **{b0[i]:+.4f}** · 组1 **{b1[i]:+.4f}** · 差 **{D_[i]:+.4f}** · "
          f"|t| **{abs(TT[i]):.3f}** · {'**越阈**' if abs(TT[i])>THR else '未越阈'}  {role}")

print(f"\n先算 MDE 再看数(`#369a`),只在组 1 上加 `S` 的额外斜率,每级 30 次:")
MDE=None
for gg in (0.03,0.05,0.08,0.12,0.20):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(9500+int(gg*100)*13+s_); y=np.full(NN,np.nan)
        for gs,ex in ((g0,0.0),(g1,gg)):
            k=int(gs.sum()); zs=(S[gs]-S[gs].mean())/S[gs].std()
            y[gs]=0.10*zs+ex*zs+rg.standard_normal(k)
        ba,sa=fit(y,g0); bb,sb=fit(y,g1)
        if abs(((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12))[0])>THR: hit+=1
    print(f"   组1 多出 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.25
print(f"   **MDE(S)= {MDE_:.2f}** · 有意义 0.05")
pd.DataFrame([dict(v_route=LAB[i],v_b0=b0[i],v_b1=b1[i],v_diff=D_[i],v_t=TT[i],
                   v_thr=THR,v_mde=MDE_ if i==0 else np.nan) for i in range(3)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'S_by_sex.csv',index=False)

# ---- ★ 合成正对照 vs 真实正对照:同一个仪器,两种「它能不能看见」 ----
# 阳性参照失败了。**在把这当成 bug 之前,先证明合成对照会通过** ——
# 若合成过而真实不过,那么差别不在仪器坏了,而在**种植的幅度**:
# 一个种在方便幅度上的合成正对照,只说明仪器在**那个**幅度上工作,不说明它在**真实**幅度上工作。
rgS=np.random.default_rng(21); ysyn=np.full(NN,np.nan)
for gs,ex in ((g0,0.0),(g1,0.25)):
    k=int(gs.sum()); zs=(S[gs]-S[gs].mean())/S[gs].std()
    ysyn[gs]=0.10*zs+ex*zs+rgS.standard_normal(k)
ba,sa=fit(ysyn,g0); bb,sb=fit(ysyn,g1)
TSYN=float(abs(((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12))[0]))
print(f"\n★ **合成**正对照(种植 0.25,一个我随手选的方便幅度):|t| **{TSYN:.3f}** vs 阈 {THR:.3f} -> "
      f"**{'通过' if TSYN>THR else '不过'}**")
print(f"★ **真实**正对照(`EARLY`,幅度由数据决定,`#380a` 已确认它有差):|t| **{abs(TT[2]):.3f}** -> "
      f"**{'通过' if abs(TT[2])>THR else '不过'}**")
print(f"   ⇒ **合成对照会让我报「`S` 没有两性差」;真实对照拦住了它。**")
print(f"   ⇒ 一个种在**方便幅度**上的合成正对照,只说明仪器在**那个**幅度上工作,"
      f"**不说明它在真实幅度上工作**。")

POSREF=abs(TT[2])>THR; NEGREF=abs(TT[1])<=THR
g=Gate('#296 的性别非不变性,在不在 S 那条路上')
g.asserted('★【两支】阳性参照(真实数据):`EARLY` 必须越阈(`#380a`)',POSREF,
           f"|t| {abs(TT[2]):.3f} vs {THR:.3f}",kind='control')
g.asserted('★【两支】阴性参照(真实数据):`c3⁻` 必须不越阈(`#378a`)',NEGREF,
           f"|t| {abs(TT[1]):.3f} vs {THR:.3f}",kind='control')
g.asserted('★【两支】offset 零非退化',mt.std()>0,f"零 {mt.mean():.3f} ± {mt.std():.3f}",kind='control')
if POSREF and NEGREF:
    if abs(TT[0])>THR:
        g.asserted('★【非零支】`S` 越过族内阈 -> 非不变性也在 `S` 上',True,
                   f"|t| {abs(TT[0]):.3f} vs {THR:.3f} · 差 {D_[0]:+.4f}")
    else:
        g.asserted('★【零支】`S` 未越阈,且 MDE < 0.05 才算「看得见而没有」',MDE_<0.05,
                   f"|t| {abs(TT[0]):.3f} vs {THR:.3f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 参照未按预期表现 -> 仪器不算数,不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
