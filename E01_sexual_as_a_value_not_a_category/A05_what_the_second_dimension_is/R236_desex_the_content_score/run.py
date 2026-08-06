import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A37 R236 -- 内容维度是不是性别的影子

`#190` 的六道 GENUINE 里,**三道是性别/身份**(`biomale` +0.0769 · `allrollidentity` −0.0423
· 想象自己是女性自慰 −0.0401),而 **`biomale` 是一个人口学变量,不是偏好**。
`#164` 曾给**位置**侧做过这个检验(去性别后保留 102%),**内容侧从没做过。**

ESTIMAND        把 `biomale` 从**每一个块层内容分**里回归掉(`#164` 同款),重算 C_desex,
                重跑那六道的 `c_j`(未残差化的原始内容相关)。
KILL            **若 `animated`/`written` 的 c 掉到全族阈值以下 -> 内容维度很大程度上就是性别,
                `#188` 的"媒介"读法要改成"性别的媒介表达"。**
POSITIVE CTRL   **`biomale` 自己**:去性别之后它的 c 必须塌到 ~0。塌不掉 = 去性别没生效。
NEGATIVE CTRL   纯**个人**内容种植(与性别无关)在去性别后必须存活 —— 否则是过度去除(`#164` 同款)。
NOISE FLOOR     人层 bootstrap 200;全族阈值用最大统计量零重算(去性别后管道变了,阈值不能沿用)。
IMPOSSIBLE      `biomale` 是自报的二元项,不是性别的全部;去掉它不等于去掉性别相关的一切。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
sex=pd.to_numeric(df['biomale'],errors='coerce').values.astype(float)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df)
rb=np.random.default_rng(20260803)
plant_u=rb.standard_normal(NN)          # 纯个人内容种植向量(与性别无关)

def build(desex, plant=0.0):
    con=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN); share=[]
    for _,q in keep.iterrows():
        s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
        ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
        if len(ppl)<1200 or len(opt)<8: continue
        pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
        M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
        if plant:
            sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
            M=M+plant*np.outer(plant_u[ppl],sub)
        Z=M-M.mean(0,keepdims=True)
        w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); sc=Z@v[:,-1]
        if desex:
            # ⚠ #191a:第一版只对 `sex` 回归 —— 那让聚合分在**边际上**与性别正交
            #   (`check_residualized` 因此通过),而下游 `rr_` 判的是**给定勾选数的偏相关**。
            #   于是正对照读到 +0.0769 -> −0.0793:不是没塌,是**在另一个条件下重新出现**。
            #   **残差化的条件集,必须与评估时的条件集相同。** 现在同时回归掉块内勾选数。
            g=sex[ppl]; kk=M.sum(1).astype(float); m=np.isfinite(g)
            if m.sum()>100:
                X=np.c_[np.ones(m.sum()),g[m],kk[m]]
                pred=X@np.linalg.lstsq(X,sc[m],rcond=None)[0]
                r2=1-((sc[m]-pred)**2).sum()/max(((sc[m]-sc[m].mean())**2).sum(),1e-9)
                share.append(r2)
                sc=sc.copy(); sc[m]=sc[m]-pred
        con[ppl]+=sc; KB[ppl]+=M.sum(1); cnt[ppl]+=1
    ok=cnt>=8
    return np.where(ok,con/np.maximum(cnt,1),np.nan), np.where(ok,KB,np.nan), share

C0,KB,_=build(False)
C1,_,share=build(True)
print(f"性别解释块层内容分方差:中位 {100*np.median(share):.2f}% · 最大 {100*max(share):.2f}%")
base=np.isfinite(C0)&np.isfinite(C1)&np.isfinite(KB); bi=np.flatnonzero(base)
check_residualized(C1[bi],np.nan_to_num(sex[bi]),'R236 去性别后的内容分',tol=0.05)
print(f"n = {len(bi):,};corr(C0, C1) = {np.corrcoef(C0[bi],C1[bi])[0,1]:+.4f}")

def rr_(y,x,ii):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    XX=np.c_[np.ones(len(jj)),KB[jj]]
    ry=y[jj]-XX@np.linalg.lstsq(XX,y[jj],rcond=None)[0]
    rx=x[jj]-XX@np.linalg.lstsq(XX,x[jj],rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])

rows=[]; nulls=[]
for c in lik:
    y=df[c].values.astype(float)
    r0=rr_(y,C0,bi); r1=rr_(y,C1,bi)
    sd=float(np.std([rr_(y,C1,rb.choice(bi,len(bi),replace=True)) for _ in range(200)]))
    ps=[]
    for _ in range(40):
        yp=y.copy(); yp[bi]=rb.permutation(y[bi]); v=rr_(yp,C1,bi)
        if np.isfinite(v): ps.append(abs(v))
    if len(ps)>=20: nulls.append(ps)
    rows.append(dict(q=c[:60],c_before=r0,c_after=r1,sd=sd,
                     retained=r1/r0 if abs(r0)>1e-6 else np.nan))
T=pd.DataFrame(rows); check_columns(T,'R236')
L=min(len(x) for x in nulls)
thr=float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nulls]),axis=0),0.95))
T['clears']=T.c_after.abs()>thr
T.to_csv(pathlib.Path(__file__).parent/'results'/'desex.csv',index=False)
print(f"\n去性别后的全族阈值 |c| = {thr:.4f}\n")
SIX=['animated','written','biomale','allrollidentity','submissive','biological \\*female\\*']
print(f"{'去性别前':>10}{'去性别后':>10}{'保留':>8}{'越阈':>6}  题")
for _,r in T.iterrows():
    if any(k.replace('\\','') .lower() in r.q.lower() for k in SIX):
        print(f"{r.c_before:>+10.4f}{r.c_after:>+10.4f}{100*r.retained:>7.0f}%{'★' if r.clears else ' ':>5}  {r.q[:50]}")

# 负对照:纯个人内容种植,去性别后必须存活
Cp0,_,_=build(False,plant=0.6); Cp1,_,_=build(True,plant=0.6)
m=np.isfinite(Cp0)&np.isfinite(Cp1)&np.isfinite(plant_u)
r_p0=float(np.corrcoef(Cp0[m],plant_u[m])[0,1]); r_p1=float(np.corrcoef(Cp1[m],plant_u[m])[0,1])
print(f"\n负对照(纯个人内容种植):去性别前 {r_p0:+.4f} -> 后 {r_p1:+.4f}  保留 {100*r_p1/r_p0:.0f}%")
ani=T[T.q.str.contains('animated')].iloc[0]; wri=T[T.q.str.contains('written')].iloc[0]
bio=T[T.q=='biomale'].iloc[0]
g=Gate('内容维度是不是性别的影子')
g.asserted('正对照:`biomale` 自己去性别后必须塌掉',abs(bio.c_after)<thr,
           f"{bio.c_before:+.4f} -> {bio.c_after:+.4f}(阈值 {thr:.4f})")
g.asserted('负对照:纯个人内容种植去性别后必须存活',r_p1/r_p0>0.7,
           f"保留 {100*r_p1/r_p0:.0f}%")
g.asserted('注册的 kill:animated/written 掉到阈值以下 -> 内容维度就是性别',
           (not ani.clears) and (not wri.clears),
           f"animated {ani.c_after:+.4f}({'越阈' if ani.clears else '未越'}) · "
           f"written {wri.c_after:+.4f}({'越阈' if wri.clears else '未越'})")
g.resolvable('去性别后的 animated',float(ani.c_after),float(ani.sd))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
