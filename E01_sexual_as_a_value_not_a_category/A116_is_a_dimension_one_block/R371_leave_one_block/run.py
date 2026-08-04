import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A116 R371 -- `S` 与 `c3` 会不会其实是某一块

`#325d`:合成分不等于共同因子 —— **「合成分与 X 的相关」可以几乎全部来自一个指标**。
`form` 只有 2 个指标所以一眼看穿;**`S`/`c3` 是 32 块的合成,看不穿,而它们是页面的主干。**

ESTIMAND        留一块重算(32 次):`corr(S_{-b}, 羞耻)` 与 `corr(c3_{-b}, 羞耻)`,
                报「去掉这一块后相关的变化」的分布,以及**最大的一块占多少**。
KILL            **若最大的一块贡献 > 30% -> 这个维度的名字该换成那一块的名字;
                若分布平 -> 它确实是 32 块共有的东西,页面上的名字站得住。**
POSITIVE CTRL   合成一个**只由一块驱动**的量 -> 留一必须把它挑出来(那一块的变化必须最大且远超其余)。
NEGATIVE CTRL   合成一个**均摊**的量 -> 留一的分布必须平。
⚠ 符号         `c3` 每次重估都要**对齐到全量版本**(guard 20 / `R210:73` 的老坑)。
IMPOSSIBLE      留一测的是**边际**贡献;若两块高度共线,各自留一都不动,而两块一起去掉会动 ——
                本轮测不到那种成对的集中。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def fit_apply')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
m=ok.copy()
def S_drop(drop):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for b,(M,ppl) in enumerate(MB):
        if b==drop: continue
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
        v=np.where(n>0,(M@rr)/np.maximum(n,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(cv>=(7 if drop>=0 else 8),ps/np.maximum(cv,1),np.nan)
def prof_(X,keep):
    F=np.isfinite(X); F=F&np.isin(np.arange(NB),keep)[:,None]
    Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in keep:
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
    return R
def c3_drop(drop,ref=None):
    keep=[b for b in range(NB) if b!=drop]
    Ra,Rb=prof_(A,keep),prof_(B,keep)
    C=np.zeros((NB,NB))
    for i in keep:
        for j in keep:
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
    if ref is not None and float(V[:,2]@ref)<0: V[:,2]=-V[:,2]   # ⚠ 符号对齐(guard 20)
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(V[:,2][:,None]*Zm).sum(0); den=(Fm*np.abs(V[:,2])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan),V[:,2]
def cor(u):
    k=np.isfinite(u)&np.isfinite(sh)&m
    return float(np.corrcoef(u[k],sh[k])[0,1]) if k.sum()>200 else np.nan
S0=S_drop(-1); r_S0=cor(S0)
c30,V0=c3_drop(-1); r_c30=cor(c30)
print(f"全量:`S ↔ 羞耻` **{r_S0:+.4f}** · `c3 ↔ 羞耻` **{r_c30:+.4f}**(块 {NB})")
dS=np.array([cor(S_drop(b))-r_S0 for b in range(NB)])
dC=np.array([cor(c3_drop(b,ref=V0)[0])-r_c30 for b in range(NB)])
def summ(nm,d0,r0):
    a=np.abs(d0); sh_=a/max(np.sum(a),1e-12)
    i=int(np.argmax(a))
    print(f"\n**{nm}**(全量 {r0:+.4f}):留一变化 |Δ| 中位 **{np.median(a):.4f}** · "
          f"最大 **{a[i]:.4f}**(块 #{i})· **最大块占总变化的 {100*sh_[i]:.1f}%**")
    top=np.argsort(-a)[:4]
    print(f"   前四块:" + ' · '.join(f"#{t}(Δ {d0[t]:+.4f})" for t in top))
    return float(sh_[i]),float(a[i]/max(abs(r0),1e-12))
shS,relS=summ('S',dS,r_S0); shC,relC=summ('c3',dC,r_c30)
rg=np.random.default_rng(21)
def synth(one_block):
    v=np.full(NN,np.nan)
    for b,(M,ppl) in enumerate(MB):
        w=(1.0 if (one_block is None or b==one_block) else 0.0)
        if w==0: continue
        z=M.sum(1).astype(float); z=(z-z.mean())/max(z.std(),1e-12)
        cur=np.where(np.isnan(v[ppl]),0.0,v[ppl]); v[ppl]=cur+w*z
    y=np.full(NN,np.nan); g=np.isfinite(v)&m
    y[g]=v[g]+rg.standard_normal(int(g.sum()))
    return y,g
def loo_on(y,g):
    def c_(drop):
        vv=np.full(NN,np.nan)
        for b,(M,ppl) in enumerate(MB):
            if b==drop: continue
            z=M.sum(1).astype(float); z=(z-z.mean())/max(z.std(),1e-12)
            cur=np.where(np.isnan(vv[ppl]),0.0,vv[ppl]); vv[ppl]=cur+z
        k=np.isfinite(vv)&np.isfinite(y)&m
        return float(np.corrcoef(vv[k],y[k])[0,1])
    base=c_(-1); dd=np.array([c_(b)-base for b in range(NB)]); a=np.abs(dd)
    return float(np.max(a)/max(np.sum(a),1e-12)),int(np.argmax(a))
y1,g1=synth(5); p_share,p_arg=loo_on(y1,g1)
y2,g2=synth(None); n_share,n_arg=loo_on(y2,g2)
print(f"\n正对照(只由块 #5 驱动):最大块占 **{100*p_share:.1f}%**,是块 **#{p_arg}**")
print(f"负对照(均摊):最大块占 **{100*n_share:.1f}%**(32 块均摊的期望 {100/NB:.1f}%)")
T=pd.DataFrame([dict(v_dim='S',full=r_S0,max_share=shS,max_rel=relS),
                dict(v_dim='c3',full=r_c30,max_share=shC,max_rel=relC)])
check_columns(T,'R371'); T.to_csv(pathlib.Path(__file__).parent/'results'/'loo.csv',index=False)
gg=Gate('`S` 与 `c3` 会不会其实是某一块')
gg.asserted('★ 正对照:只由一块驱动的量,留一必须把它挑出来',p_share>0.30 and p_arg==5,
            f"最大块占 {100*p_share:.1f}%,是 #{p_arg}(应 #5)")
gg.asserted('★ 负对照:均摊的量,留一分布必须平',n_share<0.15,
            f"最大块占 {100*n_share:.1f}%(均摊期望 {100/NB:.1f}%)")
gg.asserted('★ 注册的 kill:最大的一块贡献是否 > 30%',max(shS,shC)>0.30,
            f"S **{100*shS:.1f}%** · c3 **{100*shC:.1f}%** —— "
            f"{'某一块主导' if max(shS,shC)>0.30 else '**分布平,32 块共有的东西,页面上的名字站得住**'}")
gg.asserted('⚠ 边界:留一测的是边际贡献',True,
            '两块高度共线时各自留一都不动,而一起去掉会动 —— 本轮测不到成对的集中')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
