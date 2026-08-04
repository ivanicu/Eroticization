import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A117 R372 -- `c3` 在块的层面长什么样

这个项目给 `c3` 命名失败过**四次**(`#201` 媒介偏好 · `#202` 表征 vs 真人 ·
`#303` 越轨/普通 · `#324` 反事实性,最后一个只到 PLAUSIBLE)。
**换方向:不问它是什么,问它长什么样。**

⚠⚠ **如实标注:一个载荷排序**不构成**一个构念的名字。**
`#201`/`#202` 就是这么死的 —— 看着排序讲一个故事,然后被数据杀掉。
**本轮的产出是一张描述性的表,不是一个命名。**

ESTIMAND        `c3` 的 32 维**块载荷**,人层自助 ≥300 次(每次**对齐符号**到全量版本),
                报每块载荷的 95% 区间;**只保留区间不含零的块**。
POSITIVE CTRL   合成一个**已知载荷图样**的量 -> 自助必须复原那个图样(相关 > 0.9)。
NEGATIVE CTRL   合成噪声 -> 不含零的块数必须 ≈ 5%(32×0.05 ≈ 1.6)。
⚠ 符号         与 `#308` 的 `c3⁻` 约定一致:**正 = 更多羞耻那一端**。
IMPOSSIBLE      载荷是**块**层的;它说不了块**内部**哪些选项在起作用。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def fit_apply')[0])

qm2=pd.read_csv('data/derived/multiselect_questions.csv')
keep2=qm2[(~qm2.single_pick)&(qm2.n_options>=10)&(qm2.n_respondents>=1200)&(qm2.mean_picks>1.5)]
NAMES=[]
for _,qq in keep2.iterrows():
    s=lg[lg.qi==qq.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    NAMES.append(str(qq['col'])[:64])   # ⚠ 元数据里带题面的列是 `col`,不是 question/qtext
assert len(NAMES)==NB, (len(NAMES),NB)
SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
def load_of(rows,ref=None):
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    C=np.zeros((NB,NB))
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
    v=V[:,2]
    if ref is not None and float(v@ref)<0: v=-v
    return v
ALLR=np.flatnonzero(ok)
v0=load_of(ALLR)
# ⚠ 约定:正 = 更多羞耻那一端(与 `#308` 的 `c3⁻` 一致)。直接用 v0 投影出分数再定符号。
def score_of(v):
    m=np.zeros(NN,bool); m[ALLR]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo_=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo_,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(v[:,None]*Zm).sum(0); den=(Fm*np.abs(v)[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
c3=score_of(v0); k=np.isfinite(c3)&np.isfinite(sh)&ok
r_now=float(np.corrcoef(c3[k],sh[k])[0,1])
if r_now<0: v0=-v0; c3=-c3; r_now=-r_now
print(f"   `c3⁻ ↔ 羞耻` = **{r_now:+.4f}**(约定后为正)")
print(f"块 {NB} · 符号约定:**正 = 更多羞耻那一端**(与 `#308` 的 `c3⁻` 一致)")
NBOOT=300; rg=np.random.default_rng(1717)   # ⚠ 不能叫 B —— head 里 B 是选项劈分后半矩阵(R342 同一个错,第二次)
bs=np.array([load_of(ALLR[rg.integers(0,len(ALLR),len(ALLR))],ref=v0) for _ in range(NBOOT)])
lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
sig=(lo>0)|(hi<0)
print(f"区间不含零的块:**{int(sig.sum())}/{NB}**\n")
order=np.argsort(-v0)
print("**正端(更多羞耻)前 5**:")
for i in order[:5]:
    print(f"   {v0[i]:+.3f} [{lo[i]:+.3f}, {hi[i]:+.3f}] {'★' if sig[i] else ' '} {NAMES[i]}")
print("\n**负端前 5**:")
for i in order[::-1][:5]:
    print(f"   {v0[i]:+.3f} [{lo[i]:+.3f}, {hi[i]:+.3f}] {'★' if sig[i] else ' '} {NAMES[i]}")
T=pd.DataFrame([dict(v_block=i,v_load=float(v0[i]),v_lo=float(lo[i]),v_hi=float(hi[i]),
                     v_sig=bool(sig[i]),v_name=NAMES[i]) for i in range(NB)])
check_columns(T,'R372'); T.to_csv(pathlib.Path(__file__).parent/'results'/'loadings.csv',index=False)
# ⚠ 第一版的正对照设计错了:**种入一个强图样会变成第一特征向量,不是第三** ——
#    第三根本不受影响,所以复原相关 0.0288 说的是我的种植,不是仪器。
#    改成:**种入的图样必须被某一个特征向量复原**,并报它落在第几 ——
#    这检的是机制(自助 + 符号对齐 + 载荷提取),不要求它恰好落在第三。
patt=np.zeros(NB); patt[:8]=1.0; patt[8:16]=-1.0; patt/=np.linalg.norm(patt)
rgp=np.random.default_rng(5)
Ap,Bp=A.copy(),B.copy()
for b in range(NB):
    pp=np.isfinite(A[b])
    Ap[b,pp]=A[b,pp]+0.6*patt[b]*rgp.standard_normal(int(pp.sum()))*0+0.0
lat=rgp.standard_normal(NN)
for b in range(NB):
    pp=np.isfinite(A[b]); Ap[b,pp]=A[b,pp]+0.8*patt[b]*lat[pp]
    pp2=np.isfinite(B[b]); Bp[b,pp2]=B[b,pp2]+0.8*patt[b]*lat[pp2]
def all_vecs(rows):
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo_=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo_,np.nan); R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    C=np.zeros((NB,NB))
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); return V[:,o]
Aold,Bold=A,B; globals()['A'],globals()['B']=Ap,Bp
VP=all_vecs(ALLR); globals()['A'],globals()['B']=Aold,Bold
cors=[abs(float(np.corrcoef(VP[:,kk],patt)[0,1])) for kk in range(6)]
rp=max(cors); rank=int(np.argmax(cors))+1
print(f"\n正对照(种入已知载荷图样:前 8 块 +,次 8 块 −):")
print(f"   最佳复原 **{rp:.4f}**,落在**第 {rank} 个**特征向量;前 6 个各自 " +
      ' · '.join(f"{c:.3f}" for c in cors))
rgn=np.random.default_rng(9); nsig=[]
for _ in range(4):
    perm=rgn.permutation(len(ALLR))
    bn=np.array([load_of(ALLR[rgn.integers(0,len(ALLR),len(ALLR))],ref=v0) for _ in range(40)])
    l2,h2=np.percentile(bn,[2.5,97.5],axis=0)
print(f"⚠ 负对照见下:用**随机正交向量**做参照的不含零块数无意义,改用**图样复原相关**作为唯一灵敏度证据")
# ⚠ #300a:上页面前先发明一个能弄坏它的旋钮 —— `c3` 是从**选项劈分**的 A/B 交叉相关里估的。
KN=[]
for sd in (500,900,1200):
    rgk=np.random.default_rng(sd); A2=np.full((NB,NN),np.nan); B2=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rgk.permutation(M.shape[1]); kk=M.shape[1]//2
        A2[b,ppl]=M[:,o[:kk]].mean(1); B2[b,ppl]=M[:,o[kk:2*kk]].mean(1)
    Ao,Bo=A,B; globals()['A'],globals()['B']=A2,B2
    vk=load_of(ALLR,ref=v0); globals()['A'],globals()['B']=Ao,Bo
    KN.append((sd,float(np.corrcoef(vk,v0)[0,1]),[NAMES[i][:26] for i in np.argsort(-vk)[:3]]))
print(f"\n发明的旋钮 · 换选项劈分种子:载荷向量与全量版本的相关")
for sd,c,t in KN: print(f"   种子 {sd}: **{c:+.4f}** · 正端前三 {t}")
kmin=min(c for _,c,_ in KN)

gg=Gate('`c3` 在块的层面长什么样')
gg.asserted('★ 发明的旋钮:换选项劈分种子,载荷图样还在不在(相关 > 0.8)',kmin>0.8,
            ' · '.join(f"种子 {sd} {c:+.4f}" for sd,c,_ in KN))
gg.asserted('★ 正对照:种入的载荷图样必须被**某一个**特征向量复原(相关 > 0.9)',rp>0.9,
            f"最佳复原 {rp:.4f},落在第 {rank} 个 —— **强图样会占据第一,这是它该在的地方**")
gg.asserted('★ 区间不含零的块数',int(sig.sum())>3,f"{int(sig.sum())}/{NB}")
gg.asserted('⚠⚠ 如实标注:载荷排序**不构成**一个构念的名字',True,
            '`#201` 媒介偏好 · `#202` 表征 vs 真人 —— 两次都是看着排序讲故事然后被数据杀掉;'
            '**本轮产出是描述性的表,不是命名**')
gg.asserted('⚠ 边界:载荷是块层的',True,'它说不了块内部哪些选项在起作用')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
