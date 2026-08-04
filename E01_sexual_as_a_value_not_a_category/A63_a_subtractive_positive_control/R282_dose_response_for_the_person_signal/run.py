import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A63 R282 -- 减法式正对照:人层信号的剂量-反应曲线

`#236c①`:本项目的正对照几乎全是**加法式**(种入结构,看数变大),
而在**基线贴顶**的量上它结构上无效 —— `#281` 的 k=2 基线 0.9379、上限 1.0,
种入后给出 0.905 / 0.896 / 0.843 / 0.986,**非单调**。

⚠ **先纠正我自己在 `#236` NEXT 里写下的一个推断,它不成立**:
我写了「若曲线在小 λ 就急降 -> 结构集中在少数人」。**这不对** ——
随机替换一个 λ 比例的人,平均而言会同等比例地打掉集中型与分布型的贡献,
两者都近似线性下降。**这个操作分不开那两个世界。**
它能给的是另一样东西,而那正是 `#236c` 真正欠的:**一条真的功效曲线**。

ESTIMAND        对两个已建立的结构量,把随机 λ 比例的人替换成置换噪声,扫 λ:
                (a) 块仪器:人×块残差剖面的两半复现(`#228a` 基线 0.4290,零 +0.0043)
                (b) 起始仪器:类别×类别残差子空间 k=2 复现(`#236a` 基线 0.9379,零 0.4587)
KILL            **两端必须都对上:λ=0 落在各自基线上,λ=1 落在各自置换零上。
                任一端对不上 -> 这个减法式对照不合身,不能用它读中间的曲线。**
交付            **检测极限**:λ 多大时统计量掉到零+2sd 以内 ——
                以「人层信号被稀释到多少还看得见」为单位,而且**两把仪器可比**。
IMPOSSIBLE      λ 混合的是「谁的行数」,不是「信号强度」;
                所以曲线的**形状**只在这一种稀释方式下有意义,不能外推到别的噪声。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
OBS=np.column_stack([np.isfinite(d[c].map(BIN).values.astype(float)) for c in ons]).astype(float)
NC=OBS.shape[1]; okO=OBS.sum(1)>=8
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
okB=cov>=8
print(f"块 {NB}(n={int(okB.sum()):,})· 起始类别 {NC}(n={int(okO.sum()):,})")

def block_rep(seed, lam):
    """人×块残差剖面的两半复现;随机 lam 比例的人被块内置换。"""
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        ha=M[:,o[:k]].mean(1); hb=M[:,o[k:2*k]].mean(1)
        n=len(ppl); hit=rg.random(n)<lam
        if hit.any():
            idx=np.flatnonzero(hit)
            ha[idx]=ha[rg.permutation(idx)]; hb[idx]=hb[rg.permutation(idx)]   # 两半各自独立(#228c)
        A[b,ppl]=ha; B[b,ppl]=hb
    def prof(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof(A),prof(B); g=np.isfinite(Ra)&np.isfinite(Rb)&okB[None,:]
    return float(np.corrcoef(Ra[g],Rb[g])[0,1])
def onset_k2(seed, lam):
    rg=np.random.default_rng(seed); X=OBS.copy()
    hit=np.flatnonzero((rg.random(NN)<lam)&okO)
    for c in range(NC): X[hit,c]=X[rg.permutation(hit),c]
    rows=np.flatnonzero(okO); p=rg.permutation(rows); h=len(p)//2; Vs=[]
    for idx in (p[:h],p[h:]):
        m=np.zeros(NN,bool); m[idx]=True; Z=np.where(m[:,None],X,np.nan)
        tot=np.nansum(Z,1); ct=np.sum(np.isfinite(Z),1); R=np.full_like(Z,np.nan)
        for c in range(NC):
            lo=(tot-np.nan_to_num(Z[:,c]))/np.maximum(ct-np.isfinite(Z[:,c]),1); R[:,c]=Z[:,c]-lo
        R=R-np.nanmean(R,0,keepdims=True)
        C=np.zeros((NC,NC))
        for i in range(NC):
            for j in range(NC):
                if i==j: continue
                g=np.isfinite(R[:,i])&np.isfinite(R[:,j])
                if g.sum()>300: C[i,j]=np.corrcoef(R[g,i],R[g,j])[0,1]
        C=(C+C.T)/2; Vs.append(np.linalg.eigh(C)[1][:,::-1])
    P1=Vs[0][:,:2]@Vs[0][:,:2].T
    return float(np.linalg.norm(P1@Vs[1][:,:2],'fro')**2/2)

LAM=[0.0,0.25,0.5,0.75,0.9,1.0]
rows=[]
print(f"\n{'λ':>6}{'块:剖面复现':>16}{'起始:k=2 子空间':>18}")
for lam in LAM:
    b=[block_rep(1200+s,lam) for s in range(3)]; o=[onset_k2(1300+s,lam) for s in range(3)]
    rows.append(dict(lam=lam,block=float(np.mean(b)),block_sd=float(np.std(b)),
                     onset=float(np.mean(o)),onset_sd=float(np.std(o))))
    print(f"{lam:>6.2f}{np.mean(b):>10.4f}±{np.std(b):<5.4f}{np.mean(o):>12.4f}±{np.std(o):<5.4f}")
T=pd.DataFrame(rows); check_columns(T,'R282')
T.to_csv(pathlib.Path(__file__).parent/'results'/'dose_response.csv',index=False)
b0,b1=T.block.iloc[0],T.block.iloc[-1]; o0,o1=T.onset.iloc[0],T.onset.iloc[-1]
def limit(col,null,nsd):
    for _,r in T.iterrows():
        if r[col]<=null+2*nsd: return float(r.lam)
    return float('nan')
lb=limit('block',b1,float(T.block_sd.iloc[-1])); lo=limit('onset',o1,float(T.onset_sd.iloc[-1]))
print(f"\n检测极限(掉到自己的 λ=1 端 +2sd 以内的最小 λ):块 **{lb}** · 起始 **{lo}**")
print(f"两端:块 λ=0 {b0:.4f}(`#228a` 报 0.4290)· λ=1 {b1:.4f}(零报 +0.0043);"
      f"起始 λ=0 {o0:.4f}(`#236a` 报 0.9379)· λ=1 {o1:.4f}(零报 0.4587)")

g=Gate('减法式正对照:剂量-反应曲线')
g.asserted('★ 两端必须都对上 —— λ=0 落在基线上(块 0.4290 / 起始 0.9379)',
           abs(b0-0.4290)<0.05 and abs(o0-0.9379)<0.08, f"块 {b0:.4f} · 起始 {o0:.4f}")
g.asserted('★ 两端必须都对上 —— λ=1 落在置换零上(块 +0.0043 / 起始 0.4587)',
           abs(b1-0.0043)<0.03 and abs(o1-0.4587)<0.08, f"块 {b1:.4f} · 起始 {o1:.4f}")
g.asserted('曲线单调不增(稀释不该让结构变强)',
           all(T.block.iloc[i]>=T.block.iloc[i+1]-0.02 for i in range(len(T)-1)) and
           all(T.onset.iloc[i]>=T.onset.iloc[i+1]-0.05 for i in range(len(T)-1)),
           f"块 {T.block.round(4).tolist()};起始 {T.onset.round(4).tolist()}")
g.asserted('⚠ 已撤回的推断:随机替换一部分人分不开「集中在少数人」与「分布在全体」',
           True, '两种世界在这个操作下都近似线性下降 —— `#236` NEXT 里那句话不成立,本轮不据此读曲线')
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
