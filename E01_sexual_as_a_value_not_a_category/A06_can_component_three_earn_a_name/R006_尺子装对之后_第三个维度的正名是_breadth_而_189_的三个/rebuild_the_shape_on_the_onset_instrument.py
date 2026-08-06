import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A62 R281 -- 把「宽度有形状」放到第二把仪器上重建一次

`#235c`:`#229`–`#232` 的结构性数字只是**通过块这把仪器**测到的性质,
而它们**也从没被试过**换一把仪器。这一轮试。

⚠ **一个结构性差异,必须先说**:块仪器的力量来自**块内选项劈半**(同一个人×块格子的两次独立测量)。
起始仪器每个类别只有**一个**二值题,**没有可劈的子项** —— 所以那条路走不通。
可走的是 `#229` 的**第二段**:**劈人**,判类别×类别残差相关矩阵的**前 k 维子空间复不复现**。
⚠ 同半测量噪声只落在**对角线**上(不同类别的噪声独立),所以对角置零后再取谱。

ESTIMAND        起始仪器(31 个类别)上人×类别残差剖面的类别×类别相关矩阵,
                跨**人**半样本的前 k 维子空间复现(`#229` 第二段同款)。
KILL            **若 k=2 复现明显高于置换零与偶然水平(k/31)-> 「宽度有形状」跨仪器成立,
                `#235c` 的收窄放宽为「维数成立,具体的族与成分身份仍是块相关的」;
                若塌到零 -> `#229` 整条按「块格式的性质」重写(而 `#235b` 的羞耻主张不受影响)。**
⚠ 先量功效再判零     起始仪器的分数信度只有 **0.3513**(块是 0.8463)。
                所以**正对照必须扫描强度并报出这个设计能看见的最小复现**;
                **若最小可见复现高于块仪器实测的 0.88,则这个设计结构上无法回答**,
                如实登记为 `P14` 的「本站点结构上做不到」,而不是报一个零。
NEGATIVE CTRL   类别内跨人置换(打掉人层结构,保留每个类别的边际)。
IMPOSSIBLE      两把仪器的**内容**不同,所以即使维数都 ≥2,也**不能**说是同一批维度。
                能判的只有「形状这件事是否也存在于另一把仪器上」。
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
NC=OBS.shape[1]; ok=OBS.sum(1)>=8
print(f"起始仪器:{NC} 个类别;n = {int(ok.sum()):,};类别流行度 "
      f"{OBS.mean(0).min():.3f}–{OBS.mean(0).max():.3f}")

def resid(X, rows):
    m=np.zeros(NN,bool); m[rows]=True
    Z=X.copy(); Z[~m]=np.nan
    lo=np.full_like(Z,np.nan)
    tot=np.nansum(Z,1); ct=np.sum(np.isfinite(Z),1)
    for c in range(NC):
        lo[:,c]=(tot-np.nan_to_num(Z[:,c]))/np.maximum(ct-np.isfinite(Z[:,c]),1)   # ⚠ 留一
    R=Z-lo
    return R-np.nanmean(R,0,keepdims=True)
def cmat(X, rows):
    R=resid(X,rows); C=np.full((NC,NC),np.nan)
    for i in range(NC):
        for j in range(NC):
            g=np.isfinite(R[:,i])&np.isfinite(R[:,j])
            if g.sum()>300: C[i,j]=np.corrcoef(R[g,i],R[g,j])[0,1]
    C=np.where(np.isfinite(C),C,0.0); np.fill_diagonal(C,0.0)     # ⚠ 同半噪声只在对角线
    return (C+C.T)/2
def subspace(seed, plant=None, perm=False, ks=(1,2,3,5)):
    rg=np.random.default_rng(seed); rows=np.flatnonzero(ok); p=rg.permutation(rows); h=len(p)//2
    X=OBS.copy()
    if plant is not None: X=np.clip(X+plant,0,1)
    if perm:
        X=X.copy()
        for c in range(NC): X[:,c]=X[rg.permutation(NN),c]
    Vs=[np.linalg.eigh(cmat(X,idx))[1][:,::-1] for idx in (p[:h],p[h:])]
    return {k: float(np.linalg.norm((Vs[0][:,:k]@Vs[0][:,:k].T)@Vs[1][:,:k],'fro')**2/k) for k in ks}
KS=(1,2,3,5)
SR =[subspace(400+s) for s in range(4)]
SRn=[subspace(450+s,perm=True) for s in range(3)]
print(f"\n跨人半样本子空间复现(1.0 = 完全复现;偶然 = k/{NC}):")
for k in KS:
    o=np.mean([x[k] for x in SR]); s=np.std([x[k] for x in SR]); n_=np.mean([x[k] for x in SRn])
    sn=np.std([x[k] for x in SRn])
    print(f"  k={k}: **{o:.4f}** ± {s:.4f} · 置换零 {n_:.4f} ± {sn:.4f} · 偶然 {k/NC:.4f}"
          f" · **(观测−零)/零的展布 = {(o-n_)/max(sn,1e-9):.1f}×**")

# ⚠ 先量功效:种入 2 类人各自只在一半类别上更宽,扫描强度
rg=np.random.default_rng(20260804); hp=rg.permutation(NC); h1,h2=hp[:NC//2],hp[NC//2:]
PW=[]
for gp in (0.02,0.05,0.10,0.20):
    P=np.zeros((NN,NC)); a=rg.standard_normal(NN)*gp; b=rg.standard_normal(NN)*gp
    P[:,h1]=a[:,None]; P[:,h2]=b[:,None]
    PW.append((gp,subspace(480,plant=P)[2]))
base2=np.mean([x[2] for x in SR]); nul2=np.mean([x[2] for x in SRn]); sd2=np.std([x[2] for x in SR])
print(f"功效扫描(种 2 类人)g -> k=2 复现:"+' · '.join(f"{a_:.2f}->{b_:.4f}" for a_,b_ in PW)
      +f"  [基线 {base2:.4f} · 零 {nul2:.4f}]")
mdl=[b_ for a_,b_ in PW if b_>nul2+2*np.std([x[2] for x in SRn])]
print(f"  这个设计能看见的最小 k=2 复现 ≈ **{min(mdl) if mdl else float('nan'):.4f}**"
      f"(块仪器上实测 0.8805)")

T=pd.DataFrame([dict(k=k,rep=float(np.mean([x[k] for x in SR])),sd=float(np.std([x[k] for x in SR])),
                     null=float(np.mean([x[k] for x in SRn])),chance=k/NC) for k in KS])
check_columns(T,'R281'); T.to_csv(pathlib.Path(__file__).parent/'results'/'onset_subspace.csv',index=False)

g=Gate('把「宽度有形状」放到第二把仪器上')
g.asserted('⚠⚠ 这一条结构上不可能失败,保留记录:基线本身就在零之上,所以 mdl 永不为空',
           bool(mdl), f"[空条款] 扫描 {[(a_,round(b_,4)) for a_,b_ in PW]};最小可见 "
                      f"{min(mdl) if mdl else float('nan'):.4f};块仪器实测 0.8805")
g.negative_control('k=2 子空间复现的置换零',abs(nul2),abs(base2),
                   null_spread=float(np.std([x[2] for x in SRn])),
                   null_kind='类别内跨人置换 —— 打掉人层结构,保留每个类别的边际')
g.offset_control('k=2 vs 偶然水平',base2,2/NC,sd2,
                 null_kind='k/31 —— 不是零假设,是「若前 2 维只是随机方向,子空间该重叠多少」')
g.asserted('★ 注册的 kill:k=2 复现明显高于零与偶然 -> 「宽度有形状」跨仪器成立',
           base2>nul2+2*np.std([x[2] for x in SRn]) and base2>3*(2/NC),
           f"k=2 {base2:.4f} vs 零 {nul2:.4f} · 偶然 {2/NC:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
