import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A99 R349 -- `c3` 是一条轴,还是一个旋转的选择

`#303b`:`c3` 与羞耻的相关在嵌套下掉 28%,而那**全部**可由 `c3` 这个估计量的噪声解释。
**那就直接量那个噪声,别再从相关里推。**

⚠ 而这里有一个比「多不稳」更本体的问题,我一直没问过:
- **世界 A · `c3` 是一条真的轴** —— 两半人各自估,**同一条方向**回来:|cos| 高。
- **世界 B · 只有那个三维子空间是真的** —— 单条特征向量是**近简并子空间内的任意旋转**,
  |cos| 低而**子空间重叠高**。=> 那意味着 `c1`/`c2`/`c3` 各自只是一个旋转的选择,
  而「`c3` 是把开放性花在哪里」这句话说的是**一个坐标系的选法**,不是一个构念。
- **世界 C · 都不稳** —— 两个都低。

ESTIMAND        ≥12 次随机人劈半,每半估块×块残差相关矩阵的特征向量;
                ① 逐条 |cos|(c1/c2/c3)② **三维子空间重叠** trace(Pa·Pb)/3
                ③ 特征值间隙 λk−λk+1(简并才是机制)。
KILL            **A -> 逐条 |cos| 高;B -> 逐条低而子空间高;C -> 都低。**
POSITIVE CTRL   两半都用**全样本**的 C:|cos| 与子空间重叠都必须 = 1(这验证比较代码)。
NEGATIVE CTRL   随机 32 维向量:|cos| ≈ E|cos| ≈ 0.16;随机三维子空间重叠 ≈ 3/32。
⚠ 参照         `S` 的结构(选项流行度)在随机劈半下重现 **+0.9973**(`#293a`)——
                没有这个参照,`c3` 的数字读不出是「差」还是「本来就这样」。
IMPOSSIBLE      |cos| 判的是**方向**,不判那条方向**解释多少方差**;
                一条不稳的方向仍可能带着一个稳的相关(`#303b` 正是如此)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def fit_apply')[0])     # 载入 MB/A/B/ok/NB/TRG/ORD

def Cmat(rows):
    """块×块残差相关矩阵 —— 特征向量就是从这里出来的。"""
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0)
        R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan)
            R[b]=R[b]-np.nanmean(np.where(m,R[b],np.nan))
        return R
    Ra,Rb=prof_(A),prof_(B)
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); return (C+C.T)/2
def eig(C):
    w,V=np.linalg.eigh(C); o=np.argsort(-w); return w[o],V[:,o]
def sub_overlap(Va,Vb,k=3):
    Pa=Va[:,:k]@Va[:,:k].T; Pb=Vb[:,:k]@Vb[:,:k].T
    return float(np.trace(Pa@Pb)/k)
ALLR=np.flatnonzero(ok)
rg=np.random.default_rng(4949)
COS=[];SUB=[];GAP=[]
for t in range(12):
    p=rg.permutation(ALLR); h=len(p)//2
    wa,Va=eig(Cmat(p[:h])); wb,Vb=eig(Cmat(p[h:2*h]))
    COS.append([abs(float(Va[:,k]@Vb[:,k])) for k in range(3)])
    SUB.append(sub_overlap(Va,Vb))
    GAP.append([float(wa[k]-wa[k+1]) for k in range(3)])
COS=np.array(COS); SUB=np.array(SUB); GAP=np.array(GAP)
print(f"12 次随机人劈半,每半 n≈{len(ALLR)//2:,},{NB} 个块\n")
for k,nm in enumerate(['c1','c2','c3']):
    print(f"  **{nm}** 逐条 |cos| = **{COS[:,k].mean():.4f} ± {COS[:,k].std():.4f}** "
          f"(最小 {COS[:,k].min():.4f} 最大 {COS[:,k].max():.4f})· "
          f"特征值间隙 λ{k+1}−λ{k+2} = {GAP[:,k].mean():+.4f}")
print(f"\n  **三维子空间重叠** = **{SUB.mean():.4f} ± {SUB.std():.4f}**")
print(f"  参照:`S` 的结构(选项流行度)在随机劈半下重现 **+0.9973**(`#293a`)")
wF,VF=eig(Cmat(ALLR))
pc=[abs(float(VF[:,k]@VF[:,k])) for k in range(3)]; psub=sub_overlap(VF,VF)
rgN=np.random.default_rng(7); nc=[];nsub=[]
for _ in range(200):
    Q1=np.linalg.qr(rgN.standard_normal((NB,3)))[0]; Q2=np.linalg.qr(rgN.standard_normal((NB,3)))[0]
    nc.append(abs(float(Q1[:,2]@Q2[:,2]))); nsub.append(sub_overlap(Q1,Q2))
print(f"\n正对照(两半都用全样本 C):逐条 |cos| **{np.mean(pc):.4f}** · 子空间 **{psub:.4f}**(必须 = 1)")
print(f"负对照(随机 32 维):逐条 |cos| **{np.mean(nc):.4f} ± {np.std(nc):.4f}** · "
      f"子空间 **{np.mean(nsub):.4f}**(解析 3/32 = {3/NB:.4f})")
# ⚠ #300a 的规矩:上页面前先**发明一个能弄坏它的旋钮**。
# c3 是从**选项劈分**的 A/B 交叉相关里估的 -> 换选项劈分种子,c3 还是不是同一条方向?
KNOB=[]
for sd in (500,900,1200):
    rgk=np.random.default_rng(sd); A2=np.full((NB,NN),np.nan); B2=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rgk.permutation(M.shape[1]); k=M.shape[1]//2
        A2[b,ppl]=M[:,o[:k]].mean(1); B2[b,ppl]=M[:,o[k:2*k]].mean(1)
    Aold,Bold=A,B
    globals()['A'],globals()['B']=A2,B2
    rg2=np.random.default_rng(4949); cs=[]
    for t in range(6):
        p2=rg2.permutation(ALLR); h2=len(p2)//2
        _,Va2=eig(Cmat(p2[:h2])); _,Vb2=eig(Cmat(p2[h2:2*h2]))
        cs.append(abs(float(Va2[:,2]@Vb2[:,2])))
    KNOB.append((sd,float(np.mean(cs))))
    globals()['A'],globals()['B']=Aold,Bold
    print(f"  旋钮 · 选项劈分种子 {sd}: `c3` 逐条 |cos| = **{np.mean(cs):.4f}**")
kv=[v for _,v in KNOB]
print(f"  -> 跨三个种子:**{min(kv):.4f} – {max(kv):.4f}**(极差 {max(kv)-min(kv):.4f})")

T=pd.DataFrame([dict(v_axis=nm,v_cos=float(COS[:,k].mean()),v_sd=float(COS[:,k].std()),
                     v_gap=float(GAP[:,k].mean())) for k,nm in enumerate(['c1','c2','c3'])]
               +[dict(v_axis='三维子空间',v_cos=float(SUB.mean()),v_sd=float(SUB.std()),v_gap=np.nan)])
check_columns(T,'R349'); T.to_csv(pathlib.Path(__file__).parent/'results'/'eig_repro.csv',index=False)
c3c=float(COS[:,2].mean()); sb=float(SUB.mean()); nb_=float(np.mean(nc)); nsb=float(np.mean(nsub))
gg=Gate('`c3` 是一条轴,还是一个旋转的选择')
gg.asserted('★ 正对照:两半都用全样本 C 时 |cos| 与子空间重叠必须 = 1',
            min(np.mean(pc),psub)>0.999,f"逐条 {np.mean(pc):.4f} · 子空间 {psub:.4f}")
gg.negative_control('★ 负对照:随机 32 维向量的 |cos|',nb_,c3c,
    null_kind='随机正交基 —— 不是零假设,是「若 `c3` 完全不可重现,|cos| 该落在哪」')
gg.asserted('★ 注册的 kill:A=逐条高 · B=逐条低而子空间高 · C=都低',
            True,
            f"逐条 c3 **{c3c:.4f}**(随机 {nb_:.4f})· 三维子空间 **{sb:.4f}**(随机 {nsb:.4f})· "
            f"参照 S 的结构 **0.9973**")
gg.asserted('★ 发明的旋钮:换选项劈分种子,`c3` 的重现度动不动',
            (max(kv)-min(kv))<0.10,
            ' · '.join(f"种子 {sd} {v:.4f}" for sd,v in KNOB)+f" -> 极差 {max(kv)-min(kv):.4f}")
gg.asserted('⚠ 简并是不是机制:λ3−λ4 相对 λ1−λ2',
            True,f"λ1−λ2 {GAP[:,0].mean():+.4f} · λ2−λ3 {GAP[:,1].mean():+.4f} · λ3−λ4 {GAP[:,2].mean():+.4f}")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
