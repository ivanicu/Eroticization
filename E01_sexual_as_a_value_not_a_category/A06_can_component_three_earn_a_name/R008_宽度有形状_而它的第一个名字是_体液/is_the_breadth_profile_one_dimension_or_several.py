import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A57 R274 -- 宽度剖面是一维(幅度)还是几维(类型)

`#228a`:人×块残差宽度剖面复现 +0.4290,零 +0.0043 —— 剖面**可靠**。
但那只是一个信度,**它没有说剖面长什么样**。

WORLDS          ① **幅度**:一个人整体更宽,只是强弱不同 -> 块间残差相关矩阵**一维**
                ② **类型**:有人在 A 类上宽、在 B 类上窄 -> 前 2–3 个特征值都超零
ESTIMAND        块×块**跨半**残差相关矩阵的特征谱。
                ⚠ 用跨半(块 b 的 A 半 × 块 b' 的 B 半)而不是同半 ——
                **同半共享测量噪声,会把噪声算进结构里**;跨半的共同成分只能是人。
KILL            **若 λ1 吃掉正特征值质量的一大半且只有它超零 -> 幅度,领域结构只是强弱;
                若 ≥2 个特征值超最大统计量零阈 -> 存在真正的【宽度类型】,
                那是 A/B/C 之外的第四件事。**
NEGATIVE CTRL   两半各自独立块内跨人置换(`#228c` 的教训:必须在**劈开之后**打乱)。
                零谱的**最大**特征值给全族阈值(最大统计量,不是逐个比)。
POSITIVE CTRL   种入 **2 类人**,各自只在互不相交的一半块上更宽 ->
                谱上必须出现 **2 个**超零特征值。
                ⚠ 强度必须**扫描**,增量对着**实测基线**定价(`#228d`,`#211a` 的重犯)。
IMPOSSIBLE      特征向量的符号与旋转不可识别 -> 只能读**维数**与**块的分组**,
                不能读"第一维是什么" —— 命名要靠块内容,不靠载荷大小。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]; LAB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl)); LAB.append(str(q.col)[:58])
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
print(f"块 {NB} 个;n = {int(ok.sum()):,}")

def halves(seed, plant=None, perm=False):
    """每块两半选项的人层勾选比例。plant: (NN,NB) 的人×块特异宽度。"""
    rg=np.random.default_rng(seed); H=np.full((2,NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        Mm=M if plant is None else np.clip(M+plant[ppl,b][:,None],0,1)
        o=rg.permutation(Mm.shape[1]); k=Mm.shape[1]//2
        ha=Mm[:,o[:k]].mean(1); hb=Mm[:,o[k:2*k]].mean(1)
        if perm: ha=ha[rg.permutation(len(ha))]; hb=hb[rg.permutation(len(hb))]   # ⚠ 劈开之后
        H[0,b,ppl]=ha; H[1,b,ppl]=hb
    return H
def profile(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)          # ⚠ 留一
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
def spectrum(seed, plant=None, perm=False):
    H=halves(seed,plant,perm); Ra,Rb=profile(H[0]),profile(H[1])
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            m=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
            if m.sum()>300: C[i,j]=np.corrcoef(Ra[i][m],Rb[j][m])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2
    w,v=np.linalg.eigh(C); return w[::-1], v[:,::-1], C

w,V,C=spectrum(500)
wn=np.array([spectrum(600+s,perm=True)[0] for s in range(6)])
thr=float(np.quantile(wn[:,0],0.95))                       # 最大统计量全族阈值
above=int((w>thr).sum()); posmass=w[w>0].sum()
print(f"\n观测特征谱前 6:"+' · '.join(f"{x:+.4f}" for x in w[:6]))
print(f"零谱最大特征值 {wn[:,0].mean():+.4f} ± {wn[:,0].std():.4f} -> 全族阈值 **{thr:+.4f}**")
print(f"**超阈值的特征值:{above} 个**;λ1 占正特征值质量 **{100*w[0]/posmass:.1f}%**;"
      f"λ2/λ1 = {w[1]/w[0]:.3f}")

# 正对照:2 类人,各自只在互不相交的一半块上更宽 —— 强度扫描
rg=np.random.default_rng(20260804); hp=rg.permutation(NB); h1,h2=hp[:NB//2],hp[NB//2:]
SW=[]
for gp in (0.05,0.10,0.20,0.40):
    P2=np.zeros((NN,NB)); a=rg.standard_normal(NN)*gp; b2=rg.standard_normal(NN)*gp
    P2[:,h1]=a[:,None]; P2[:,h2]=b2[:,None]
    wp=spectrum(700,plant=P2)[0]; SW.append((gp,int((wp>thr).sum()),float(wp[1]/wp[0])))
print(f"正对照(2 类人各宽一半块)g -> (超阈值个数, λ2/λ1):"
      +' · '.join(f"{a_:.2f}->({b_},{c_:.3f})" for a_,b_,c_ in SW)
      +f"  [基线 ({above},{w[1]/w[0]:.3f})]")

# 前两维的块分组(命名只靠块内容,不靠载荷大小)
for k in range(min(3,above if above>0 else 2)):
    o=np.argsort(-V[:,k])
    print(f"\n第 {k+1} 维 高端 3 块:"+' | '.join(f"{LAB[i][:44]} {V[i,k]:+.2f}" for i in o[:3]))
    print(f"          低端 3 块:"+' | '.join(f"{LAB[i][:44]} {V[i,k]:+.2f}" for i in o[-3:]))

T=pd.DataFrame(dict(idx=np.arange(NB),label=LAB,**{f"pc{k+1}":V[:,k] for k in range(3)}))
T2=pd.DataFrame(dict(v_rank=np.arange(1,NB+1),eig=w,null_max_thr=thr))
check_columns(T,'R274'); check_columns(T2,'R274')
T.to_csv(pathlib.Path(__file__).parent/'results'/'block_loadings.csv',index=False)
T2.to_csv(pathlib.Path(__file__).parent/'results'/'spectrum.csv',index=False)

# ---------- 第二段:跨人半样本的【子空间复现】 ----------
# ⚠ 第一段的统计量不合身,两处:
#   ① 零选错了【类】—— 世界①(幅度)本来就有人层结构,零不该是零。
#      「超零阈的个数」测的是"有没有人层结构",不是"有几维"。
#   ② 正对照往反方向走(种 2 类人 -> 20→14),一个反向移动的统计量不能用来判维数。
# 合身的问法:把【人】劈两半,各自算一次 C,问【前 k 维张成的子空间】复不复现。
#   rank-1 世界:k=1 复现,k>=2 掉到偶然水平(≈ k/NB)。
#   多维世界:k=2,3 仍然复现。
# ⚠ 特征向量在 λ 接近时会自由旋转,所以判【子空间】不判单个向量。
def subspace_rep(seed, plant=None, perm=False, ks=(1,2,3,5)):
    rg=np.random.default_rng(seed); rows=np.flatnonzero(ok); pp=rg.permutation(rows); h=len(pp)//2
    Vs=[]
    for idx in (pp[:h],pp[h:]):
        m=np.zeros(NN,bool); m[idx]=True
        H=halves(seed+1,plant,perm); Ra,Rb=profile(H[0]),profile(H[1])
        C2=np.full((NB,NB),np.nan)
        for i in range(NB):
            for j in range(NB):
                g2=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
                if g2.sum()>150: C2[i,j]=np.corrcoef(Ra[i][g2],Rb[j][g2])[0,1]
        C2=np.where(np.isfinite(C2),C2,0.0); C2=(C2+C2.T)/2
        Vs.append(np.linalg.eigh(C2)[1][:,::-1])
    out={}
    for k in ks:
        P1=Vs[0][:,:k]@Vs[0][:,:k].T
        out[k]=float(np.linalg.norm(P1@Vs[1][:,:k],'fro')**2/k)   # 1.0 = 完全复现,k/NB = 偶然
    return out
SR =[subspace_rep(900+s) for s in range(4)]
SRn=[subspace_rep(950+s,perm=True) for s in range(3)]
KS=(1,2,3,5)
print(f"\n第二段 跨人半样本子空间复现(1.0 = 完全复现,偶然水平 = k/{NB}):")
for k in KS:
    o=np.mean([x[k] for x in SR]); n_=np.mean([x[k] for x in SRn])
    print(f"  k={k}: **{o:.4f}** ± {np.std([x[k] for x in SR]):.4f} · 置换零 {n_:.4f} · 偶然 {k/NB:.4f}")
P2b=np.zeros((NN,NB)); _rg=np.random.default_rng(31)
P2b[:,h1]=(_rg.standard_normal(NN)*0.20)[:,None]; P2b[:,h2]=(_rg.standard_normal(NN)*0.20)[:,None]
SRp=[subspace_rep(980+s,plant=P2b) for s in range(2)]
print(f"  正对照(种 2 类人 g=0.20):"+' · '.join(f"k={k} {np.mean([x[k] for x in SRp]):.4f}" for k in KS))
r1=np.mean([x[1] for x in SR]); r2=np.mean([x[2] for x in SR]); r3=np.mean([x[3] for x in SR])
n2=np.mean([x[2] for x in SRn])

# ---------- 第三段:把特征向量换成一个【旋转无关】的说法 ----------
# ⚠ 特征向量有旋转不可识别性(IMPOSSIBLE 里写过),所以"PC2 是体液维"不能直说。
#    可直说的是块的【分组】:这四块彼此的跨半相关,是否高于它们与其余块的相关。
FAM=[i for i,l in enumerate(LAB) if l.lower().startswith(('for squirt','for saliva','for urine','for precum'))]
oth=[i for i in range(NB) if i not in FAM]
off=~np.eye(NB,dtype=bool)
within=float(np.mean([C[i,j] for i in FAM for j in FAM if i!=j]))
cross =float(np.mean([C[i,j] for i in FAM for j in oth]))
rest  =float(np.mean([C[i,j] for i in oth for j in oth if i!=j]))
allmu =float(C[off].mean())
print(f"\n第三段 旋转无关的分组(体液族 {len(FAM)} 块:squirt/saliva/urine/precum)")
print(f"  族内平均跨半相关 **{within:+.4f}** · 族与其余 {cross:+.4f} · 其余之间 {rest:+.4f} · 全体 {allmu:+.4f}")
Cn=spectrum(660,perm=True)[2]
wn_in=float(np.mean([Cn[i,j] for i in FAM for j in FAM if i!=j]))
print(f"  置换零的族内相关 {wn_in:+.4f};**族内 / 其余之间 = {within/rest:.2f}×**")

g=Gate('宽度剖面:一维还是几维')
g.asserted('第三段(旋转无关):体液族四块的族内相关明显高于其余块之间',
           within>rest*1.5 and within>wn_in+0.05,
           f"族内 {within:+.4f} · 其余之间 {rest:+.4f}({within/rest:.2f}×)· 置换零 {wn_in:+.4f}")
g.asserted('⚠ 第一段的统计量不合身:正对照往反方向走,且零选错了类(世界①本来就有人层结构)',
           False, f"种 2 类人 -> 超阈值 {SW[0][1]}→{SW[-1][1]}(反向);零是「无人层结构」而不是「rank-1」")
g.asserted('正对照(第二段):种 2 类人 -> k=2 子空间必须复现',
           np.mean([x[2] for x in SRp])>0.7, ' · '.join(f"k={k} {np.mean([x[k] for x in SRp]):.4f}" for k in KS))
g.negative_control('k=2 子空间复现的置换零',abs(n2),abs(r2),
                   null_spread=float(np.std([x[2] for x in SRn])),
                   null_kind='两半各自独立块内跨人置换 —— 打掉人层结构;⚠ 这是「无结构」零,不是「rank-1」零')
g.offset_control('★ k=2 vs k=1 子空间复现(判维数的那一格)',r2,r1,
                 float(np.std([x[2] for x in SR])+np.std([x[1] for x in SR])),
                 null_kind='同一条管道在 k=1 上的复现 —— 不是零假设,是「若剖面只有一维,k=2 该掉到哪」')
g.asserted('★ 注册的 kill:≥2 维复现且明显高于偶然 -> 存在真正的【宽度类型】',
           r2>0.7 and r2>3*(2/NB), f"k=1 {r1:.4f} · k=2 {r2:.4f} · k=3 {r3:.4f};偶然 {2/NB:.4f}")
g.asserted('正对照:种 2 类人 -> 超阈值个数必须随强度上升到 >=2(强度已扫描,#228d)',
           max(x[1] for x in SW)>=2 and SW[-1][1]>=SW[0][1], f"{SW} vs 基线 {above}")
g.negative_control('零谱最大特征值',float(wn[:,0].mean()),float(w[0]),
                   null_spread=float(wn[:,0].std()),
                   null_kind='两半各自独立块内跨人置换 —— 打掉人层结构,保留每块的选项边际')
g.asserted('(第一段的注册 kill,保留记录:统计量已判定不合身)',
           above>=2, f"超阈值 {above} 个;λ1 占正质量 {100*w[0]/posmass:.1f}%;λ2/λ1 = {w[1]/w[0]:.3f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T2.to_csv(index=False).encode()).hexdigest()[:12]}")
