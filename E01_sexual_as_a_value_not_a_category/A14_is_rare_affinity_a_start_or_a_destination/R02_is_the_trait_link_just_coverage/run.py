import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A14 R02 -- corr(rho,S) 掉的那一半,是心理学还是仪器?

#128c:「越是最终口味罕见的人,罕见兴趣提前得越多」= -0.0502,但**扣掉年龄与类别数后只剩
-0.0244**。掉的这一半由两个**仪器量**携带,而 A05 的覆盖度定律(#5,corr = +0.815)说:
本 release 上任何未按块数匹配的分组比较,有一大半在测问卷覆盖度而不是性癖差异。

ESTIMAND        corr(rho_i, S_i),在按**答题块数**做卡钳 1:1 匹配之后。
IDENTIFICATION  匹配把覆盖度这条通路**在设计上**关掉,而不是靠回归假设线性。
SCOPE           报告 >=8 个类别起始年龄且 S 可算的人。
WORLDS          psych  匹配后仍在 -> 特质链接是心理学的,128c 可留在 README
                instr  匹配后归零 -> 「越罕见的人越早」只是「答题多的人报得细」,128c 撤回
                partial 介于之间 -> 报匹配后的数,并把范围写清
KILL            条件式:匹配必须真的把块数差压到 <0.1 sd(检验它,不假设),
                且置换零必须为零,才读阈值。
POSITIVE CTRL   种植一个已知的 corr(rho,S):必须在匹配后仍被回收。
NEGATIVE CTRL   题内跨人置换,与 R01 同一个零。
NOISE FLOOR     200 次按人自助 + 5 个匹配种子。
MULTIPLICITY    3 个协变量集 x 5 seeds,整格发表。
IMPOSSIBLE      覆盖度与特质若真的是同一件事,匹配会把信号一起拿掉,而这一轮分不出
                「过度控制」与「本来就没有」。所以只在**匹配后仍在**的方向下结论。
"""
import pandas as pd, numpy as np, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage

# 跨轮依赖显式声明(P16):复用 R01 的加载器与 S 的构造,到 '# ---- beta_i' 为止。
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R01_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('# ---- beta_i')[0])

def demean_conv(Vm,tol=1e-10,cap=500):
    D=np.where(obs,Vm,np.nan)
    for k in range(cap):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-a
        b=np.nanmean(D,axis=1,keepdims=True); D=D-b
        if max(np.nanmax(np.abs(a)),np.nanmax(np.abs(b)))<tol: break
    return D
def rho_of(D):
    rl=np.full(len(D),np.nan)
    for i in np.flatnonzero(KEEP):
        m_=obs[i]; x=rar[m_]-rar[m_].mean(); y=D[i,m_]
        sy=np.sqrt(np.nansum((y-np.nanmean(y))**2)); v=np.sqrt((x*x).sum())
        if sy>1e-9 and v>1e-9: rl[i]=np.nansum(y*x)/(v*sy)
    return rl
def perm_null(Vm,rng):
    Wm=Vm.copy()
    for j in range(Vm.shape[1]):
        idx=np.flatnonzero(obs[:,j]); Wm[idx,j]=Vm[rng.permutation(idx),j]
    return Wm

rho=rho_of(demean_conv(V))
M=np.isfinite(rho)&KEEP&np.isfinite(S)&np.isfinite(age)
ii=np.flatnonzero(M)
NC=NCAT.astype(float)
print(f"可用 {len(ii):,} 人",flush=True)
cS=np.corrcoef(NC[ii],S[ii])[0,1]
print(f"corr(类别数, S) = {cS:+.4f}   corr(类别数, rho) = {np.corrcoef(NC[ii],rho[ii])[0,1]:+.4f}",flush=True)
pk=np.nan_to_num(PK,nan=np.nanmean(PK))
print(f"corr(类别数, 勾选数) = {np.corrcoef(NC[ii],pk[ii])[0,1]:+.4f}   "
      f"corr(勾选数, S) = {np.corrcoef(pk[ii],S[ii])[0,1]:+.4f}(S 已对勾选数残差化)",flush=True)
print(f"S 与「报告了几个类别」的相关 {cS:+.3f}"
      f"{'  ⚠ 覆盖度污染,必须匹配后才可读' if abs(cS)>0.2 else '  (与覆盖度基本正交)'}",flush=True)

def matched_corr(seed,cov):
    """按 cov 做卡钳 1:1 匹配:S 高半 vs S 低半,每个高的人配一个 cov 最近的低的人。"""
    rg=np.random.default_rng(seed)
    med=np.median(S[ii]); hi=ii[S[ii]>med]; lo=ii[S[ii]<=med]
    C=np.c_[tuple((cov[:,k]-cov[ii,k].mean())/cov[ii,k].std() for k in range(cov.shape[1]))]
    used=np.zeros(len(S),bool); pairs=[]
    order=rg.permutation(len(hi))
    lo_c=C[lo]
    for a in hi[order]:
        d=np.abs(lo_c-C[a]).max(1)
        d[used[lo]]=np.inf
        j=int(np.argmin(d))
        if d[j]<0.25: used[lo[j]]=True; pairs.append((a,lo[j]))
    if len(pairs)<200: return np.nan,np.nan,0
    P=np.array(pairs); sel=np.r_[P[:,0],P[:,1]]
    bal=abs(cov[P[:,0],0].mean()-cov[P[:,1],0].mean())/cov[ii,0].std()
    return float(np.corrcoef(rho[sel],S[sel])[0,1]), float(bal), len(P)

COVSETS={'块数':np.c_[NC],'块数+年龄':np.c_[NC,age],'块数+年龄+勾选数':np.c_[NC,age,np.nan_to_num(PK,nan=np.nanmean(PK))]}
raw=float(np.corrcoef(rho[ii],S[ii])[0,1])
rb=np.random.default_rng(4321)
boot=float(np.std([np.corrcoef(rho[s_],S[s_])[0,1] for s_ in
                   (ii[rb.integers(0,len(ii),len(ii))] for _ in range(200))]))
print(f"\n未匹配 corr(rho,S) = {raw:+.4f}   自助展布 {boot:.4f} -> {abs(raw)/boot:.1f}x\n")
print(f"  {'协变量集':<18} {'匹配后':>9} {'保留':>7} {'块数残差':>9} {'配对数':>7}")
rows=[]
for nm,cov in COVSETS.items():
    vals=[matched_corr(900+s,cov) for s in range(5)]
    v=np.array([x[0] for x in vals]); b=np.mean([x[1] for x in vals]); n=int(np.mean([x[2] for x in vals]))
    rows.append(dict(cov=nm,corr=float(np.nanmean(v)),sd=float(np.nanstd(v)),bal=b,n=n))
    print(f"  {nm:<18} {np.nanmean(v):>+9.4f} {100*np.nanmean(v)/raw:>6.0f}% {b:>9.3f} {n:>7,}")

# 零:同一匹配管线跑在置换数据上
rho_n=rho_of(demean_conv(perm_null(V,np.random.default_rng(6600))))
rho_save=rho.copy(); rho=rho_n
nullv=float(np.nanmean([matched_corr(900+s,COVSETS['块数+年龄'])[0] for s in range(5)]))
rho=rho_save
# 正对照:种植一个已知的 corr(rho,S) 再匹配
x=rar-rar.mean()
Vp=V+1.5*np.outer(np.nan_to_num(S),x)*obs
rho_p=rho_of(demean_conv(Vp)); rho_save=rho; rho=rho_p
plantv=float(np.nanmean([matched_corr(900+s,COVSETS['块数+年龄'])[0] for s in range(5)]))
rho=rho_save
main=[r for r in rows if r['cov']=='块数+年龄'][0]
print(f"\n  置换零 {nullv:+.4f}   种植正对照 {plantv:+.4f}")

g=Gate('corr(rho,S) 掉的一半是心理学还是仪器')
g.asserted('匹配真的把块数差压下去了',main['bal']<0.1,f"块数残差 {main['bal']:.3f} sd")
g.asserted('正对照在匹配后仍被回收',abs(plantv)>3*abs(nullv),
           f"种植 {plantv:+.4f} vs 零 {nullv:+.4f}")
g.require_resolvable_first('匹配后 corr(rho,S) 可分辨',abs(main['corr']-nullv),boot)
g.negative_control('匹配后 corr(rho,S) 对置换零',nullv,main['corr'],null_spread=boot)
g.no_sign_crossing('匹配前后不换号',[raw,main['corr']])
print(g)
D_=pd.DataFrame(rows); D_.to_csv(pathlib.Path(__file__).parent/'results'/'matched.csv',index=False)
print(f"\nartifact sha1 {hashlib.sha1(D_.to_csv(index=False).encode()).hexdigest()[:12]}")
