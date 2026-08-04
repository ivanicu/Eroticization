import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A54 R269 -- 用更灵敏的那把刀,重切那条被撤回过的关系

`#223b`:本项目从 `#128` 起一直用 `rho_i`,而 `u_i` 在同一批结局上普遍大 20%–100%。
`#210` 正是**用 `rho_i`** 判定「年龄关系是截断假象」的。
**一个被撤回的关系,值得用更灵敏的仪器复查一次** —— 两种结果都有意义:

    复现   受限臂仍反号 -> `#210` 的撤回被**独立复现**,可升 D8
    重开   受限臂在 `u_i` 上不再反号 -> 那次撤回是**仪器不够灵敏**造成的,要重开

ESTIMAND        `r(u_i, age)` 在 ① 全部起始年龄 ② `onset ≤ 14` 两档上,以及偏移。
KILL            **若受限臂在 `u_i` 上不再反号 -> 重开 `#210`。**
可判前提        ① 全样本上 `u_i × age` 复现 `#223a` 的 +0.1818;
                ② `u_i` 在受限臂上确实比 `rho_i` 更灵敏(否则"更灵敏"这个前提不成立)。
NEGATIVE CTRL   人内跨人置换。
IMPOSSIBLE      受限后剩下的人不是随机子集(`#210c`)—— 足以证伪,不足以立反向声明。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=d['age'].map(AGE).values.astype(float)

def demean_np(Aa,iters=200,tol=1e-10):
    D=np.where(np.isfinite(Aa),Aa,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def both_scores(Aa,seed=0,iters=300):
    """返回 (u_i, rho_i) —— 同一个矩阵、同一个 KEEP,两把刀并排。"""
    keep=(np.isfinite(Aa).sum(1)>=8)
    D=demean_np(Aa); W=np.isfinite(D)&keep[:,None]; Z=np.where(W,D,0.0)
    # rho_i(固定 x = 稀有度)
    k=W.sum(1); rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); Y0=np.where(W,D,0.0)
    yb=np.where(k>0,Y0.sum(1)/np.maximum(k,1),0.0); Yc=W*(Y0-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    rho=np.full(len(D),np.nan); ok=(k>=8)&(den>1e-12)&keep; rho[ok]=num[ok]/den[ok]
    # u_i(x 自由)
    rng=np.random.default_rng(seed); x=rng.standard_normal(Mc)
    for _ in range(iters):
        Xf=W*x[None,:]; dn=(Xf*Xf).sum(1); u=np.where(dn>1e-12,(Z*Xf).sum(1)/np.maximum(dn,1e-12),0.0)
        Uc=W*u[:,None]; dn2=(Uc*Uc).sum(0); x=np.where(dn2>1e-12,(Z*Uc).sum(0)/np.maximum(dn2,1e-12),0.0)
        n=np.linalg.norm(x)
        if n>0: x=x/n
    Xf=W*x[None,:]; dn=(Xf*Xf).sum(1)
    U=np.where(dn>1e-12,(Z*Xf).sum(1)/np.maximum(dn,1e-12),np.nan); U[~(keep&(dn>1e-12))]=np.nan
    if np.corrcoef(x,rar0)[0,1]<0: U=-U
    return U,rho,int(keep.sum())

def rr(x,y):
    m=np.isfinite(x)&np.isfinite(y); return float(np.corrcoef(x[m],y[m])[0,1]), int(m.sum())
rng=np.random.default_rng(20260803)
rows=[]
for tag,Vm in (('全部起始年龄',V0),('仅 onset ≤ 14',np.where(V0<=14,V0,np.nan))):
    U,R,nk=both_scores(Vm)
    ru,nu=rr(U,age); rr_,_=rr(R,age)
    m=np.isfinite(U)&np.isfinite(age); idx=np.flatnonzero(m)
    sd=float(np.std([rr(U[i],age[i])[0] for i in
        [rng.choice(idx,len(idx),replace=True) for _ in range(200)]]))
    rows.append(dict(arm=tag,keep=nk,n=nu,r_u=ru,sd_u=sd,r_rho=rr_))
    print(f"{tag:<16} KEEP {nk:>6,}  r(u_i, age) = {ru:+.4f} ± {sd:.4f}  ·  r(rho_i, age) = {rr_:+.4f}")
Un,Rn,_=both_scores(perm_null(V0,np.random.default_rng(31)))
r_null,_=rr(Un,age)
print(f"{'负对照 题内跨人':<16} {'':>6}   r(u_i, age) = {r_null:+.4f}")

T=pd.DataFrame(rows); check_columns(T,'R269'); T.to_csv(pathlib.Path(__file__).parent/'results'/'recheck.csv',index=False)
full=T.iloc[0]; cap=T.iloc[1]
g=Gate('用更灵敏的刀,那条撤回还站不站得住')
g.asserted('可判前提一:全样本上 u_i × age 复现 `#223a` 的 +0.1818',abs(full.r_u-0.1818)<0.02,
           f"{full.r_u:+.4f}")
g.asserted('可判前提二:u_i 在两档上都比 rho_i 灵敏(否则"更灵敏"不成立)',
           abs(full.r_u)>abs(full.r_rho) and abs(cap.r_u)>=abs(cap.r_rho)*0.9,
           f"全样本 {abs(full.r_u):.4f} vs {abs(full.r_rho):.4f} · 受限 {abs(cap.r_u):.4f} vs {abs(cap.r_rho):.4f}")
g.negative_control('题内跨人置换',abs(r_null),abs(float(full.r_u)),null_kind='题内跨人置换')
g.resolvable('受限臂的 r(u_i, age)',float(cap.r_u),float(cap.sd_u))
g.offset_control('受限 vs 全样本(u_i)',float(cap.r_u),float(full.r_u),
                 float(np.hypot(cap.sd_u,full.sd_u)),
                 null_kind='同一条管道在全部起始年龄上的 r(u_i, age) —— 若截断无关,受限臂该落在哪')
g.asserted('注册的 kill:受限臂在 u_i 上不再反号 -> 重开 `#210`',
           np.sign(cap.r_u)==np.sign(full.r_u),
           f"全样本 {full.r_u:+.4f} -> 受限 {cap.r_u:+.4f}")
g.asserted('⚠ 受限后的人不是随机子集 —— 足以证伪,不足以立反向声明(`#210c`)',True,
           f"KEEP {int(full.keep):,} -> {int(cap.keep):,}")
print(g)
print(f"\n  => {'重开 `#210`' if np.sign(cap.r_u)==np.sign(full.r_u) else '`#210` 的撤回被独立复现'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
