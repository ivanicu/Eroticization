import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A42 R244 -- 「何时」的人际差异,是不是「什么」的影子

`#198a`:`rho_i` 有真实的人际方差(超额 sd 0.176,49×),**而它从没被拿去对过任何外部变量。**
`#159` 说「何时」那点结构几乎全在类别上(题目侧 0.96 vs 人侧 0.105)——
**但那是在说 Δ 的平均,不是在说它的人际差异。**

    SHADOW  rho_i 只与 S(位置分)相关,20 道结局全无 -> 「何时」的人际差异是「什么」的影子
    THIRD   rho_i 自己越全族阈值 -> 这份数据里存在**第三个**独立的人层维度

ESTIMAND        rho_i 对 20 道 Likert 的相关(最大统计量零给全族阈值),以及 r(rho_i, S)。
                **同时报去衰减值** —— `#198a` 给了 rho_i 一个真信度:
                `rel = var_超额 / var_总 = 0.176²/0.357² ≈ 0.242`。
KILL            条件式:先要**正对照开火**(种一个人层信号进 rho_i,面板必须测到);
                再判:**rho_i 自己越阈值的题数 == 0 -> SHADOW;>0 -> THIRD。**
NEGATIVE CTRL   每题在分析样本内打乱(`#184b` 的教训)。
NOISE FLOOR     人层 bootstrap 200。
IMPOSSIBLE      rel ≈ 0.24 —— **rho_i 是一个信度很低的量**,所以"越不过阈值"有很大一部分
                是衰减造成的。**去衰减值必须并排报,而它带着 rel 估计本身的误差。**
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_coverage

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df_raw.columns if df_raw[c].dtype!=object and
     set(pd.Series(df_raw[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df_raw[c].notna().sum()>10000]
print(f"载入:KEEP {int(KEEP.sum()):,} 人 · Likert {len(lik)} 道",flush=True)

_,RHO=betas(V)
rng=np.random.default_rng(20260803)
var_null=float(np.mean([np.nanvar(betas(perm_null(V,np.random.default_rng(4200+s)))[1][np.isfinite(RHO)&KEEP])
                        for s in range(5)]))
m0=np.isfinite(RHO)&KEEP
var_tot=float(np.var(RHO[m0])); rel=max(var_tot-var_null,0)/var_tot
print(f"rho_i 信度 rel = var_超额/var_总 = {rel:.4f}  -> 去衰减因子 1/√rel = {1/np.sqrt(rel):.2f}×")

base=m0&np.isfinite(S)
r_S=float(np.corrcoef(RHO[base],S[base])[0,1])
print(f"r(rho_i, S) = {r_S:+.4f}(n={int(base.sum()):,})· 去衰减 {r_S/np.sqrt(rel*0.432):+.4f}"
      f"   [S 的信度 0.432,`#100`]")

def cr(y,x,ii):
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    return float(np.corrcoef(y[jj],x[jj])[0,1]), len(jj)
bi=np.flatnonzero(m0)
rows=[]; nulls=[]
for c in lik:
    y=df_raw[c].values.astype(float)
    r,n=cr(y,RHO,bi)
    sd=float(np.std([cr(y,RHO,rng.choice(bi,len(bi),replace=True))[0] for _ in range(200)]))
    ps=[]
    for _ in range(40):
        yp=y.copy(); yp[bi]=rng.permutation(y[bi]); v,_=cr(yp,RHO,bi)
        if np.isfinite(v): ps.append(abs(v))
    if len(ps)>=20: nulls.append(ps)
    rows.append(dict(q=c[:60],n=n,r_rho=r,r_disatt=r/np.sqrt(rel),sd=sd,ratio=abs(r)/sd))
T=pd.DataFrame(rows); check_columns(T,'R244'); check_coverage(len(T),len(lik),'R244 面板',tol=0.10)
L=min(len(x) for x in nulls)
thr=float(np.nanquantile(np.nanmax(np.array([x[:L] for x in nulls]),axis=0),0.95))
T=T.sort_values('r_rho',key=abs,ascending=False)
T.to_csv(pathlib.Path(__file__).parent/'results'/'panel.csv',index=False)
print(f"\n全族阈值 |r| = {thr:.4f}\n")
print(f"{'r(rho,·)':>10}{'去衰减':>10}{'比':>7}  题")
for _,r in T.head(8).iterrows():
    print(f"{r.r_rho:>+10.4f}{r.r_disatt:>+10.4f}{r.ratio:>7.1f}{'★' if abs(r.r_rho)>thr else ' '}  {r.q[:54]}")
n_pass=int((T.r_rho.abs()>thr).sum()); print(f"\n越阈值 {n_pass}/{len(T)}")

# 正对照:种一个人层信号进 rho_i
u=rng.standard_normal(len(V)); Vp=V+2.0*np.outer(np.abs(u)+0.5,rar-rar.mean())*obs
_,RHOp=betas(Vp); mp=np.isfinite(RHOp)&KEEP
r_plant=float(np.corrcoef(RHOp[mp],(np.abs(u)+0.5)[mp])[0,1])
print(f"正对照(种一个人层信号):r(rho_种植, u) = {r_plant:+.4f}")

g=Gate('「何时」的人际差异是不是「什么」的影子')
g.asserted('正对照:种进去的人层信号必须被 rho_i 测到',abs(r_plant)>0.3,f"{r_plant:+.4f}")
g.asserted('⚠ rho_i 信度很低,已量化并并排报去衰减值',rel<0.5,
           f"rel = {rel:.4f} -> 去衰减因子 {1/np.sqrt(rel):.2f}×")
g.resolvable('r(rho_i, S)',r_S,float(np.std([np.corrcoef(RHO[b],S[b])[0,1] for b in
             [rng.choice(np.flatnonzero(base),int(base.sum()),replace=True) for _ in range(200)]])))
g.asserted('注册的 kill:rho_i 自己越阈值的题数 > 0 -> 第三个维度',n_pass>0,
           f"越阈值 {n_pass}/{len(T)}(阈值 {thr:.4f})")
print(g)
print(f"\n  => {'THIRD —— 存在第三个人层维度' if n_pass>0 else 'SHADOW —— 「何时」的人际差异对外部无话可说'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
