import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A85 R327 -- 被吸收剩下的那一成,和年龄

`#281b`:`rate⊥rho` 仍打中 5/29,最强格是 **`age` −0.072**。
而年龄在这条线上有一段历史:`#206b`/`#207c` 报过 `rho_i × age = +0.1563`,
**`#210` 把它撤回为截断假象**,`#246b` 又实测「当前年龄在五个钟的联合模型里只剩 +0.0226」。
**现在 `rate⊥rho × age` 是 −0.072,符号与 `#206b` 相反。**

ESTIMAND        ① `corr(rate⊥rho, age)` 的自助展布与两层区间;
                ② 用 `#210` 的**截断控制**(`V → where(V ≤ 14)`)判它是不是同一个截断假象,
                   **并按守卫 12 在交集样本上重报**(`#240a` 的做法);
                ③ 同报 `corr(rate, age)` 与 `corr(rho_i, age)`,三者一张表。
KILL            **若截断控制之后 −0.072 消失 -> 它与 `#206b` 是同一个假象,只是符号相反;
                若存活 -> 这条线上第一个在截断控制后仍站着的年龄关系。**
⚠ 守卫 12        截断控制**天然改变纳入**,必须在交集样本上重报。
NEGATIVE CTRL   跨人置换年龄(**只在有限值内**)。
POSITIVE CTRL   复现 `#210` 的两个已发表读数(`rho×age` 全样本 +0.1532 / 截断 −0.0741)——
                否则不是同一条管道。
IMPOSSIBLE      `age` 是**当前**年龄,分箱后取中点;它同时携带世代与观测窗口。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_of(Vm):
    Dm=np.where(np.isfinite(Vm),Vm,np.nan)
    for _ in range(300):
        a=np.nanmean(Dm,0,keepdims=True); Dm=Dm-np.where(np.isfinite(a),a,0)
        b=np.nanmean(Dm,1,keepdims=True); Dm=Dm-np.where(np.isfinite(b),b,0)
    W=np.isfinite(Dm); Z=np.where(W,Dm,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); okm=(k>=8)&(den>1e-12); out[okm]=num[okm]/den[okm]; return out
def quantities(cap=None):
    Vm=V0 if cap is None else np.where(V0<=cap,V0,np.nan)
    rho=rho_of(Vm); nc=np.isfinite(Vm).sum(1)
    Vs=np.sort(np.where(np.isfinite(Vm),Vm,np.inf),axis=1)
    a1=Vs[:,0]; acl=np.nanmax(np.where(np.isfinite(Vm),Vm,np.nan),axis=1)
    rate=np.where(np.isfinite(rho),(nc-1)/np.maximum(acl-a1,0.5),np.nan)
    m=np.isfinite(rate)&np.isfinite(rho)
    X=np.column_stack([np.ones(m.sum()),rho[m]]); res=np.full(N,np.nan)
    res[m]=rate[m]-X@np.linalg.lstsq(X,rate[m],rcond=None)[0]
    return rho,rate,res
rngB=np.random.default_rng(20260804)
def cr(a,b,mask=None):
    m=np.isfinite(a)&np.isfinite(b)
    if mask is not None: m&=mask
    if m.sum()<300: return np.nan,np.nan,0
    r=float(np.corrcoef(a[m],b[m])[0,1])
    sd=float(np.std([np.corrcoef(a[i],b[i])[0,1] for i in
        (rngB.choice(np.flatnonzero(m),int(m.sum()),True) for _ in range(300))]))
    return r,sd,int(m.sum())
RHO,RATE,RES=quantities()
RHOc,RATEc,RESc=quantities(cap=14)
print(f"正对照:复现 `#210` 的两个读数 —— ")
r_full,_,n_full=cr(RHO,age); r_cap,_,n_cap=cr(RHOc,age)
print(f"   `rho×age` 全样本 **{r_full:+.4f}**(n={n_full:,};`#210` 报 +0.1532)· "
      f"截断 ≤14 **{r_cap:+.4f}**(n={n_cap:,};报 −0.0741)")
rows=[]
for nm,v,vc in (('rho_i',RHO,RHOc),('rate',RATE,RATEc),('**rate⊥rho**',RES,RESc)):
    r0,s0,n0=cr(v,age); r1,s1,n1=cr(vc,age)
    inter=np.isfinite(v)&np.isfinite(vc)&np.isfinite(age)
    r0i,_,ni=cr(v,age,inter); r1i,_,_=cr(vc,age,inter)
    rows.append(dict(quantity=nm,r_full=r0,sd_full=s0,n_full=n0,r_cap=r1,sd_cap=s1,n_cap=n1,
                     r_full_inter=r0i,r_cap_inter=r1i,n_inter=ni))
    print(f"\n{nm}:")
    print(f"   全样本 **{r0:+.4f} ± {s0:.4f}**(n={n0:,})· 截断 ≤14 **{r1:+.4f} ± {s1:.4f}**(n={n1:,})")
    print(f"   ⚠ 交集样本(n={ni:,}):未截断 **{r0i:+.4f}** -> 截断 **{r1i:+.4f}**")
T=pd.DataFrame(rows); check_columns(T,'R327')
T.to_csv(pathlib.Path(__file__).parent/'results'/'age_table.csv',index=False)
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cr(RES,perm_finite(age,100+i))[0] for i in range(20)]
print(f"\n负对照(置换年龄,只在有限值内):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
r=T[T.quantity=='**rate⊥rho**'].iloc[0]
g=Gate('被吸收剩下的那一成和年龄')
g.asserted('正对照:复现 `#210` 的两个已发表读数(否则不是同一条管道)',
           abs(r_full-0.1532)<0.02 and abs(r_cap+0.0741)<0.02,
           f"全样本 {r_full:+.4f} vs +0.1532;截断 {r_cap:+.4f} vs −0.0741")
g.negative_control('置换年龄',abs(float(np.mean(nul))),abs(float(r.r_full)),
                   null_spread=float(np.std(nul)),null_kind='跨人置换年龄(只在有限值内)—— 只打掉配对')
g.has_error_bar('`corr(rate⊥rho, age)` 全样本',float(r.r_full),float(r.sd_full),'bootstrap_人层')
g.control_kept_the_sample('★ 截断控制(`rate⊥rho × age`)',
                          before=float(r.r_full),after=float(r.r_cap),
                          n_before=int(r.n_full),n_after=int(r.n_cap),
                          before_common=float(r.r_full_inter),after_common=float(r.r_cap_inter),
                          n_common=int(r.n_inter))
g.asserted('★ 注册的 kill:截断控制后 −0.072 消失 -> 同一个假象;存活 -> 第一个站住的年龄关系',
           abs(float(r.r_cap_inter))>2*float(r.sd_cap),
           f"交集样本上 未截断 {r.r_full_inter:+.4f} -> 截断 {r.r_cap_inter:+.4f} "
           f"(展布 ±{r.sd_cap:.4f})")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
