import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A48 R254 -- 不劈半的问法:每个类别的平均起始年龄,随受访者年龄怎么漂

`#208c`:失败在"劈成两半"(每半只剩 6 个类别),不在 `rho_i`。
**同一个问题有一个不需要劈半的问法:**
    DEVELOPMENT  **罕见类别**的平均起始年龄随受访者年龄上升得**更快**
    RECALL       所有类别**同步上升**(仅整体平移)

单位是**类别 × 年龄段**(31 × 5 = 155 格),**不依赖任何人内劈分**。

⚠ **最强混杂,跑之前写下:截断。**
一个 **15 岁**的人不可能报 **25 岁**的起始年龄。所以任何类别的平均起始年龄都会随受访者年龄
**机械上升**,而**典型获得得晚的类别(恰好也是罕见的)升得更多** ——
**这会单独造出 DEVELOPMENT 的图样。**
**控制(同一迭代内)**:平行跑一条**只取起始年龄 ≤ 14** 的臂 ——
14 岁在最年轻的段(14–17)里也可观测,所以那条臂**没有截断差异**。

ESTIMAND        每类别的漂移 = 平均起始年龄对年龄段中点的斜率;判 `corr(漂移, 稀有度)`,n=31。
KILL            **若受限臂(≤14)里 `corr(漂移, 稀有度)` 掉到不可分辨 ->
                全样本的图样是截断造的,DEVELOPMENT 不成立。**
NEGATIVE CTRL   跨类别打乱稀有度标签(2000 次)。
POSITIVE CTRL   人为给罕见类别加一个随受访者年龄增长的漂移 -> 必须被测到。
NOISE FLOOR     类别层 bootstrap 400。
IMPOSSIBLE      n=31 个类别。|r| < 0.36 在 n=31 上不显著。
"""
import numpy as np, pandas as pd, warnings, hashlib, re
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
O=pd.read_csv('data/derived/onset.csv')
ons=[c for c in O.columns if re.search(r'How old were you when you first',c)]
A_=O[ons].apply(pd.to_numeric,errors='coerce').values
A_=np.where((A_>=2)&(A_<=60),A_,np.nan)
assert np.isfinite(A_).sum()>10000
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
band=d['age'].map(AGE).values.astype(float)
have=np.isfinite(A_); rar=-np.log(np.clip(have.mean(0),1e-4,1.))
print(f"类别 {len(ons)} · 年龄段 {sorted(set(AGE.values()))}")

def drift(Amat, cap=None, jitter=None):
    """每类别:平均起始年龄对年龄段中点的斜率。cap = 只取 onset ≤ cap。"""
    out=np.full(Amat.shape[1],np.nan)
    bs=np.array(sorted(set(AGE.values())))
    for j in range(Amat.shape[1]):
        ys=[]
        for b in bs:
            m=np.isfinite(Amat[:,j])&(band==b)
            if cap is not None: m&= (Amat[:,j]<=cap)
            if m.sum()<40: ys.append(np.nan); continue
            v=Amat[m,j]
            if jitter is not None: v=v+jitter*(rar[j]-rar.mean())*(b-bs.mean())
            ys.append(float(np.mean(v)))
        ys=np.array(ys); ok=np.isfinite(ys)
        if ok.sum()>=4: out[j]=float(np.polyfit(bs[ok],ys[ok],1)[0])
    return out

D_all=drift(A_); D_cap=drift(A_,cap=14)
rng=np.random.default_rng(20260803)
def rr(x,y):
    m=np.isfinite(x)&np.isfinite(y); return float(np.corrcoef(x[m],y[m])[0,1]), int(m.sum())
r_all,n_all=rr(D_all,rar); r_cap,n_cap=rr(D_cap,rar)
print(f"\n全样本:平均漂移 {np.nanmean(D_all):+.4f} 岁/受访年 · corr(漂移, 稀有度) = **{r_all:+.4f}**(n={n_all})")
print(f"受限臂(onset ≤ 14,无截断差):平均漂移 {np.nanmean(D_cap):+.4f} · "
      f"corr = **{r_cap:+.4f}**(n={n_cap})")

def boot(D):
    out=[]
    for _ in range(400):
        i=rng.choice(np.flatnonzero(np.isfinite(D)&np.isfinite(rar)),n_all,replace=True)
        out.append(float(np.corrcoef(D[i],rar[i])[0,1]))
    return float(np.std(out))
sd_all,sd_cap=boot(D_all),boot(D_cap)
null=[float(np.corrcoef(D_all[np.isfinite(D_all)],rng.permutation(rar[np.isfinite(D_all)]))[0,1])
      for _ in range(2000)]
print(f"bootstrap sd:全样本 {sd_all:.4f} · 受限 {sd_cap:.4f}")
print(f"负对照(打乱稀有度标签,2000 次):{np.mean(null):+.4f} ± {np.std(null):.4f}")
D_pl=drift(A_,jitter=0.25); r_pl,_=rr(D_pl,rar)
print(f"正对照(给罕见类别加随受访年龄增长的漂移):corr = {r_pl:+.4f}")

T=pd.DataFrame(dict(cat_q=[c[:52] for c in ons],rarity=rar,drift_all=D_all,drift_cap14=D_cap))
check_columns(T,'R254'); T.to_csv(pathlib.Path(__file__).parent/'results'/'drift.csv',index=False)
g=Gate('发育还是回忆窗口(不劈半的问法)')
g.asserted('正对照:人为的稀有度相关漂移必须被测到',r_pl>r_all+0.15,f"{r_all:+.4f} -> {r_pl:+.4f}")
g.negative_control('打乱稀有度标签',float(abs(np.mean(null))),abs(r_all),null_spread=float(np.std(null)))
g.resolvable('全样本 corr(漂移, 稀有度)',r_all,sd_all)
g.asserted('⚠ 截断混杂已在同一迭代内控制',True,
           f"受限臂 onset ≤ 14 —— 14 岁在最年轻的段(14–17)里也可观测,无截断差")
g.asserted('注册的 kill:受限臂里 corr 掉到不可分辨 -> 全样本图样是截断造的',
           abs(r_cap)<=2*sd_cap,f"受限 {r_cap:+.4f} ± {sd_cap:.4f} = {abs(r_cap)/sd_cap:.1f}×")
print(g)
print(f"\n  => {'截断解释了它 —— DEVELOPMENT 不成立' if abs(r_cap)<=2*sd_cap else 'DEVELOPMENT 在无截断差的臂上仍然成立'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 这直接威胁 `#207c`,当场查 -----------------------------------------------
# 受限臂里**平均漂移**从 +0.203 塌到 −0.010 —— 整个"随受访年龄上升"就是截断。
# 而 `#207c` 的 `rho_i × age = +0.153` 是在**全部**起始年龄上算的。
# **同一个截断能不能解释它?** 用只取 onset ≤ 14 的数据重算 rho_i。
print("\n---- 截断能不能解释 `#207c` 的 rho_i × age ----")
def rho_of(Amat, need=6):
    D=np.where(np.isfinite(Amat),Amat,np.nan)
    for _ in range(200):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<1e-10 and np.nanmax(np.abs(b))<1e-10: break
    out=np.full(len(D),np.nan)
    for i in range(len(D)):
        idx=np.flatnonzero(np.isfinite(D[i]))
        if len(idx)<need: continue
        x=rar[idx]-rar[idx].mean(); y=D[i,idx]
        if x.std()<1e-9 or np.nanstd(y)<1e-9: continue
        out[i]=float(np.corrcoef(y,x)[0,1])
    return out
R_full=rho_of(A_)
A_cap=np.where(A_<=14,A_,np.nan)
R_cap=rho_of(A_cap)
for nm,R_ in (('全部起始年龄',R_full),('仅 onset ≤ 14',R_cap)):
    m=np.isfinite(R_)&np.isfinite(band)
    r=float(np.corrcoef(R_[m],band[m])[0,1])
    bs=[float(np.corrcoef(R_[i],band[i])[0,1]) for i in
        [rng.choice(np.flatnonzero(m),int(m.sum()),replace=True) for _ in range(300)]]
    print(f"  {nm:<16} n={int(m.sum()):>6,}  r(rho_i, age) = {r:+.4f} ± {np.std(bs):.4f}  "
          f"{abs(r)/np.std(bs):.1f}×")
    if nm.startswith('全部'): r_f,sd_f=r,float(np.std(bs))
    else: r_c,sd_c=r,float(np.std(bs))
g2=Gate('截断能不能解释 #207c')
g2.asserted('可判前提:全部起始年龄上复现 `#207c` 的 +0.153',abs(r_f-0.153)<0.03,f"{r_f:+.4f}")
g2.resolvable('仅 onset ≤ 14 时的 r(rho_i, age)',r_c,sd_c)
g2.offset_control('受限 vs 全样本',r_c,r_f,float(np.hypot(sd_c,sd_f)),
                  null_kind='同一批人在全部起始年龄上的 r(rho_i, age) —— 不是零假设,'
                            '是"若截断无关,受限臂该落在哪"')
g2.asserted('若受限臂塌到不可分辨 -> `#207c` 是截断假象',abs(r_c)<=2*sd_c,
            f"受限 {r_c:+.4f} ± {sd_c:.4f}")
print(g2)
