import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A48 R255 -- 用 `#207` 自己的代码,在 onset ≤ 14 下重跑

`#209b`:我的**重实现**给出 `r(rho_i, age)`:全部 +0.1158 -> 仅 ≤14 **−0.0330**(反号)。
`#209c`:但重实现给出的是 **+0.1158**,而 `#206b`/`#207c` 报的是 **+0.153** ——
**一个用重实现得到的反例,只能证伪"这个现象",不能证伪"那个数"。**
本轮还这笔债:**载入 `R173` 的 `betas`/`demean_conv`(与 `#207` 同一条管道),
唯一改动是把 `V` 换成 `where(V ≤ 14, V, nan)`。**

ESTIMAND        用原代码算 `mean rho_i` 与 `r(rho_i, age)`,全样本 vs onset ≤ 14。
KILL            **若用原代码也反号 -> `#206b`/`#207c` 的 `+0.153` 正式降级为截断假象,
                而 `#199` 的「第三个维度」需重查它有多少建立在那条关系上。**
可判前提        **未截断时必须复现 `+0.153`** —— 复现不了就说明我载入的不是同一条管道,整轮不可读。
POWER           `KEEP` 要求 ≥8 个起始年龄;截断后大量人不满足 —— **必须报剩多少人**,
                并报受限臂的 MDE。
NEGATIVE CTRL   人内打乱起始年龄。
IMPOSSIBLE      截断后剩下的人**不是随机子集** —— 他们是「14 岁前就获得 ≥8 个类别」的人,
                本身更早熟。所以受限臂测的是**这群人**,外推受限。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
age=d['age'].map(AGE).values.astype(float)
V_full=V.copy(); obs_full=obs.copy(); KEEP_full=KEEP.copy()
print(f"原管道载入:V {V.shape} · KEEP {int(KEEP.sum()):,}")

def run(Vm, tag):
    global V, obs, KEEP
    V=Vm; obs=np.isfinite(Vm); KEEP=(obs.sum(1)>=8)
    _,rho=betas(V)
    m=np.isfinite(rho)&KEEP&np.isfinite(age)
    r=float(np.corrcoef(rho[m],age[m])[0,1])
    rb=np.random.default_rng(20260803)
    bs=[float(np.corrcoef(rho[i],age[i])[0,1]) for i in
        [rb.choice(np.flatnonzero(m),int(m.sum()),replace=True) for _ in range(300)]]
    out=dict(arm=tag,n_keep=int(KEEP.sum()),n_used=int(m.sum()),
             mean_rho=float(np.mean(rho[m])),r_age=r,sd=float(np.std(bs)))
    V=V_full; obs=obs_full; KEEP=KEEP_full
    return out, rho

o_full,rho_full=run(V_full,'全部起始年龄')
V_cap=np.where(V_full<=14,V_full,np.nan)
o_cap,rho_cap=run(V_cap,'仅 onset ≤ 14')
rows=[o_full,o_cap]
# 负对照
rb=np.random.default_rng(31)
o_null,_=run(perm_null(V_full,rb),'负对照:人内跨人置换(全样本)')
rows.append(o_null)
T=pd.DataFrame(rows); check_columns(T,'R255'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rerun.csv',index=False)
print(f"\n{'臂':<26}{'KEEP':>8}{'用到':>8}{'mean rho':>11}{'r(rho,age)':>12}{'sd':>9}{'比':>7}")
for _,r in T.iterrows():
    print(f"{r.arm:<26}{r.n_keep:>8,}{r.n_used:>8,}{r.mean_rho:>+11.4f}{r.r_age:>+12.4f}"
          f"{r.sd:>9.4f}{abs(r.r_age)/r.sd:>7.1f}")
full=T.iloc[0]; cap=T.iloc[1]; nul=T.iloc[2]
print(f"\n受限臂 MDE = {2*cap.sd:.4f};截断后 KEEP 从 {full.n_keep:,} 掉到 {cap.n_keep:,}"
      f"({100*cap.n_keep/full.n_keep:.0f}%)")

g=Gate('用原代码,截断能不能解释 #207c')
g.asserted('可判前提:未截断时复现 `#206b`/`#207c` 的 +0.153',abs(full.r_age-0.153)<0.02,
           f"{full.r_age:+.4f}")
g.negative_control('人内跨人置换(全样本)',abs(float(nul.r_age)),abs(float(full.r_age)))
g.resolvable('受限臂的 r(rho_i, age)',float(cap.r_age),float(cap.sd))
g.offset_control('受限 vs 全样本',float(cap.r_age),float(full.r_age),float(np.hypot(cap.sd,full.sd)),
                 null_kind='同一条管道在全部起始年龄上的 r(rho_i, age) —— 不是零假设,'
                           '是"若截断无关,受限臂该落在哪"')
g.asserted('注册的 kill:用原代码也反号 -> `+0.153` 降级为截断假象',
           np.sign(cap.r_age)!=np.sign(full.r_age),
           f"全样本 {full.r_age:+.4f} -> 受限 {cap.r_age:+.4f}")
g.asserted('⚠ 受限后剩下的人不是随机子集(他们 14 岁前就获得 ≥8 个类别,更早熟)',True,
           f"KEEP {full.n_keep:,} -> {cap.n_keep:,}({100*cap.n_keep/full.n_keep:.0f}%)")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
