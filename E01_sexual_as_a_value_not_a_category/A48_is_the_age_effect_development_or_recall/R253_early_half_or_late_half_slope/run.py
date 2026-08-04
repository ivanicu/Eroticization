import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A48 R253 -- 年纪越大罕见兴趣来得越晚:发育,还是回忆窗口

`#207c`:`rho_i × age = +0.153` 是真的,而它**不中介** rho_i 对 20 道题的预测。
但那句观察本身有两种完全不同的读法:
    DEVELOPMENT  年长者确实在**更晚的年纪**才获得罕见兴趣 -> 只有**晚期一半**与 age 相关
    RECALL       年长者回溯更远,**整条曲线被重排** -> **两半同等**与 age 相关
                 -> `#207c` 那句独立观察降级成测量假象

`#162` 已证回忆偏差是**人群规律**;`#183` 已证 age 不与**截距**相关(+0.006)——
**所以两种读法在截距上没有区别,只在斜率上有。**

ESTIMAND        人内按起始年龄中位数劈半,**每半各算一次「稀有度 × 起始年龄残差」的相关**;
                判 `r(age, 早半)` 与 `r(age, 晚半)` 的**配对差**。
KILL            条件式:先要**两半信度可比**(否则"哪半更强"说的是仪器,`#180` 的教训);
                再判:**|配对差| > 2× 自身展布且晚半更强 -> DEVELOPMENT;
                两半不可分辨 -> RECALL,`#207c` 降级。**
NEGATIVE CTRL   人内打乱起始年龄(保留曲目与年龄分布,毁掉配对)。
POSITIVE CTRL   只往**晚期一半**种一个与 age 相关的信号 -> 配对差必须指向晚半。
NOISE FLOOR     配对 bootstrap 400。
IMPOSSIBLE      两半各只有 4–15 个类别,**每半的斜率噪声都很大**;
                所以"不可分辨"很可能是功效不足 —— 必须报 MDE。
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
print(f"V {V.shape} · KEEP {int(KEEP.sum()):,}")

def half_rhos(Vm, shuffle=False, rng=None, plant_late=0.0, u=None):
    """人内按年龄中位数劈半,每半各算一次 corr(残差, 稀有度)。"""
    D=demean_conv(Vm)
    n=len(D); e=np.full(n,np.nan); l=np.full(n,np.nan); ne=np.zeros(n); nl=np.zeros(n)
    for i in np.flatnonzero(KEEP):
        idx=np.flatnonzero(obs[i]); a=Vm[i,idx].astype(float)
        if shuffle: a=a[rng.permutation(len(a))]
        o=np.argsort(a,kind='stable'); h=len(o)//2
        if h<3: continue
        E=idx[o[:h]]; L=idx[o[-h:]]
        for arr,sel,cnt in ((e,E,ne),(l,L,nl)):
            y=D[i,sel].astype(float).copy()
            if plant_late and arr is l and u is not None: y=y+plant_late*u[i]*(rar[sel]-rar[sel].mean())
            x=rar[sel]-rar[sel].mean()
            if x.std()<1e-9 or np.nanstd(y)<1e-9: continue
            arr[i]=float(np.corrcoef(y,x)[0,1]); cnt[i]=len(sel)
    return e,l,ne,nl

E,L,NE,NL=half_rhos(V)
m=np.isfinite(E)&np.isfinite(L)&np.isfinite(age)&KEEP; ii=np.flatnonzero(m)
print(f"两半都算得出的人 {len(ii):,};每半类别数 中位 早 {np.median(NE[ii]):.0f} · 晚 {np.median(NL[ii]):.0f}")

def rr(x,ii_): return float(np.corrcoef(age[ii_],x[ii_])[0,1])
r_e,r_l=rr(E,ii),rr(L,ii)
rng=np.random.default_rng(20260803)
dif=[]
for _ in range(400):
    s=rng.choice(ii,len(ii),replace=True); dif.append(rr(L,s)-rr(E,s))
dm,dsd=float(np.mean(dif)),float(np.std(dif))
print(f"\nr(age, 早半) = {r_e:+.4f}   r(age, 晚半) = {r_l:+.4f}")
print(f"配对差(晚 − 早) = {dm:+.4f} ± {dsd:.4f}   |Δ|/sd = {abs(dm)/dsd:.1f}   MDE = {2*dsd:.4f}")

# 两半信度(置换零给噪声宽度,`#198a` 同款)
En,Ln,_,_=half_rhos(V,shuffle=True,rng=np.random.default_rng(77))
mn=np.isfinite(En)&np.isfinite(Ln)&KEEP
rel_e=max(np.var(E[ii])-np.var(En[mn]),0)/np.var(E[ii])
rel_l=max(np.var(L[ii])-np.var(Ln[mn]),0)/np.var(L[ii])
print(f"两半信度(超额方差/总方差):早 {rel_e:.4f} · 晚 {rel_l:.4f}")
# ⚠ #208a:第一版在这里无条件去衰减,rel≈0 时打印出 +53.4 与 +111.3 这种荒谬数。
#   **信度为零时不存在"去衰减值"** —— 那是除以零,不是一个更准的估计。
if min(rel_e,rel_l)<0.02:
    print("  ⚠ **两半的信度都 ≈ 0(超额方差为零)** —— 每半的斜率在人层上是纯噪声。")
    print("     **不打印去衰减值**:rel≈0 时那只是除以零。")
else:
    print(f"去衰减后 r(age,·):早 {r_e/np.sqrt(rel_e):+.4f} · 晚 {r_l/np.sqrt(rel_l):+.4f}")

# 正对照:只往晚半种一个与 age 相关的信号
uA=(age-np.nanmean(age[ii]))/np.nanstd(age[ii]); uA=np.nan_to_num(uA)
Ep,Lp,_,_=half_rhos(V,plant_late=0.8,u=uA)
mp=np.isfinite(Ep)&np.isfinite(Lp)&np.isfinite(age)&KEEP; ip=np.flatnonzero(mp)
print(f"正对照(只往晚半种 age 信号):早 {rr(Ep,ip):+.4f} · 晚 {rr(Lp,ip):+.4f} · "
      f"差 {rr(Lp,ip)-rr(Ep,ip):+.4f}")
# 负对照
mnn=np.isfinite(En)&np.isfinite(Ln)&np.isfinite(age)&KEEP; inn=np.flatnonzero(mnn)
print(f"负对照(人内打乱):早 {rr(En,inn):+.4f} · 晚 {rr(Ln,inn):+.4f} · 差 {rr(Ln,inn)-rr(En,inn):+.4f}")

T=pd.DataFrame([dict(arm='真实',r_early=r_e,r_late=r_l,d=r_l-r_e),
                dict(arm='正对照(晚半种植)',r_early=rr(Ep,ip),r_late=rr(Lp,ip),d=rr(Lp,ip)-rr(Ep,ip)),
                dict(arm='负对照(人内打乱)',r_early=rr(En,inn),r_late=rr(Ln,inn),d=rr(Ln,inn)-rr(En,inn))])
check_columns(T,'R253'); T.to_csv(pathlib.Path(__file__).parent/'results'/'halves.csv',index=False)
g=Gate('年龄效应是发育还是回忆窗口')
g.asserted('正对照:只往晚半种 age 信号,配对差必须指向晚半',
           (rr(Lp,ip)-rr(Ep,ip))>3*dsd,f"{rr(Lp,ip)-rr(Ep,ip):+.4f} vs 3×sd {3*dsd:.4f}")
g.negative_control('人内打乱的配对差',abs(rr(Ln,inn)-rr(En,inn)),abs(dm),null_spread=dsd)
g.asserted('可判前提:两半信度可比(否则"哪半更强"说的是仪器,`#180` 的教训)',
           abs(rel_e-rel_l)<0.15,f"早 {rel_e:.3f} vs 晚 {rel_l:.3f},差 {rel_e-rel_l:+.3f}")
g.resolvable('配对差(晚 − 早)',dm,dsd)
g.asserted('注册的判定:晚半明显更强 -> DEVELOPMENT;不可分辨 -> RECALL',
           abs(dm)>2*dsd and dm>0,f"差 {dm:+.4f} ± {dsd:.4f},MDE {2*dsd:.4f}")
print(g)
# ⚠ #208b:第一版在这里**无条件**打出判定,而负对照刚刚失败(打乱给出的差是真实的 79%),
#   两半信度又都是 0 —— verdict 必须受同一个条件式管辖(`#181a` 的 print 侧同一个 bug,第二次)。
null_d=abs(rr(Ln,inn)-rr(En,inn)); precond = (null_d < 0.5*abs(dm)) and (min(rel_e,rel_l)>=0.02)
if not precond:
    print("\n  => UNVERIFIED。**不判。** 负对照的配对差是真实值的 "
          f"{100*null_d/abs(dm):.0f}%,而两半信度都 ≈ 0 ——")
    print("     「哪一半与年龄绑得更紧」在本设计里比的是两个噪声量。**发育 vs 回忆窗口,本轮答不了。**")
else:
    print(f"\n  => {'DEVELOPMENT' if dm>0 else '晚半更弱 —— 两种读法都不符'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
