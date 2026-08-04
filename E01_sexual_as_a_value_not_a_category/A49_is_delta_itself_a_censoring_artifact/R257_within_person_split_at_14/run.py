import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A49 R257 -- 按绝对年龄 14 人内劈开:发育窗口,还是截断假象

`#211d`:两种读法绑在一起,而**它们可以被分开** ——
    ② 发育窗口:同一批人里,「14 岁后那一半」内部**仍有** Δ,「14 岁前那一半」内部**没有**
    ① 截断假象:两半内部**都没有**
**人内劈分不需要跨人比较,所以不受截断影响。**

⚠ `#208b`③ 的教训直接烤进前提:**两半的信度必须都 ≥ 下限**,
**不能只要求"可比"** —— 两个 0.000 也是可比的,而那次它 PASS 了。

ESTIMAND        按**绝对年龄 14** 劈(不是中位数),只取两半各 ≥6 个类别的人;
                每半内部各算一次 `mean rho_i`。
KILL            条件式:先要**两半的超额方差都 ≥ 0.02**(否则比的是噪声,不判);
                再判:**14 岁后那一半有 Δ 而 14 岁前没有 -> 发育窗口;两半都没有 -> 截断假象。**
NEGATIVE CTRL   人内打乱起始年龄(保留曲目与年龄分布,毁掉配对)。
POSITIVE CTRL   只往「14 岁后那一半」种一个 Δ,强度扫到推得动为止(`#211a` 的教训)。
IMPOSSIBLE      **「14 岁后那一半」主要属于年长受访者**(15 岁的人几乎没有 14 岁后的条目)——
                所以这一半本身是年龄选择的。必须报两半贡献者的年龄分布。
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
V0=V.copy(); rar0=rar.copy()

def halves(Vm, shuffle=False, rng=None, plant_late=0.0, u=None, need=6):
    """按绝对年龄 14 劈;每半内部各做一次双向去均值再算 rho。"""
    out={}
    for tag,mask in (('≤14',Vm<=14),('>14',Vm>14)):
        A=np.where(np.isfinite(Vm)&mask,Vm,np.nan)
        if shuffle:
            for i in range(len(A)):
                idx=np.flatnonzero(np.isfinite(A[i]))
                if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
        if plant_late and tag=='>14' and u is not None:
            A=A+plant_late*np.outer(u,rar0-rar0.mean())*np.isfinite(A)
        D=np.where(np.isfinite(A),A,np.nan)
        for _ in range(200):
            a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
            b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        rho=np.full(len(D),np.nan)
        for i in range(len(D)):
            idx=np.flatnonzero(np.isfinite(D[i]))
            if len(idx)<need: continue
            x=rar0[idx]-rar0[idx].mean(); y=D[i,idx]
            if x.std()<1e-9 or np.nanstd(y)<1e-9: continue
            rho[i]=float(np.corrcoef(y,x)[0,1])
        out[tag]=rho
    return out

rng=np.random.default_rng(20260803)
H=halves(V0); Hn=halves(V0,shuffle=True,rng=np.random.default_rng(55))
both=np.isfinite(H['≤14'])&np.isfinite(H['>14'])
print(f"两半各 ≥6 个类别的人:{int(both.sum()):,}")
for tag in ('≤14','>14'):
    m=np.isfinite(H[tag])
    aa=age[m&np.isfinite(age)]
    print(f"  {tag:<4} n={int(m.sum()):>6,}  贡献者年龄中位 {np.median(aa):.1f} · "
          f"14–17 段占 {100*np.mean(aa==15.5):.0f}%")

rows=[]
for tag in ('≤14','>14'):
    m=np.isfinite(H[tag]); mn=np.isfinite(Hn[tag])
    rel=max(np.var(H[tag][m])-np.var(Hn[tag][mn]),0)/np.var(H[tag][m])
    bs=[float(np.mean(H[tag][i])) for i in
        [rng.choice(np.flatnonzero(m),int(m.sum()),replace=True) for _ in range(300)]]
    rows.append(dict(half=tag,n=int(m.sum()),delta=float(np.mean(H[tag][m])),sd=float(np.std(bs)),
                     rel=float(rel),null_delta=float(np.mean(Hn[tag][mn]))))
T=pd.DataFrame(rows); check_columns(T,'R257'); T.to_csv(pathlib.Path(__file__).parent/'results'/'split14.csv',index=False)
print(f"\n{'半':<6}{'n':>8}{'Δ':>11}{'sd':>9}{'比':>7}{'信度(超额方差比)':>16}{'打乱 Δ':>10}")
for _,r in T.iterrows():
    print(f"{r.half:<6}{r.n:>8,}{r.delta:>+11.4f}{r.sd:>9.4f}{abs(r.delta)/r.sd:>7.1f}{r.rel:>16.4f}{r.null_delta:>+10.4f}")

# 正对照:只往 >14 半种,强度扫到推得动(`#211a`)
u=np.abs(rng.standard_normal(len(V0)))+0.5; G=None
base=float(T[T.half=='>14'].delta.iloc[0]); sdb=float(T[T.half=='>14'].sd.iloc[0])
for g_ in [1,3,10,30,100]:
    Hp=halves(V0,plant_late=g_,u=u); mp=np.isfinite(Hp['>14'])
    if abs(np.mean(Hp['>14'][mp])-base)>3*sdb: G=g_; break
print(f"\n正对照强度扫描 -> g = {G}")
Hp=halves(V0,plant_late=G,u=u) if G else None
if Hp is not None:
    for tag in ('≤14','>14'):
        mp=np.isfinite(Hp[tag]); print(f"  种植后 {tag}: Δ = {np.mean(Hp[tag][mp]):+.4f}")

lo=T[T.half=='≤14'].iloc[0]; hi=T[T.half=='>14'].iloc[0]
g=Gate('发育窗口还是截断假象')
g.asserted('可判前提(`#208b`③):两半的信度**都** ≥ 0.02 —— 不能只要求"可比"',
           (lo.rel>=0.02) and (hi.rel>=0.02),
           f"≤14 {lo.rel:.4f} · >14 {hi.rel:.4f}")
g.asserted('正对照:只往 >14 半种植,必须推得动',G is not None,f"g = {G}")
g.negative_control('人内打乱(>14 半)',abs(float(hi.null_delta)),abs(float(hi.delta)))
g.resolvable('>14 半的 Δ',float(hi.delta),float(hi.sd))
g.resolvable('≤14 半的 Δ',float(lo.delta),float(lo.sd))
g.asserted('⚠ >14 半本身是年龄选择的(15 岁的人几乎没有 14 岁后的条目)',True,
           f"贡献者年龄中位 ≤14 半 vs >14 半,见上表")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
