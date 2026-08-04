import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A49 R258 -- 转折点在哪里,还是根本没有转折点

`#212b` 把 Δ 定位在 14 岁之后 —— **而 14 是我从截断需要里挑的刻度**
(`#211` 选它是因为最年轻的段是 14–17)。**那个转折点本身可以被估计。**

    INFLECTION  「>c 半的 Δ」随 c 的曲线上有拐点 -> 那是数据自己给出的转折年龄
    SMOOTH      曲线平滑单调 -> 「十四岁」只是我挑的刻度,
                应改写成「越晚的那部分越强」,而不是「十四岁之后才有」

ESTIMAND        劈分点 c ∈ 10…20,每点报两半的 Δ · 信度 · n。
KILL            **若「>c 半的 Δ」随 c 平滑单调、无拐点 -> `#212` 的措辞要从
                「十四岁之后才有」改成「越晚的那部分越强」。**
判据            拐点 = 二阶差分的最大绝对值是否明显高于其余(对着自身展布)。
NEGATIVE CTRL   每点各跑一次人内打乱。
可判前提        每个 c 上两半的信度都要 ≥0.02,否则那一点不入曲线(`#208b`③ / `#212a`)。
IMPOSSIBLE      c 越大,「>c 半」的人越少、越偏向年长者 —— 曲线右端本身在变样本。
                必须同时画 n 与贡献者年龄中位。
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

def one_half(Vm, keep_mask, shuffle=False, rng=None, need=6):
    A=np.where(np.isfinite(Vm)&keep_mask,Vm,np.nan)
    if shuffle:
        for i in range(len(A)):
            idx=np.flatnonzero(np.isfinite(A[i]))
            if len(idx)>1: A[i,idx]=A[i,rng.permutation(idx)]
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
    return rho

rng=np.random.default_rng(20260803); rows=[]
for c in range(10,21):
    out={'c':c}
    for tag,msk in ((f'≤{c}',V0<=c),(f'>{c}',V0>c)):
        r=one_half(V0,msk); rn=one_half(V0,msk,shuffle=True,rng=np.random.default_rng(900+c))
        m=np.isfinite(r); mn=np.isfinite(rn)
        if m.sum()<400: out[tag[0]+'_n']=int(m.sum()); continue
        rel=max(np.var(r[m])-np.var(rn[mn]),0)/np.var(r[m])
        bs=[float(np.mean(r[i])) for i in
            [rng.choice(np.flatnonzero(m),int(m.sum()),replace=True) for _ in range(120)]]
        k='lo' if tag.startswith('≤') else 'hi'
        out[k+'_delta']=float(np.mean(r[m])); out[k+'_sd']=float(np.std(bs))
        out[k+'_rel']=float(rel); out[k+'_n']=int(m.sum())
        out[k+'_age']=float(np.median(age[m&np.isfinite(age)]))
    rows.append(out)
T=pd.DataFrame(rows); check_columns(T,'R258'); T.to_csv(pathlib.Path(__file__).parent/'results'/'sweep.csv',index=False)
print(f"\n{'c':>4}{'≤c 的 Δ':>11}{'信度':>8}{'>c 的 Δ':>11}{'信度':>8}{'>c 比':>8}{'>c n':>8}{'>c 年龄中位':>11}")
for _,r in T.iterrows():
    print(f"{int(r.c):>4}{r.lo_delta:>+11.4f}{r.lo_rel:>8.3f}{r.hi_delta:>+11.4f}{r.hi_rel:>8.3f}"
          f"{abs(r.hi_delta)/r.hi_sd:>8.1f}{int(r.hi_n):>8,}{r.hi_age:>11.1f}")

ok=T[(T.lo_rel>=0.02)&(T.hi_rel>=0.02)]
print(f"\n两半信度都 ≥0.02 的 c:{list(ok.c.astype(int))}")
y=ok.hi_delta.values; cc=ok.c.values
d2=np.diff(y,2) if len(y)>=3 else np.array([np.nan])
sd_d2=float(np.std(d2)) if np.isfinite(d2).all() and len(d2)>1 else np.nan
print(f">c 的 Δ 曲线:{' '.join(f'{v:+.3f}' for v in y)}")
print(f"二阶差分:{' '.join(f'{v:+.4f}' for v in d2)}   最大 |Δ²| = {np.nanmax(np.abs(d2)):.4f} · sd {sd_d2:.4f}")
mono=bool(np.all(np.diff(y)<=1e-9)) or bool(np.all(np.diff(y)>=-1e-9))
g=Gate('转折点在哪里,还是没有转折点')
g.asserted('可判前提:至少 5 个 c 上两半信度都 ≥0.02',len(ok)>=5,f"{len(ok)} 个:{list(ok.c.astype(int))}")
g.resolvable('c=14 处 >c 半的 Δ(复现 `#212b` 的 −0.0704)',
             float(T[T.c==14].hi_delta.iloc[0]),float(T[T.c==14].hi_sd.iloc[0]))
g.no_sign_crossing('>c 半的 Δ 在可判区间内同号',[float(v) for v in y])
g.asserted('注册的 kill:曲线平滑单调无拐点 -> 措辞改成「越晚的那部分越强」',
           mono and (np.nanmax(np.abs(d2))<=2*sd_d2 if np.isfinite(sd_d2) else False),
           f"单调={mono} · 最大二阶差分 {np.nanmax(np.abs(d2)):.4f} vs 2×sd {2*sd_d2:.4f}")
g.asserted('⚠ c 越大 >c 半的人越少越偏年长 —— 曲线右端在变样本',True,
           f"n {int(ok.hi_n.iloc[0]):,} -> {int(ok.hi_n.iloc[-1]):,} · 年龄中位 "
           f"{ok.hi_age.iloc[0]:.1f} -> {ok.hi_age.iloc[-1]:.1f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
