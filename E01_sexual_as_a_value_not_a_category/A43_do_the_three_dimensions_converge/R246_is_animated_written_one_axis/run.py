import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A43 R246 -- `animated` 与 `written` 是不是一条轴的两端

`#188` 起,六轮把 `animated`/`written` 读成**「媒介偏好」**(画的 vs 写的)。
**而从没量过它们是不是一条轴的两端。**
`#237` 已给出一个警号:去性别后 `animated` 保留 80%,`written` 只保留 **31%** ——
一条轴的两端不该有这么不同的命运。

    ONE_AXIS  corr(animated, written) 显著为**负** -> "二选一",媒介偏好这个名字成立
    TWO_THING corr 为**正**或近零 -> 它们不是一条轴的两端,`#188` 起的那个名要重估

ESTIMAND        ① `corr(animated, written)` 直接量;② `written` 跑 `#245` 的同一套三预测子分解。
KILL            **若 corr ≥ 0 -> 「媒介偏好」这个名字错了。**
POSITIVE CTRL   一对**确实二选一**的题必须读到负 —— 用 `#239` 数据自己挑出的最强反向对
                (`biomale` × `written`,r = −0.271)作为"负相关长什么样"的刻度。
NEGATIVE CTRL   一道题与自己 = +1(管道自检)。
IMPOSSIBLE      两道题都是"我觉得 X 有情色感"的**正向**李克特,一个人可以两个都高 ——
                **李克特上的"二选一"本来就不该期待强负相关**;所以本轮判的是**符号**,
                不是"够不够负"。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
_,RHO=betas(V)
df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
ANI='animated'; WRI='written'; BIO='biomale'
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df_raw); con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    con[ppl]+=Z@v[:,-1]; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
Cb=np.where(ok,con/np.maximum(cnt,1),np.nan); Sb=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
base=np.isfinite(Sb)&np.isfinite(Cb)&np.isfinite(KB)&np.isfinite(RHO)&KEEP
bi=np.flatnonzero(base)
X0=np.c_[np.ones(len(bi)),Sb[bi]]
Cr=np.full(NN,np.nan); Cr[bi]=Cb[bi]-X0@np.linalg.lstsq(X0,Cb[bi],rcond=None)[0]
check_residualized(Cr[bi],Sb[bi],'R246 内容残差')
P={'S(位置)':Sb,'Cres(内容)':Cr,'rho_i(何时)':RHO}; ks=list(P)

a=df_raw[ANI].values.astype(float); w_=df_raw[WRI].values.astype(float); b_=df_raw[BIO].values.astype(float)
mm=np.isfinite(a)&np.isfinite(w_)
r_aw=float(np.corrcoef(a[mm],w_[mm])[0,1])
mb=np.isfinite(b_)&np.isfinite(w_)
r_bw=float(np.corrcoef(b_[mb],w_[mb])[0,1])
rng=np.random.default_rng(20260803)
sd_aw=float(np.std([np.corrcoef(a[i],w_[i])[0,1] for i in
                    [rng.choice(np.flatnonzero(mm),int(mm.sum()),replace=True) for _ in range(300)]]))
print(f"corr(animated, written) = **{r_aw:+.4f}** ± {sd_aw:.4f}  (n={int(mm.sum()):,})")
print(f"正对照刻度 corr(biomale, written) = {r_bw:+.4f}   <- 「确实二选一」长这样(`#239`)")
print(f"负对照 corr(written, written) = {np.corrcoef(w_[mm],w_[mm])[0,1]:+.4f}")

def partials(y,ii):
    out={}
    for k in ks:
        oth=[P[o][ii] for o in ks if o!=k]
        X=np.c_[np.ones(len(ii)),KB[ii],*oth]
        ry=y[ii]-X@np.linalg.lstsq(X,y[ii],rcond=None)[0]
        rx=P[k][ii]-X@np.linalg.lstsq(X,P[k][ii],rcond=None)[0]
        out[k]=float(np.corrcoef(ry,rx)[0,1])
    return out
rows=[]
for nm,y in [('animated',a),('written',w_)]:
    ii=bi[np.isfinite(y[bi])]; d=partials(y,ii)
    sds={k:float(np.std([partials(y,rng.choice(ii,len(ii),replace=True))[k] for _ in range(150)])) for k in ks}
    rows.append(dict(outcome=nm,n=len(ii),**{k:d[k] for k in ks},**{f'sd_{k}':sds[k] for k in ks}))
T=pd.DataFrame(rows); check_columns(T,'R246'); T.to_csv(pathlib.Path(__file__).parent/'results'/'axis.csv',index=False)
print(f"\n{'结局':<12}{'S(位置)':>12}{'Cres(内容)':>13}{'rho_i(何时)':>13}")
for _,r in T.iterrows():
    st=lambda k:'★' if abs(r[k])>2*r[f'sd_{k}'] else ' '
    print(f"{r.outcome:<12}{r[ks[0]]:>+11.4f}{st(ks[0])}{r[ks[1]]:>+12.4f}{st(ks[1])}{r[ks[2]]:>+12.4f}{st(ks[2])}")
same=[k for k in ks if np.sign(T.iloc[0][k])==np.sign(T.iloc[1][k])]
print(f"\n两道题在三个维度上**同号**的:{len(same)}/3  {same}")

g=Gate('animated 与 written 是不是一条轴的两端')
g.asserted('负对照:一道题与自己 = +1(管道自检)',abs(np.corrcoef(w_[mm],w_[mm])[0,1]-1)<1e-9,'+1.0000')
g.asserted('正对照刻度:确实二选一的一对读到负',r_bw<-0.15,f"corr(biomale, written) = {r_bw:+.4f}")
g.resolvable('corr(animated, written)',r_aw,sd_aw)
g.asserted('注册的 kill:corr >= 0 -> 「媒介偏好」这个名字错了',r_aw>=0,f"{r_aw:+.4f}")
g.asserted('⚠ 两道都是正向李克特,不该期待强负相关 —— 本轮判**符号**不判"够不够负"',True,
           '所以 kill 的判据写的是 `>= 0`,不是 `> −0.15`')
print(g)
print(f"\n  => {'TWO_THING —— 不是一条轴的两端,`#188` 起的名要重估' if r_aw>=0 else 'ONE_AXIS'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
