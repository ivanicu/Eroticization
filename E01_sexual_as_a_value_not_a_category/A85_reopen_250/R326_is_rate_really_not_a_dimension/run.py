import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A85 R326 -- `#250a` 那一步站不站得住

`#280b`:`#250a` 用「剖面相关 −0.8499 贴着上限 0.8779」得出
**「`rate` 不是一个独立的人格维度,它就是 Δ 反过来写」**,并据此关掉整条 A67–A69 线。
**而它从没报过 `corr(rate, rho_i)` 的分数层值** —— 那正是守卫 15 挡的方向。

ESTIMAND        ① `corr(rate, rho_i)` 的**分数层**值(全仪器口径);
                ② 用守卫 15 判 `#250a` 的那一步;
                ③ 补 `#250a` 当时没做的检验:把 `rho_i` 从 `rate` 里回归掉,
                   残差 `rate⊥rho` 还打中几个结局(`#305` 给 `c3⊥D` 做过的同一件事)。
KILL            **若分数层相关也高(|r| ≥ 0.3)-> `#250a` 成立,守卫 15 放行,那条线关得对;
                若分数层低而剖面高 -> `#250a` 是被禁方向上的结论,必须重开,
                而 `rate` 可能真的是第七个维度。**
POSITIVE CTRL   守卫 15 的四端已在 `#325` 验过;本轮另加:
                把 `rho_i` 从**它自己的带噪声复制**里回归掉,残差必须什么都不剩(`#305` 同款)。
NEGATIVE CTRL   跨人置换(**只在有限值内**,`#264b`/`#278b`)。
IMPOSSIBLE      `rate` 与 `rho_i` 算自**同一批起始年龄**(`#248d` 的那条边仍开着)。
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
RHO=rho_of(V0); NCAT=np.isfinite(V0).sum(1)
Vs=np.sort(np.where(np.isfinite(V0),V0,np.inf),axis=1)
A1=Vs[:,0]; ACL=np.nanmax(np.where(np.isfinite(V0),V0,np.nan),axis=1)
RATE=np.where(np.isfinite(RHO),(NCAT-1)/np.maximum(ACL-A1,0.5),np.nan)
m0=np.isfinite(RHO)&np.isfinite(RATE)
r_score=float(np.corrcoef(RATE[m0],RHO[m0])[0,1])
print(f"n = {int(m0.sum()):,}")
print(f"**`corr(rate, rho_i)` 分数层 = {r_score:+.4f}**  ——  `#250a` 报的剖面层是 **−0.8499**")
def resid(x,c):
    m=np.isfinite(x)&np.isfinite(c); X=np.column_stack([np.ones(m.sum()),c[m]])
    o=np.full(N,np.nan); o[m]=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]; return o
RES=resid(RATE,RHO)
lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
EX={'age':df['age'].map(AGE),'openness':pd.to_numeric(df['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(df['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(df['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(df['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(df['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(df['powerlessnessvariable'],errors='coerce'),
 '关系风格':df['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':df['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '成长期性开放度':df['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,df[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
rngB=np.random.default_rng(20260804)
def panel(x,reps=12,tag=None):
    bi=np.flatnonzero(np.isfinite(x)); r=[]; nl=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rngB.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(v) for v in nl); Am=np.array([v[:L] for v in nl]); r=np.array(r); c=[]
    for _ in range(reps):
        i2=rngB.choice(L,L,True)
        c.append(int(np.nansum(np.abs(r)>float(np.nanquantile(np.nanmax(Am[:,i2],0),0.95)))))
    thr=float(np.nanquantile(np.nanmax(Am,0),0.95))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
                   key=lambda t:-abs(t[1]))[:5]
        print(f"   {tag:<16} {np.mean(c):.1f}±{np.std(c):.1f}/{len(OUT)} "
              +' · '.join(f"{n[:20]} {v:+.3f}" for n,v in top))
    return float(np.mean(c)),float(np.std(c))
print(f"\n③ 补 `#250a` 当时没做的检验:")
n_rate,_=panel(RATE,tag='rate 原始')
n_res,sd_res=panel(RES,tag='rate⊥rho')
def noisy(x,r_,seed):
    m=np.isfinite(x); zz=np.full(N,np.nan); v=(x[m]-np.nanmean(x))/np.nanstd(x)
    zz[m]=np.sqrt(r_)*v+np.sqrt(1-r_)*np.random.default_rng(seed).standard_normal(int(m.sum())); return zz
n_pc,_=panel(resid(noisy(RHO,0.7,41),RHO))
print(f"   正对照(把 rho 从它自己的带噪声复制里回归掉):残差越阈 **{n_pc:.1f}**(必须 ≈0)")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2)); z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
n_nul=np.mean([panel(perm_finite(RES,90+i))[0] for i in range(3)])
print(f"   负对照(置换 `rate⊥rho`,只在有限值内):{n_nul:.1f}")
T=pd.DataFrame([dict(quantity='rate',n_hit=n_rate),dict(quantity='rate⊥rho',n_hit=n_res,sd=sd_res),
                dict(quantity='corr_score',n_hit=r_score)])
check_columns(T,'R326'); T.to_csv(pathlib.Path(__file__).parent/'results'/'reopen_250.csv',index=False)

g=Gate('`#250a` 那一步站不站得住')
g.profile_similarity_is_not_identity('★ 守卫 15 判 `#250a` 的那一步(剖面 −0.8499)',
                                     -0.8499, r_score)
g.asserted('正对照:把 rho 从它自己的带噪声复制里回归掉,残差必须什么都不剩',
           n_pc<=2.0, f"残差越阈 {n_pc:.1f}")
g.negative_control('置换 `rate⊥rho`',abs(n_nul),abs(n_res),null_spread=None,
                   null_kind='跨人置换(只在有限值内)—— 只打掉配对')
g.count_needs_interval('`rate⊥rho` 的越阈计数',int(round(n_res)),len(OUT),sd_res,
                       'threshold_resample_阈值重抽样',n_resamples=12,seed_spread=sd_res)
g.asserted('★ 注册的 kill:分数层 |r| ≥ 0.3 -> `#250a` 成立;分数层低而剖面高 -> 必须重开',
           abs(r_score)>=0.3,
           f"分数层 {r_score:+.4f} vs 剖面层 −0.8499;`rate⊥rho` 越阈 {n_res:.1f}±{sd_res:.1f}/{len(OUT)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
