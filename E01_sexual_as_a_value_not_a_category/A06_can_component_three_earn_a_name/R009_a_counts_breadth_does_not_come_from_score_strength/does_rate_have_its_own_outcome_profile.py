import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A69 R295 -- 「积累得快」是一个独立的人格维度,还是 Δ 的另一种写法

`#247a`:Δ 沿积累速率从 +0.088 走到 −0.361;`#248b`:跨半设计下保住 56%。
**但整条线到现在为止,Δ 只跟「什么时候」绑在一起,从没跟「这个人是谁」绑过。**

WORLDS          ① **独立维度**:`rate` 有自己的结局剖面,与 `Δ` 的不重合 ->
                   继位置分(`#179`)、宽度类型(`#229`)之后的第三个人层维度
                ② **同一件事的两种写法**:两条剖面几乎相同 -> 这条线到此为止
ESTIMAND        `rate` 与 `rho_i` 各跑 31 个结局(逐个重算最大统计量阈值),判两条剖面的相关。
⚠ 零应该是零吗?     **不应该**(`#230d` 同款):两者算自同一批起始年龄,共享噪声。
                判据是 **offset**:「若是同一件事的两次带噪声读数,剖面该有多像」,
                噪声按**各自的分半信度**校准。
KILL            **若剖面明显低于上限 且 `rate` 自己有越阈结局 -> 世界①;
                若接近上限 -> 世界②。**
POSITIVE CTRL   两端:已知不同(age vs openness)剖面必须低;
                同一变量的两次带噪声复制必须高。
NEGATIVE CTRL   置换结局。
IMPOSSIBLE      `rate` 与 `rho_i` 算自同一批起始年龄 —— 本轮判的是**剖面是否可互相替代**,
                不是它们在人身上是否同源。跨仪器版本这份数据没有(块仪器没有时间)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
def rho_sub(Vm,cols,need=8):
    Vc=np.full_like(Vm,np.nan); Vc[:,cols]=Vm[:,cols]
    D=np.where(np.isfinite(Vc),Vc,np.nan)
    for _ in range(250):
        a=np.nanmean(D,0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
    W=np.isfinite(D); Z=np.where(W,D,0.0); k=W.sum(1)
    rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    out=np.full(N,np.nan); ok=(k>=need)&(den>1e-12); out[ok]=num[ok]/den[ok]; return out
def rate_sub(Vm,cols,need=3):
    Vc=np.where(np.isfinite(Vm[:,cols]),Vm[:,cols],np.nan)
    k=np.isfinite(Vc).sum(1); lo=np.nanmin(Vc,1); hi=np.nanmax(Vc,1)
    return np.where(k>=need,(k-1)/np.maximum(hi-lo,0.5),np.nan)
ALL=np.arange(Mc)
RHO=rho_sub(V0,ALL); RATE=rate_sub(V0,ALL,need=8)
rng=np.random.default_rng(20260804)
def rel(fn,need):
    vs=[]
    for s in range(4):
        p=np.random.default_rng(60+s).permutation(Mc); h=Mc//2
        a,b=fn(V0,p[:h],need),fn(V0,p[h:],need); m=np.isfinite(a)&np.isfinite(b)
        r=float(np.corrcoef(a[m],b[m])[0,1]); vs.append(2*r/(1+r) if r<0.999 else np.nan)
    return float(np.nanmean(vs))
RRHO=rel(rho_sub,5); RRATE=rel(rate_sub,3)
m0=np.isfinite(RHO)&np.isfinite(RATE)
print(f"n = {int(m0.sum()):,};分半信度:rho_i **{RRHO:+.4f}** · rate **{RRATE:+.4f}**;"
      f"corr(rate, rho_i) = **{np.corrcoef(RATE[m0],RHO[m0])[0,1]:+.4f}**")

lik=[c for c in df.columns if df[c].dtype!=object and
     set(pd.Series(df[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EX={'age':df['age'].map(AGE),'openness':pd.to_numeric(df['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(df['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(df['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(df['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(df['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(df['powerlessnessvariable'],errors='coerce'),
 '关系风格':df['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':df['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':df['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':df['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,df[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
def prof(x,tag=None):
    bi=np.flatnonzero(np.isfinite(x)); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),0),0.95))
    r=np.array(r); h=int(np.nansum(np.abs(r)>thr))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
                   key=lambda t:-abs(t[1]))[:5]
        print(f"   {tag:<10} {h:>2}/{len(OUT)}(阈值 {thr:.4f}) "+' · '.join(f"{n[:20]} {v:+.3f}" for n,v in top))
    return r,h
print("\n两个量各自的结局剖面:")
pR,hR=prof(np.where(m0,RATE,np.nan),'rate'); pD,hD=prof(np.where(m0,RHO,np.nan),'rho_i')
mm=np.isfinite(pR)&np.isfinite(pD); obs=float(np.corrcoef(pR[mm],pD[mm])[0,1])
def noisy(x,r_,seed):
    m=np.isfinite(x); z=np.full(N,np.nan); v=(x[m]-x[m].mean())/x[m].std()
    z[m]=np.sqrt(max(r_,1e-3))*v+np.sqrt(max(1-r_,0))*np.random.default_rng(seed).standard_normal(m.sum())
    return z
OFF=[]
for t in range(4):
    a_,b_=prof(noisy(np.where(m0,RATE,np.nan),RRATE,5000+2*t))[0],prof(noisy(np.where(m0,RATE,np.nan),RRHO,5001+2*t))[0]
    q=np.isfinite(a_)&np.isfinite(b_); OFF.append(abs(float(np.corrcoef(a_[q],b_[q])[0,1])))
print(f"   -> 剖面相关 观测 **{obs:+.4f}** · offset(同一件事的两次带噪声读数)"
      f"**{np.mean(OFF):.4f} ± {np.std(OFF):.4f}**")
ka,kb=prof(EX['age'].values.astype(float))[0],prof(EX['openness'].values.astype(float))[0]
q=np.isfinite(ka)&np.isfinite(kb); kd=abs(float(np.corrcoef(ka[q],kb[q])[0,1]))
sa,sb=prof(noisy(EX['age'].values.astype(float),0.6,71))[0],prof(noisy(EX['age'].values.astype(float),0.6,72))[0]
q=np.isfinite(sa)&np.isfinite(sb); ks=abs(float(np.corrcoef(sa[q],sb[q])[0,1]))
print(f"正对照两端:已知不同(age vs openness)**{kd:.4f}** · 同一变量两次带噪声复制 **{ks:.4f}**")

T=pd.DataFrame([dict(quantity='rate',rel=RRATE,n_hit=hR),dict(quantity='rho_i',rel=RRHO,n_hit=hD)])
check_columns(T,'R295'); T.to_csv(pathlib.Path(__file__).parent/'results'/'rate_vs_rho_profiles.csv',index=False)

g=Gate('积累速率是不是一个独立的人格维度')
g.asserted('正对照两端:已知不同必须低、已知相同必须高',kd<ks-0.15,f"不同 {kd:.4f} vs 相同 {ks:.4f}")
g.offset_control('★ 剖面相关 vs「同一件事的两次带噪声读数」',abs(obs),float(np.mean(OFF)),
                 float(np.std(OFF)),
                 null_kind='rate 的得分 + 校准噪声 —— 不是零假设,是「若两者是同一件事,剖面该有多像」')
g.count_needs_interval('rate 的越阈计数',hR,len(OUT),1.2,'threshold_resample_阈值重抽样',n_resamples=20)
g.asserted('★ 注册的 kill:剖面明显低于上限 且 rate 自己有越阈结局 -> 独立维度',
           abs(obs)<np.mean(OFF)-0.15 and hR>0,
           f"观测 {obs:+.4f} vs 上限 {np.mean(OFF):.4f};rate 越阈 {hR}/{len(OUT)} · rho_i {hD}/{len(OUT)}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
