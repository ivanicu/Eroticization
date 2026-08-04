import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A57 R275 -- 宽度类型预测的结局是不是不一样的(把第四件事接回 A/B/C)

`#229a`/`#229b` 给了一个 A/B/C 之外的结构:宽度有形状,至少三维,第一个族是体液。
**但它还没跟主问题接上** —— 一个漂亮的聚类,若三个类型预测同一批东西,那它只是幅度的分解。

WORLDS          ① **幅度的分解**:三个成分是同一件事的三次带噪声读数 ->
                   三条结局剖面**高度共线**(接近"同一件事"的上限)
                ② **真正不同的类型**:剖面明显不共线,而且 2/3 号成分**自己**能打中结局 ->
                   「被什么打开」与「因此感到什么」是绑定的,
                   **那是 C(递归)第一个真正的立足点**
ESTIMAND        前 3 个宽度成分的人层得分 × 31 个结局 = 3 条剖面;判剖面间的两两相关。
KILL            **若两两相关接近「同一件事」的上限 -> 幅度的分解,`#229` 只是结构不是内容;
                若明显低于上限、且成分 2 或 3 自己有越阈结局 -> 真正的类型。**
⚠ 零应该是零吗?        **不应该。** 三个成分共享同一批人、同一批块、同一条管道,
                即使世界②为真,剖面也不会独立。所以判据是 **offset_control**:
                「若它们是同一件事的三次带噪声读数,剖面相关该是多少」——
                用**成分 1 的得分 + 噪声**,噪声量校准到各成分**自己的分半信度**。
POSITIVE CTRL   已知不同的三个真变量(age / openness / neuroticism)-> 剖面相关必须低;
                同一变量的三次带噪声复制 -> 剖面相关必须高。**两端都要过,才证明这把尺子能分辨。**
IMPOSSIBLE      成分 2/3 的方向在旋转下不唯一(`#229` 已记),
                所以只能判**「这三条剖面加起来是不是一件事」**,不能判「第 2 类型是什么」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
def halves(seed):
    rg=np.random.default_rng(seed); H=np.full((2,NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB):
        o=rg.permutation(M.shape[1]); k=M.shape[1]//2
        H[0,b,ppl]=M[:,o[:k]].mean(1); H[1,b,ppl]=M[:,o[k:2*k]].mean(1)
    return H
def profile(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
H=halves(500); Ra,Rb=profile(H[0]),profile(H[1])
C=np.full((NB,NB),np.nan)
for i in range(NB):
    for j in range(NB):
        m=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
        if m.sum()>300: C[i,j]=np.corrcoef(Ra[i][m],Rb[j][m])[0,1]
C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2
V=np.linalg.eigh(C)[1][:,::-1]
def sc(R,k):
    F=np.isfinite(R); Z=np.where(F,R,0.0); num=(V[:,k][:,None]*Z).sum(0); den=(F*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
SA=np.array([sc(Ra,k) for k in range(3)]); SB=np.array([sc(Rb,k) for k in range(3)])
REL=[]
for k in range(3):
    m=np.isfinite(SA[k])&np.isfinite(SB[k])&ok; REL.append(float(np.corrcoef(SA[k][m],SB[k][m])[0,1]))
S=np.array([(SA[k]+SB[k])/2 for k in range(3)])
print(f"块 {NB};n = {int(ok.sum()):,};三个成分的人层得分分半信度 "
      + ' · '.join(f"c{k+1} {REL[k]:+.4f}" for k in range(3)))

lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EXTRA={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EXTRA.items()]
rng=np.random.default_rng(777)
def outprofile(x,tag=None):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl)
    thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),axis=0),0.95))
    r=np.array(r); nhit=int(np.nansum(np.abs(r)>thr))
    if tag:
        top=sorted([(OUT[i][0],r[i]) for i in range(len(r)) if np.isfinite(r[i]) and abs(r[i])>thr],
                   key=lambda t:-abs(t[1]))[:4]
        print(f"   {tag:<10} {nhit:>2}/{len(OUT)}(阈值 {thr:.4f}) "+' · '.join(f"{n[:20]} {v:+.3f}" for n,v in top))
    return r,nhit,thr
def pairs(P):
    out=[]
    for i in range(3):
        for j in range(i+1,3):
            m=np.isfinite(P[i])&np.isfinite(P[j]); out.append(float(np.corrcoef(P[i][m],P[j][m])[0,1]))
    return out
print("\n三个宽度成分各自的结局剖面:")
PR=[outprofile(S[k],f"成分 {k+1}") for k in range(3)]
obs=pairs([p[0] for p in PR])
print(f"   -> 剖面两两相关(观测):"+' · '.join(f"{v:+.4f}" for v in obs)+f"  |mean| = {np.mean(np.abs(obs)):.4f}")

# offset:若它们是同一件事的三次带噪声读数
def noisy(x,rel,seed):
    m=np.isfinite(x); z=np.full(NN,np.nan); v=x[m]; v=(v-v.mean())/v.std()
    n=np.random.default_rng(seed).standard_normal(m.sum())
    z[m]=np.sqrt(max(rel,1e-3))*v+np.sqrt(max(1-rel,0))*n; return z
OFF=[]
for t in range(3):
    P=[outprofile(noisy(S[0],REL[k],1000+10*t+k))[0] for k in range(3)]; OFF.append(pairs(P))
off=np.mean(np.abs(np.array(OFF)),axis=0)
print(f"   -> offset(同一件事的三次带噪声读数,噪声按各成分自己的信度校准):"
      +' · '.join(f"{v:.4f}" for v in off)+f"  |mean| = {off.mean():.4f}")

# 正对照:已知不同 vs 已知相同
known=[EXTRA['age'].values.astype(float),EXTRA['openness'].values.astype(float),
       EXTRA['neuroticism'].values.astype(float)]
kd=np.mean(np.abs(pairs([outprofile(k_)[0] for k_ in known])))
ks=np.mean(np.abs(pairs([noisy(known[0],0.6,2000+i) for i in range(3)] and
                        [outprofile(noisy(known[0],0.6,2000+i))[0] for i in range(3)])))
print(f"\n正对照:已知不同的三个真变量(age/openness/neuroticism)剖面 |r| = **{kd:.4f}** · "
      f"同一变量的三次带噪声复制 = **{ks:.4f}**")

# ---------- 它是不是位置分换了张皮?(#179 的位置分是这个项目的主claim) ----------
posv=np.zeros(NN); cntv=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); posv[ppl]+=(M@rr)/np.maximum(M.sum(1),1); cntv[ppl]+=1
Spos=np.where(ok,posv/np.maximum(cntv,1),np.nan)
Kn=np.zeros(NN)
for M,ppl in MB: Kn[ppl]+=M.sum(1)
Kn=np.where(ok,Kn,np.nan)
SH=[nm for nm,_ in OUT if nm.lower().startswith('"i am ashamed')]
print(f"\n它是不是位置分换了张皮?")
for k in range(3):
    m=np.isfinite(S[k])&np.isfinite(Spos)
    print(f"  成分 {k+1}: corr(位置分 S) = {np.corrcoef(S[k][m],Spos[m])[0,1]:+.4f} · "
          f"corr(勾选数) = {np.corrcoef(S[k][m],Kn[m])[0,1]:+.4f}")
if SH:
    y=dict(OUT)[SH[0]]
    def cr(x):
        m=np.isfinite(x)&np.isfinite(y)&ok; return float(np.corrcoef(x[m],y[m])[0,1])
    print(f"  羞耻题「{SH[0][:56]}」:位置分 S {cr(Spos):+.4f} · "
          + ' · '.join(f"c{k+1} {cr(S[k]):+.4f}" for k in range(3)))
    import numpy.linalg as _la
    m=np.isfinite(S[2])&np.isfinite(Spos)&np.isfinite(y)&ok
    X=np.column_stack([np.ones(m.sum()),(S[2][m]-S[2][m].mean())/S[2][m].std(),
                       (Spos[m]-Spos[m].mean())/Spos[m].std()])
    b=_la.lstsq(X,(y[m]-y[m].mean())/y[m].std(),rcond=None)[0]
    print(f"  同时放进去(标准化 beta):成分 3 = **{b[1]:+.4f}** · 位置分 S = **{b[2]:+.4f}**(n = {m.sum():,})")

T=pd.DataFrame(dict(comp=[1,2,3],rel=REL,n_hit=[p[1] for p in PR],thr=[p[2] for p in PR]))
T2=pd.DataFrame(dict(pair=['1-2','1-3','2-3'],observed=obs,offset_same_thing=off))
check_columns(T,'R275'); check_columns(T2,'R275')
T.to_csv(pathlib.Path(__file__).parent/'results'/'components.csv',index=False)
T2.to_csv(pathlib.Path(__file__).parent/'results'/'profile_pairs.csv',index=False)

g=Gate('宽度类型预测的结局是不是不一样的')
g.asserted('正对照两端:已知不同必须低、已知相同必须高',
           kd<ks-0.15, f"已知不同 {kd:.4f} vs 已知相同 {ks:.4f}")
g.offset_control('★ 剖面两两相关 vs「同一件事的三次带噪声读数」',
                 float(np.mean(np.abs(obs))),float(off.mean()),
                 float(np.std(np.abs(np.array(OFF)).mean(axis=1))),
                 null_kind='成分 1 的得分 + 校准噪声 —— 不是零假设,是「若三个成分是同一件事,剖面该有多像」')
g.asserted('★ 注册的 kill:剖面明显低于上限,且成分 2 或 3 自己有越阈结局 -> 真正的类型',
           np.mean(np.abs(obs))<off.mean()-0.15 and max(PR[1][1],PR[2][1])>0,
           f"观测 {np.mean(np.abs(obs)):.4f} vs 上限 {off.mean():.4f};"
           f"越阈结局 c1 {PR[0][1]} · c2 {PR[1][1]} · c3 {PR[2][1]}")
print(g)
print(f"\nsha1 {hashlib.sha1(T2.to_csv(index=False).encode()).hexdigest()[:12]}")
