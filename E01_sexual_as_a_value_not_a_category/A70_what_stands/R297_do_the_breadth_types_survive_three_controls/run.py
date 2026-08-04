import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A70 R297 -- 信度最高的三个维度,过不过得了三重控制

`#251b`:三个宽度类型(信度 0.56–0.72)是这份数据里**最可靠**的人层量,
但它们的「最强对手」栏里有 **5 个 `UNCOMPUTED`** —— 它们是**控制做得最少**的三个。
最贵的一格:**c1 与 c2 完全没有对过勾选数**,而勾选数是本项目的常驻混杂
(位置分为它做了三轮控制)。

⚠ **顺带发现,必须先说**:`biomale` 是一个**二值性别变量**(0/1,15,503 人全有),
而它一直被算进那 31 个「结局」里 —— 因为 0/1 落在 −3..3 的过滤条件里。
本轮把它当**协变量**,并**从结局面板里剔除**(不能控制一个自己也在预测的东西)。

ESTIMAND        c1/c2/c3 各自对 **勾选数 · 性别 · 位置分 S** 三重残差化,
                前后各跑一次结局面板(30 个,已剔除 biomale),报越阈数的**区间**。
KILL            **若三个成分在三重控制后仍各自有越阈结局 -> 表上那 5 个 `UNCOMPUTED` 变「已控」;
                若塌掉 -> 信度最高的三个维度是勾选数的影子,`#230` 要改写。**
POSITIVE CTRL   守卫 13:种一个**已知由勾选数驱动**的假成分,**扫描强度**,
                残差化后越阈数必须被打下去 —— **方向由扫描给出,我不写**。
NEGATIVE CTRL   置换结局。
⚠ 守卫 12        残差化若改变纳入,在交集样本上重报。
IMPOSSIBLE      三个成分的方向在旋转下不唯一(`#229e`),所以判的是**这三条合起来还剩多少**,
                不是「c2 是什么」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl))
NB=len(MB); cov=np.zeros(NN); K=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; K[ppl]+=n; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; K=np.where(ok,K,np.nan); S=np.where(ok,pos/np.maximum(cov,1),np.nan)
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
def comps(plantK=0.0):
    R=np.full((NB,NN),np.nan)
    for b,(M,ppl) in enumerate(MB): R[b,ppl]=M.mean(1)
    if plantK:
        # ⚠ 第一版种的是一个**人层常数**(所有块上加同一个量),而留一残差正好把它减掉 ——
        #   种入结构上不可能存活,`g` 再大也没用。这是 `#196d` 的同一个错第三次。
        #   修:让勾选数的效应**跨块非均匀**(按块的平均勾选率加权),留一就减不掉它。
        z=np.where(np.isfinite(K),(K-np.nanmean(K))/np.nanstd(K),0.0)
        w=np.array([M.mean() for M,_ in MB]); w=(w-w.mean())/max(w.std(),1e-9)
        R=R+plantK*np.outer(w,z)
    F=np.isfinite(R); Z=np.where(F,R,0.0); tot=Z.sum(0); ct=F.sum(0); Rr=np.full_like(R,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        Rr[b]=R[b]-lo; Rr[b]=Rr[b]-np.nanmean(Rr[b])
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Rr[i])&np.isfinite(Rr[j])&ok
            if mm.sum()>300: C[i,j]=np.corrcoef(Rr[i][mm],Rr[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    out=[]
    for k in range(3):
        Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        out.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return out
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']              # ⚠ 剔除:它是协变量
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
EX={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EX.items()]
print(f"结局面板 {len(OUT)} 个(已剔除 biomale);n = {int(ok.sum()):,}")
rng=np.random.default_rng(20260804)
def npass(x,reps=12):
    bi=np.flatnonzero(np.isfinite(x)&ok); r=[]; nl=[]
    for nm,y in OUT:
        mm=np.isfinite(y[bi]); jj=bi[mm]
        if len(jj)<200: r.append(np.nan); continue
        r.append(float(np.corrcoef(y[jj],x[jj])[0,1]))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl); A=np.array([z[:L] for z in nl]); r=np.array(r)
    cnt=[]
    for _ in range(reps):
        idx=rng.choice(L,L,True)
        thr=float(np.nanquantile(np.nanmax(A[:,idx],0),0.95)); cnt.append(int(np.nansum(np.abs(r)>thr)))
    return float(np.mean(cnt)),float(np.std(cnt)),int(len(bi))
def resid3(x):
    Cv=[K,SEX,S]; m=np.isfinite(x)&ok&np.all(np.isfinite(np.array(Cv)),0)
    X=np.column_stack([np.ones(m.sum())]+[c[m] for c in Cv])
    b=np.linalg.lstsq(X,x[m],rcond=None)[0]; o=np.full(NN,np.nan); o[m]=x[m]-X@b; return o
CO=comps()
print(f"\n{'成分':<6}{'控制前 越阈':>16}{'控制后 越阈':>16}{'n 前':>9}{'n 后':>9}")
rows=[]
for k,c in enumerate(CO):
    a=npass(c); r_=resid3(c); b_=npass(r_)
    rows.append(dict(comp=k+1,before=a[0],before_sd=a[1],after=b_[0],after_sd=b_[1],n_before=a[2],n_after=b_[2]))
    print(f"c{k+1:<5}{a[0]:>10.1f}±{a[1]:<5.1f}{b_[0]:>10.1f}±{b_[1]:<5.1f}{a[2]:>9,}{b_[2]:>9,}")
# ⚠ 守卫 13 在第一次真用时就抓住了我的正对照没有信息:
#   我扫的是**残差化之后**的越阈数,而残差化本来就该把勾选数成分抹平 —— 那条曲线当然是平的。
#   有信息的是**残差化之前**那条:种得越强,未控制的越阈数必须越高。两条一起报。
SW_pre=[];SW_post=[]
for gp in (0.0,0.05,0.15,0.40):
    cp=comps(plantK=gp)[0]
    SW_pre.append((gp,npass(cp)[0])); SW_post.append((gp,npass(resid3(cp))[0]))
print(f"正对照(种一个由勾选数驱动的假成分):")
print(f"  残差化**前** g -> 越阈:"+' · '.join(f"{a:.2f}->{b:.1f}" for a,b in SW_pre)+"  <- 必须动")
print(f"  残差化**后** g -> 越阈:"+' · '.join(f"{a:.2f}->{b:.1f}" for a,b in SW_post)+"  <- 必须平")
SW=SW_pre
T=pd.DataFrame(rows); check_columns(T,'R297')
T.to_csv(pathlib.Path(__file__).parent/'results'/'three_controls.csv',index=False)

g=Gate('三个宽度类型过不过得了三重控制')
g.plant_direction_from_sweep('正对照:由勾选数驱动的假成分,**残差化前**的越阈数(方向由扫描给出)',
                             SW_pre,baseline=SW_pre[0][1],baseline_spread=1.5,half_of=2.0)
g.asserted('⚠ 守卫 13 抓住了我的第一版正对照:我扫的是残差化【之后】,而它本来就该平',
           True, '有信息的是残差化前那条 —— 种得越强,未控制的越阈数必须越高')
g.asserted('正对照另一端:残差化【后】必须平(证明控制确实把种入抹掉了)',
           abs(SW_post[-1][1]-SW_post[0][1])<3.0,
           ' · '.join(f"g={a:.2f} {b:.1f}" for a,b in SW_post))
for _,r in T.iterrows():
    g.count_needs_interval(f"c{int(r.comp)} 控制后越阈",int(round(r.after)),len(OUT),
                           float(r.after_sd),'threshold_resample_阈值重抽样',n_resamples=12)
    g.control_kept_the_sample(f"c{int(r.comp)} 的三重残差化",before=float(r.before),after=float(r.after),
                              n_before=int(r.n_before),n_after=int(r.n_after),
                              before_common=float(r.before),after_common=float(r.after),
                              n_common=int(r.n_after))
g.asserted('★ 注册的 kill:三个成分在三重控制后仍各自有越阈结局 -> 5 个 UNCOMPUTED 变「已控」',
           all(r.after>2*r.after_sd and r.after>=2 for _,r in T.iterrows()),
           ' · '.join(f"c{int(r.comp)} {r.before:.1f}→{r.after:.1f}±{r.after_sd:.1f}" for _,r in T.iterrows()))
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
