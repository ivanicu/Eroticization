import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A82 R321 -- `animated` · `written` · `allrollidentity` 是不是同一件事

`#275b`:六坐标最擅长的前几格里,这三个都是**形式/身份**类;
而 `#260b`/`#263b` 发现 `c3⊥D` 最强地指着 `animated` 与 `written`,**四次给不出名字**。
**换角度:不问「c3⊥D 是什么」,问「这三格是不是同一件事」。**

ESTIMAND        把六坐标从三个结局里各回归掉 -> 三条残差;判三条残差之间的两两相关。
⚠ 零不该是零      三者本来就都是情色 Likert 题,共享一个大共同因子(`#267a②`:λ1 占 7.3%)。
                **offset = 从 19 道 Likert 题里随机抽三道、同一条管道处理后的两两相关分布。**
KILL            **若三条残差的相关明显高于随机三元组的分布(> 其 95 分位)->
                六坐标之外还有一个专门管「形式」的人层量,它是第七维的候选;
                若落在分布内 -> 这三格各自独立,`#275b` 的「形式」只是我的归类。**
⚠ 与 `#266a` 不矛盾  那一轮判的是**全部 29 个结局**的残差谱;
                **一个只在三格上出现的因子会被全局谱淹没。**
POSITIVE CTRL   两端:① 三道**已知近重复**的题(两道 "existing as biological …" + 一道
                "masturbating alone as …")必须落在分布高端;② 三道随机题必须落在分布内。
NEGATIVE CTRL   跨人置换其中一条残差。
IMPOSSIBLE      「形式」是我给这三格起的名字;本轮判的是**它们是否共享残差**,不是这个名字对不对。
"""
import numpy as np, pandas as pd, warnings, hashlib, itertools
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
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S=np.where(ok,pos/np.maximum(cov,1),np.nan)
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
rg=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    o=rg.permutation(M.shape[1]); k=M.shape[1]//2
    A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
def prof(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
Ra,Rb=prof(A),prof(B)
st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
for i in np.flatnonzero(ok):
    mm=G[:,i]
    if mm.sum()<8: continue
    x,y_=Ra[mm,i],Rb[mm,i]
    if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
def z(v):
    m=np.isfinite(v); w=np.full(NN,np.nan); w[m]=(v[m]-v[m].mean())/v[m].std(); return w
def hs(R,cols):
    sub=R[cols]; F2=np.isfinite(sub)
    return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
C=np.full((NB,NB),np.nan)
for i in range(NB):
    for j in range(NB):
        mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&ok
        if mm.sum()>300: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); cs=[]
for k in range(3):
    num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
    cs.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
Q=[S,D]+cs+[st]
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
base=np.ones(NN,bool)&ok
for q_ in Q: base&=np.isfinite(q_)
X6=np.column_stack([np.ones(int(base.sum()))]+[(q_[base]-q_[base].mean())/q_[base].std() for q_ in Q])
def resid(col):
    y=pd.to_numeric(d[col],errors='coerce').values.astype(float)[base]
    f=np.isfinite(y); out=np.full(int(base.sum()),np.nan)
    if f.sum()<400: return out
    b=np.linalg.lstsq(X6[f],(y[f]-y[f].mean())/y[f].std(),rcond=None)[0]
    out[f]=(y[f]-y[f].mean())/y[f].std()-X6[f]@b; return out
RES={c:resid(c) for c in lik}
def pw(cols):
    vs=[]
    for a,b in itertools.combinations(cols,2):
        m=np.isfinite(RES[a])&np.isfinite(RES[b])
        if m.sum()>400: vs.append(float(np.corrcoef(RES[a][m],RES[b][m])[0,1]))
    return float(np.mean(vs)) if vs else np.nan
def find(sub): return next(c for c in lik if sub in str(c))
FORM=['animated','written','allrollidentity']
obs=pw(FORM)
rng=np.random.default_rng(20260804)
NULLS=[pw(list(rng.choice(lik,3,replace=False))) for _ in range(400)]
NULLS=[v for v in NULLS if np.isfinite(v)]
q95=float(np.quantile(NULLS,0.95)); q50=float(np.quantile(NULLS,0.50))
print(f"n = {int(base.sum()):,};19 道 Likert;随机三元组 {len(NULLS)} 个")
print(f"\n**三个「形式」结局的残差两两相关均值 = {obs:+.4f}**")
print(f"随机三元组分布:中位 {q50:+.4f} · 95 分位 **{q95:+.4f}** · 最大 {max(NULLS):+.4f}")
print(f"  -> 观测在分布的第 **{100*np.mean([v<obs for v in NULLS]):.1f}** 百分位")
DUP=[find('as a biological *female*'),find('as a biological *male*'),
     find('masturbating alone as a biological female')]
p1=pw(DUP); RND=list(rng.choice([c for c in lik if c not in FORM+DUP],3,replace=False)); p2=pw(RND)
print(f"\n正对照两端:① 已知近重复三道 **{p1:+.4f}**(第 {100*np.mean([v<p1 for v in NULLS]):.1f} 百分位)"
      f" · ② 随机三道 {p2:+.4f}(第 {100*np.mean([v<p2 for v in NULLS]):.1f} 百分位)")
m=np.isfinite(RES['animated'])&np.isfinite(RES['written'])
nul=[float(np.corrcoef(rng.permutation(RES['animated'][m]),RES['written'][m])[0,1]) for _ in range(30)]
print(f"负对照(置换 `animated` 的残差):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
prs={f"{a[:14]}×{b[:14]}":float(np.corrcoef(RES[a][np.isfinite(RES[a])&np.isfinite(RES[b])],
     RES[b][np.isfinite(RES[a])&np.isfinite(RES[b])])[0,1]) for a,b in itertools.combinations(FORM,2)}
print(f"  三对分别:"+' · '.join(f"{k} {v:+.4f}" for k,v in prs.items()))
# ⚠ 三元组均值不是这个问题的正确统计量:`allrollidentity` 与两者都无关,把均值拖下去了。
#    真正的对象是**一对**:`animated × written`。所以补一个**成对**的零 ——
#    19 道题全部 171 个配对的残差相关分布。
PAIRS={}
for a,b in itertools.combinations(lik,2):
    m=np.isfinite(RES[a])&np.isfinite(RES[b])
    if m.sum()>400: PAIRS[(a,b)]=float(np.corrcoef(RES[a][m],RES[b][m])[0,1])
pv=np.array(list(PAIRS.values()))
aw=prs[[k for k in prs if k.startswith('animated×written')][0]]
print(f"\n⚠ 成对的零(全部 {len(pv)} 个配对的残差相关):"
      f"中位 {np.median(pv):+.4f} · 95 分位 **{np.quantile(pv,0.95):+.4f}** · 最大 {pv.max():+.4f}")
top=sorted(PAIRS.items(),key=lambda kv:-kv[1])[:3]
print(f"   最高的三对:"+' · '.join(f"{a[:16]}×{b[:16]} {v:+.4f}" for (a,b),v in top))
print(f"   **`animated × written` = {aw:+.4f} -> 第 {100*np.mean(pv<aw):.1f} 百分位**")

T=pd.DataFrame([dict(pair=k,r=v) for k,v in prs.items()]+
               [dict(pair='均值',r=obs),dict(pair='随机95分位',r=q95)])
check_columns(T,'R321'); T.to_csv(pathlib.Path(__file__).parent/'results'/'form_residuals.csv',index=False)

g=Gate('三个「形式」结局是不是同一件事')
g.asserted('正对照两端:已知近重复必须落在高端,随机三道必须落在分布内',
           p1>q95 and abs(p2-q50)<0.05, f"① {p1:+.4f}(>95分位 {q95:+.4f})· ② {p2:+.4f}(中位 {q50:+.4f})")
g.negative_control('置换 `animated` 的残差',abs(float(np.mean(nul))),abs(obs),
                   null_spread=float(np.std(nul)),null_kind='跨人置换一条残差 —— 只打掉配对')
g.offset_control('★ 三个「形式」结局 vs 随机三元组的 95 分位',obs,q95,
                 float(np.std(NULLS)),
                 null_kind='19 道 Likert 里随机三道、同一条管道处理后的两两相关 —— '
                           '不是零假设,是「任意三道情色题本来就有多像」')
g.asserted('★ 注册的 kill(三元组):三条残差的相关明显高于随机三元组 -> 一个专管「形式」的人层量',
           obs>q95, f"观测 {obs:+.4f} vs 95 分位 {q95:+.4f};第 {100*np.mean([v<obs for v in NULLS]):.1f} 百分位")
g.asserted('⚠ 三元组均值不是正确的统计量:`allrollidentity` 与两者都无关,把均值拖下去了',
           abs(prs[[k for k in prs if 'allroll' in k][0]])<0.05,
           ' · '.join(f"{k} {v:+.4f}" for k,v in prs.items()))
g.offset_control('★★ 真正的对象:`animated × written` vs 全部 171 个配对的 95 分位',
                 aw,float(np.quantile(pv,0.95)),float(np.std(pv)),
                 null_kind='19 道 Likert 全部配对的残差相关 —— 不是零假设,'
                           '是「任意两道情色题去掉六坐标之后本来就有多像」')
g.asserted('★★ 修正后的 kill:`animated × written` 高于全部配对的 95 分位 -> 一个六坐标之外的双题因子',
           aw>float(np.quantile(pv,0.95)),
           f"{aw:+.4f} vs 95 分位 {np.quantile(pv,0.95):+.4f};第 {100*np.mean(pv<aw):.1f} 百分位;"
           f"而近重复正对照是 {p1:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
