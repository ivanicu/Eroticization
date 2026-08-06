import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A76 R306 -- `c3⊥D` 是不是「对非真人的表征更容易被唤起」

⚠ **分组与预测已在跑之前冻结并提交**:见同目录 `PREREGISTRATION.md`(独立 commit)。
⚠ **污染声明**:`#260` 已让我看到 12 道里 4 道的值,所以**以 8 道保留集为准**。

KILL            **若两组的平均相关反号且差 > 2×展布(在 8 道保留集上)-> 名字挣到第一份证据;
                若两组相当 -> 第三次命名失败,而那本身值得记。**
POSITIVE CTRL   两端:① 只贴 A 组的合成分必须被分开;② 与两组等相关的必须不被分开。
NEGATIVE CTRL   跨人置换 `c3⊥D`。
IMPOSSIBLE      「需要真实在场他人」由我读题目文本判定,是**一个人的编码**;
                `#203c` 要求的独立编码者本会话不可得,如实登记。
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
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
ok=cov>=8
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
def build(seed):
    rg=np.random.default_rng(seed); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
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
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(Vv[:,2][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,2])[:,None]).sum(0)
    c3=np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
    m=np.isfinite(c3)&np.isfinite(D)
    X=np.column_stack([np.ones(m.sum()),D[m]]); o2=np.full(NN,np.nan)
    o2[m]=c3[m]-X@np.linalg.lstsq(X,c3[m],rcond=None)[0]
    return o2
RES=build(500)
cols=list(d.columns)
def find(sub): return next(c for c in cols if sub in str(c))
A_ALL=['animated','written',find('as a biological *female*'),find('as a biological *male*'),
       find('masturbating alone as a biological female'),find('masturbating alone as a biological male')]
B_ALL=[find("aren't aroused by it themselves"),find('two people of the opposite gender'),
       'knowwhatarousesyou',find('If my partner is aroused by something'),
       find('aroused by being dominant'),find('aroused by being submissive')]
SEEN={'animated','written',find('aroused by being dominant'),find('two people of the opposite gender')}
A_HO=[c for c in A_ALL if c not in SEEN]; B_HO=[c for c in B_ALL if c not in SEEN]
print(f"A 组(不需要真实在场他人){len(A_ALL)} 道 · B 组(需要){len(B_ALL)} 道")
print(f"⚠ 保留集(#260 未看过):A {len(A_HO)} 道 · B {len(B_HO)} 道")
rng=np.random.default_rng(20260804)
def grp(x,names):
    rs=[]
    for c in names:
        y=pd.to_numeric(d[c],errors='coerce').values.astype(float)
        m=np.isfinite(x)&np.isfinite(y)&ok
        if m.sum()>200: rs.append(float(np.corrcoef(x[m],y[m])[0,1]))
    return float(np.mean(rs)),rs
def contrast(x,Ag,Bg):
    a,_=grp(x,Ag); b,_=grp(x,Bg); return a-b,a,b
def boot(x,Ag,Bg,n=200):
    idx=np.flatnonzero(np.isfinite(x)&ok); out=[]
    for _ in range(n):
        s=rng.choice(idx,len(idx),True); xx=np.full(NN,np.nan); xx[idx]=x[s]
        out.append(contrast(xx,Ag,Bg)[0])
    return float(np.std(out))
for tag,Ag,Bg in (('全部 12 道',A_ALL,B_ALL),('★ 保留集 8 道',A_HO,B_HO)):
    g_,a_,b_=contrast(RES,Ag,Bg); sd=boot(RES,Ag,Bg)
    print(f"\n{tag}:A 组平均 **{a_:+.4f}** · B 组平均 **{b_:+.4f}** · "
          f"**A−B = {g_:+.4f} ± {sd:.4f}({abs(g_)/max(sd,1e-9):.1f}×)**")
    if tag.startswith('★'): GAP,GSD,AA,BB=g_,sd,a_,b_
nul=[contrast(rng.permutation(RES),A_HO,B_HO)[0] for _ in range(30)]
print(f"  置换 `c3⊥D` 的零(保留集):{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
za=np.zeros(NN)
for c in A_HO:
    y=pd.to_numeric(d[c],errors='coerce').values.astype(float); m=np.isfinite(y)
    za[m]+=(y[m]-np.nanmean(y))/np.nanstd(y)
n_=rng.standard_normal(NN)
y1=np.where(ok,0.5*za/len(A_HO)+n_,np.nan); y2=np.where(ok,n_,np.nan)
p1=contrast(y1,A_HO,B_HO)[0]; p2=contrast(y2,A_HO,B_HO)[0]
print(f"\n正对照两端:只贴 A 组 -> A−B **{p1:+.4f}** · 纯噪声 -> **{p2:+.4f}**")
T=pd.DataFrame([dict(set_name='全部12',a=contrast(RES,A_ALL,B_ALL)[1],b=contrast(RES,A_ALL,B_ALL)[2]),
                dict(set_name='保留集8',a=AA,b=BB)])
check_columns(T,'R306'); T.to_csv(pathlib.Path(__file__).parent/'results'/'groups.csv',index=False)

g=Gate('`c3⊥D` 是不是「对非真人的表征更容易被唤起」')
g.asserted('⚠ 分组与预测已在跑之前冻结并提交(见 PREREGISTRATION.md,独立 commit)',
           (pathlib.Path(__file__).parent/'PREREGISTRATION.md').exists(),'预注册文件存在')
g.asserted('正对照两端:只贴 A 组必须分开、纯噪声必须不分开',
           p1>2*GSD and abs(p2)<abs(p1)/2, f"① {p1:+.4f} · ② {p2:+.4f} · 2×展布 {2*GSD:.4f}")
g.negative_control('置换 `c3⊥D`(保留集)',abs(float(np.mean(nul))),abs(GAP),
                   null_spread=float(np.std(nul)),null_kind='跨人置换 —— 只打掉配对')
g.asserted('★ 注册的 kill(以保留集为准):两组反号且差 > 2×展布 -> 名字挣到第一份证据',
           AA>0 and BB<0 and GAP>2*GSD,
           f"保留集 A {AA:+.4f} · B {BB:+.4f} · 差 {GAP:+.4f} ± {GSD:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
