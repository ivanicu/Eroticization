import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A74 R303 -- c3 的名字,换一把仪器还挣不挣得到

`#257c`:c3 是这个项目现在**最结实却最没被解释**的一个 ——
29 个结局里打中 **17–21**,信度 **0.5631**,三重控制后仍然如此,对羞耻 **+0.1286**(比位置分还大),
**而它至今没有名字**(`#230c`)。
`#274` 的块载荷:高端 nonconsent · urine · transformations · abnormal body states,
低端 sensations · sex acts —— 读法是「**对越轨类别宽,而不是对普通性行为宽**」。

⚠ 这个读法是**从块载荷看出来的**,所以在块上再量一次不是独立检验(`#231b` 的教训)。
**换一把仪器**:起始仪器的 31 个类别与块选项**完全不相交**(`#235` 已断言)。

ESTIMAND        在起始仪器上按**类别流行度**切「普通/越轨」两半,各建一个背书率分,
                判**哪一半与 c3 相关**。
KILL            **若只有越轨半与 c3 相关(差 > 2×展布)-> 名字挣到跨仪器的第一份证据;
                若两半相当 -> 名字被杀,c3 仍是一个没有名字的真实维度。**
⚠ 常驻混杂        **起始仪器上的总类别数** —— 必须先报 `corr(c3, 总类别数)`,并控制后重报。
POSITIVE CTRL   两端(`#276` 同款):① 只贴越轨半的合成分必须被分开;
                ② 与两半等相关的必须不被分开。
NEGATIVE CTRL   跨人置换 c3。
IMPOSSIBLE      「越轨」= **类别流行度低**,是粗代理,不区分「社会禁忌」与「单纯少见」(`#231c` 同款)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
OBS=np.column_stack([np.isfinite(d[c].map(BIN).values.astype(float)) for c in ons]).astype(float)
NC=OBS.shape[1]; PREV=OBS.mean(0); okO=OBS.sum(1)>=8
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MB=[]; BCOL=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p_:i for i,p_ in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    MB.append((M,ppl)); BCOL.append(str(q.col))
NB=len(MB); cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
okB=cov>=8
assert not (set(BCOL)&set(ons)), "两把仪器共享了列"
def c3_of(seed=500):
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
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&okB
            if mm.sum()>300: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0)
    num=(Vv[:,2][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,2])[:,None]).sum(0)
    return np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan)
oo=np.argsort(-PREV); COM,TRG=oo[:NC//2],oo[NC//2:]
NCAT=np.where(okO,OBS.sum(1),np.nan)
SC={'普通半':np.where(okO,OBS[:,COM].mean(1),np.nan),'越轨半':np.where(okO,OBS[:,TRG].mean(1),np.nan)}
print(f"起始仪器 {NC} 类别(与块交集 0);普通半流行度 {PREV[COM].mean():.3f} · "
      f"越轨半 {PREV[TRG].mean():.3f}")
rng=np.random.default_rng(20260804)
C3=c3_of()
m0=np.isfinite(C3)&okO
print(f"⚠ 常驻混杂先报:corr(c3, 起始仪器总类别数) = "
      f"**{np.corrcoef(C3[m0],NCAT[m0])[0,1]:+.4f}**(n = {int(m0.sum()):,})")
def cr(x,y,ctrl=None):
    m=np.isfinite(x)&np.isfinite(y)&m0
    if ctrl is None:
        r=float(np.corrcoef(x[m],y[m])[0,1])
    else:
        m&=np.isfinite(ctrl)
        X=np.column_stack([np.ones(m.sum()),(ctrl[m]-ctrl[m].mean())/ctrl[m].std()])
        rx=x[m]-X@np.linalg.lstsq(X,x[m],rcond=None)[0]; ry=y[m]-X@np.linalg.lstsq(X,y[m],rcond=None)[0]
        r=float(np.corrcoef(rx,ry)[0,1])
    sd=float(np.std([np.corrcoef(x[i],y[i])[0,1] for i in
        (rng.choice(np.flatnonzero(m),int(m.sum()),True) for _ in range(200))]))
    return r,sd,int(m.sum())
print(f"\nc3 与两半的相关:")
res={};resc={}
for nm,v in SC.items():
    res[nm]=cr(C3,v); resc[nm]=cr(C3,v,ctrl=NCAT)
    print(f"  {nm}  **{res[nm][0]:+.4f} ± {res[nm][1]:.4f}**(n={res[nm][2]:,})· "
          f"控制类别数后 **{resc[nm][0]:+.4f}**")
gap=res['越轨半'][0]-res['普通半'][0]; gsd=float(np.hypot(res['越轨半'][1],res['普通半'][1]))
gapc=resc['越轨半'][0]-resc['普通半'][0]
print(f"  **越轨−普通 = {gap:+.4f} vs 2×展布 {2*gsd:.4f};控制类别数后 {gapc:+.4f}**")
nul=[cr(rng.permutation(C3),SC['越轨半'])[0] for _ in range(30)]
print(f"  置换 c3 的零:{np.mean(nul):+.4f} ± {np.std(nul):.4f}")
zt=(SC['越轨半']-np.nanmean(SC['越轨半']))/np.nanstd(SC['越轨半'])
zo=(SC['普通半']-np.nanmean(SC['普通半']))/np.nanstd(SC['普通半'])
n_=rng.standard_normal(NN)
y1=np.where(np.isfinite(zt),0.4*zt+n_,np.nan); y2=np.where(np.isfinite(zt)&np.isfinite(zo),0.28*(zt+zo)+n_,np.nan)
p1=cr(y1,SC['越轨半'])[0]-cr(y1,SC['普通半'])[0]; p2=cr(y2,SC['越轨半'])[0]-cr(y2,SC['普通半'])[0]
print(f"\n正对照两端:只贴越轨半 -> 差 **{p1:+.4f}** · 与两半等相关 -> **{p2:+.4f}**")
T=pd.DataFrame([dict(half=nm,r=res[nm][0],sd=res[nm][1],r_ctrl=resc[nm][0],n=res[nm][2]) for nm in SC])
check_columns(T,'R303'); T.to_csv(pathlib.Path(__file__).parent/'results'/'c3_vs_onset_halves.csv',index=False)

g=Gate('c3 的名字换一把仪器还挣不挣得到')
g.asserted('两把仪器的题目不相交',True,f"块 {NB} 列 / 起始 {NC} 列;交集 0")
g.asserted('正对照两端:只贴越轨半必须分开,两半等相关必须不分开',
           p1>2*gsd and abs(p2)<p1/2, f"① {p1:+.4f} · ② {p2:+.4f} · 2×展布 {2*gsd:.4f}")
g.negative_control('置换 c3',abs(float(np.mean(nul))),abs(res['越轨半'][0]),
                   null_spread=float(np.std(nul)),null_kind='跨人置换 c3 —— 只打掉配对')
g.offset_control('★ 越轨半 vs 普通半(c3 与它们的相关)',res['越轨半'][0],res['普通半'][0],gsd,
                 null_kind='同一条管道在普通半上的相关 —— 不是零假设,是「若名字不对,越轨半该落在哪」')
g.asserted('★ 注册的 kill:只有越轨半与 c3 相关 -> 名字挣到跨仪器证据;两半相当 -> 名字被杀',
           gap>2*gsd and gapc>0,
           f"越轨 {res['越轨半'][0]:+.4f} · 普通 {res['普通半'][0]:+.4f} · 差 {gap:+.4f};"
           f"控制类别数后差 {gapc:+.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
