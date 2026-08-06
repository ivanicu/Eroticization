import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A44 R248 -- 去看那个维度本身:载荷最高与最低的**选项**

`#202c`:内容维度有完整的行为画像,**没有名字**;而两次起名都被数据否掉。
`#202` 的 NEXT:**名字撑不住可能因为我一直在结局侧找它** ——
六道结局是它的**影子**,不是它本身。**它本身是 32 个块的第一主成分,而从没人看过那个载荷。**

**这一轮不做统计,只做展示**(door ①:去看对象本身)。
但它有**一个真实的方法步骤**必须验:**每块 PC1 的符号是任意的**(特征向量符号不定),
不对齐就没法跨块看载荷。

ESTIMAND        无。**这是一次展示。**
唯一的可判项    **符号对齐是否成功**:对齐后,每一块的 PC 分数与全局 `Cres` 的相关必须**全为正**。
                对齐失败 -> 载荷表不可读,整轮作废。
IMPOSSIBLE      载荷是**块内**的相对权重,跨块的绝对大小不可比;
                所以表里按块内 z 分数排,**不做跨块排名**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df_raw)
BL=[]
con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    load=v[:,-1]; sc=Z@load
    BL.append(dict(qi=q.qi,ppl=ppl,opt=opt,load=load,score=sc,rate=M.mean(0),
                   ev=float(w[-1]/w.sum())))
    con[ppl]+=sc; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
Cb=np.where(ok,con/np.maximum(cnt,1),np.nan); Sb=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
base=np.isfinite(Sb)&np.isfinite(Cb)&np.isfinite(KB); bi=np.flatnonzero(base)
X0=np.c_[np.ones(len(bi)),Sb[bi]]
Cres=np.full(NN,np.nan); Cres[bi]=Cb[bi]-X0@np.linalg.lstsq(X0,Cb[bi],rcond=None)[0]
print(f"块 {len(BL)} · n {len(bi):,} · 块内 PC1 解释率 中位 {100*np.median([b['ev'] for b in BL]):.1f}%")

# ---- 符号对齐(唯一的可判项)-------------------------------------------------
befores=[];afters=[]
for b in BL:
    m=np.isin(b['ppl'],bi)
    r=float(np.corrcoef(b['score'][m],Cres[b['ppl'][m]])[0,1]); befores.append(r)
    if r<0: b['load']=-b['load']; b['score']=-b['score']
    afters.append(abs(r) if r<0 else r)
print(f"对齐前与 Cres 相关为负的块:{sum(1 for r in befores if r<0)}/{len(BL)}")
print(f"对齐后全为正:{all(a>0 for a in afters)} · 相关范围 {min(afters):.3f}..{max(afters):.3f}")

rows=[]
for b in BL:
    z=(b['load']-b['load'].mean())/max(b['load'].std(),1e-9)
    for o,l,zz,rt in zip(b['opt'],b['load'],z,b['rate']):
        rows.append(dict(qi=b['qi'],option=str(o)[:60],loading=float(l),z_in_block=float(zz),
                         endorse_rate=float(rt)))
T=pd.DataFrame(rows); check_columns(T,'R248'); T.to_csv(pathlib.Path(__file__).parent/'results'/'loadings.csv',index=False)
print(f"\n{'='*74}\n载荷最**高**的 18 个选项(块内 z,对齐后)\n{'='*74}")
for _,r in T.sort_values('z_in_block',ascending=False).head(18).iterrows():
    print(f"  z {r.z_in_block:>+5.2f}  勾选率 {100*r.endorse_rate:>5.1f}%  {r.option[:56]}")
print(f"\n{'='*74}\n载荷最**低**的 18 个选项\n{'='*74}")
for _,r in T.sort_values('z_in_block').head(18).iterrows():
    print(f"  z {r.z_in_block:>+5.2f}  勾选率 {100*r.endorse_rate:>5.1f}%  {r.option[:56]}")

g=Gate('载荷表可不可读')
g.asserted('唯一可判项:符号对齐后每块与 Cres 全为正',all(a>0 for a in afters),
           f"{len(BL)}/{len(BL)} 为正,范围 {min(afters):.3f}..{max(afters):.3f}")
g.asserted('本轮不做统计,只做展示 —— 不注册 kill',True,
           'door ①:在影子上猜了两轮(#201/#202)之后,去看对象本身')
g.asserted('⚠ 载荷是块内相对权重,跨块绝对大小不可比',True,'表里按块内 z 排,不做跨块排名')
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
