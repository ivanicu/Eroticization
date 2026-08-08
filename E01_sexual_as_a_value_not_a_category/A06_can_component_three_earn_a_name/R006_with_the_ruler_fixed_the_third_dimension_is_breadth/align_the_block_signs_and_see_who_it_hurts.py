import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A56 R272 -- 把块的主成分符号对齐,然后当场判它伤到了谁

`#226b`:`np.linalg.eigh` 的特征向量符号任意,而 `R210:73` 的 `con[ppl] += Z@pc`
把 32 个块的主成分人分**直接相加**。实测 17 个块载荷为正、15 个为负 ——
**"内容分"是一个 ±breadth 的相消混合。**

ESTIMAND        用**共识法留一对齐**块符号后的内容人分,与未对齐版本并排,回答三问:
                ① 对齐后分半信度高多少;② `corr(内容分, 勾选数)` 是否说明它该叫 breadth;
                ③ 用对齐分数重跑 31 个结局,`#189` 那一支还在不在。
KILL(逐问预注册)
                Q1 若对齐后信度**没有**明显高于未对齐 -> `#226b` 只是一个无后果的洁癖,
                   "尺子装反"这个说法要收窄成"符号未对齐,但不影响结果"。
                Q2 若 `|corr(内容分, 勾选数)| > 0.5` -> 这个维度的正名是 **breadth**,
                   `#188`/`#189` 的"第三个维度"要按这个名字重写。
                Q3 若对齐后越阈值结局数**不增反降** -> 对齐不是修复,是另一种破坏。
POSITIVE CTRL   (这一轮唯一合法的正对照)种入一个内容信号,**再人为打乱块符号**;
                **对齐法必须把它救回来,未对齐法必须救不回来。**
                只有这两条同时成立,"对齐"才被证明是在做它声称的事。
⚠ 留一,不然是一个不可能失败的检查
                共识分若包含块 i 自己,块 i 与它必然正相关,符号永远不翻 ——
                所以每个块只与**其余块**的共识分比对。
IMPOSSIBLE      每一半内独立对齐后仍有一个整体符号自由度 -> 相关取绝对值,
                这会让分半信度**偏高**(选择在样本相关上),对齐/未对齐**同样偏**,
                所以可比的是**差**,不是各自的绝对值。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); NN=len(d)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BL=[]; pos=np.zeros(NN); K=np.zeros(NN); cnt=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    BL.append(dict(M=M,ppl=ppl,rar=rr))
    pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); K[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8; S=np.where(ok,pos/np.maximum(cnt,1),np.nan); K=np.where(ok,K,np.nan)
print(f"块 {len(BL)} 个;n = {int(ok.sum()):,}")

def block_z(blocks, plant=0.0, u=None, scramble_seed=None):
    """每块的主成分人分(NN 长,缺失为 nan)。scramble_seed 不为 None 时人为打乱块符号。"""
    Zs=[]
    for j in blocks:
        b=BL[j]; M=b['M']; ppl=b['ppl']
        if plant and u is not None:
            sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
            M=M+plant*np.outer(u[ppl],sub)
        Z=M-M.mean(0,keepdims=True)
        w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); pc=v[:,-1]
        if scramble_seed is not None:
            if np.random.default_rng(zlib.crc32(f'{j}{scramble_seed}'.encode())%(1<<30)).random()<0.5: pc=-pc
        z=np.full(NN,np.nan); z[ppl]=Z@pc; Zs.append(z)
    return np.array(Zs)

def combine(Zs, align=True, iters=5):
    """align=False 复现 R210:73 的直接相加;align=True 用共识法【留一】对齐符号。"""
    sg=np.ones(len(Zs)); F=np.isfinite(Zs)
    if align:
        for _ in range(iters):
            tot=np.nansum(Zs*sg[:,None],0); ct=F.sum(0); flipped=0
            for i in range(len(Zs)):
                lo_tot=tot-np.where(F[i],Zs[i]*sg[i],0.0); lo_ct=ct-F[i]      # ⚠ 留一
                m=np.where(lo_ct>=3, lo_tot/np.maximum(lo_ct,1), np.nan)
                g=F[i]&np.isfinite(m)
                if g.sum()>100:
                    c=np.corrcoef(Zs[i][g],m[g])[0,1]
                    if np.isfinite(c) and c<0: sg[i]=-sg[i]; flipped+=1
                    tot=np.nansum(Zs*sg[:,None],0)
            if flipped==0: break
    tot=np.nansum(Zs*sg[:,None],0); ct=F.sum(0)
    return np.where(ct>=3, tot/np.maximum(ct,1), np.nan), sg

def reliability(align, seeds=10, plant=0.0, u=None, scramble_seed=None):
    vs=[]
    for s in range(seeds):
        rg=np.random.default_rng(700+s); p=rg.permutation(len(BL)); h=len(BL)//2
        a,_=combine(block_z(p[:h],plant,u,scramble_seed),align)
        b,_=combine(block_z(p[h:2*h],plant,u,scramble_seed),align)
        m=np.isfinite(a)&np.isfinite(b)&ok
        if m.sum()<500: continue
        r=abs(float(np.corrcoef(a[m],b[m])[0,1])); vs.append(2*r/(1+r) if r<0.999 else np.nan)
    return float(np.nanmean(vs)), float(np.nanstd(vs))

# ---------- 正对照:种入内容信号 + 人为打乱块符号 ----------
rg=np.random.default_rng(20260804); u=rg.standard_normal(NN)   # 全员种入,信度本来就只在 ok 上算
pa=reliability(True ,seeds=5,plant=0.5,u=u,scramble_seed=77)
pu=reliability(False,seeds=5,plant=0.5,u=u,scramble_seed=77)
print(f"\n正对照(种入内容信号 g=0.5 + 人为打乱块符号):"
      f"对齐 {pa[0]:.4f} ± {pa[1]:.4f} · 未对齐 {pu[0]:.4f} ± {pu[1]:.4f}")

# ---------- Q1 真实数据上的信度 ----------
ra=reliability(True); ru=reliability(False)
print(f"\nQ1 内容分半信度:**对齐 {ra[0]:.4f} ± {ra[1]:.4f}** · 未对齐 {ru[0]:.4f} ± {ru[1]:.4f}"
      f" · 差 **{ra[0]-ru[0]:+.4f}**")

# ---------- Q2 它是不是 breadth ----------
Za=block_z(range(len(BL))); Ca,sg=combine(Za,True); Cu,_=combine(Za,False)
mm=np.isfinite(Ca)&np.isfinite(K)&ok
cKa=float(np.corrcoef(Ca[mm],K[mm])[0,1]); cSa=float(np.corrcoef(Ca[mm],S[mm])[0,1])
mu=np.isfinite(Cu)&np.isfinite(K)&ok
cKu=float(np.corrcoef(Cu[mu],K[mu])[0,1])
print(f"\nQ2 翻了 {int((sg<0).sum())}/{len(BL)} 个块的符号")
print(f"   corr(对齐内容分, 勾选数) = **{cKa:+.4f}** · corr(未对齐, 勾选数) = {cKu:+.4f}")
print(f"   corr(对齐内容分, 位置分 S) = {cSa:+.4f}")

# ---------- Q3 31 个结局 ----------
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
def resid(x,z):
    m=np.isfinite(x)&np.isfinite(z); o=np.full(NN,np.nan)
    b=np.polyfit(z[m],x[m],1); o[m]=x[m]-np.polyval(b,z[m]); return o
def npass(x,rng,tag):
    bi=np.flatnonzero(np.isfinite(x)); res=[]; nl=[]
    for nm,y in OUT:
        m=np.isfinite(y[bi]); jj=bi[m]
        if len(jj)<200: continue
        res.append((nm,float(np.corrcoef(y[jj],x[jj])[0,1])))
        nl.append([abs(float(np.corrcoef(rng.permutation(y[jj]),x[jj])[0,1])) for _ in range(20)])
    L=min(len(z) for z in nl)
    thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl]),axis=0),0.95))
    hit=[(n,v) for n,v in res if abs(v)>thr]
    print(f"   {tag:<22} {len(hit):>2}/{len(res)}(阈值 {thr:.4f}) "
          + ' · '.join(f"{n[:22]} {v:+.3f}" for n,v in sorted(hit,key=lambda t:-abs(t[1]))[:5]))
    D=dict(res)
    print(f"     `#189` 的三个具名结局:"+' · '.join(
        f"{k} {D.get(k,float('nan')):+.4f}{'✓' if abs(D.get(k,0))>thr else '✗'}"
        for k in ('animated','written','biomale')))
    return len(hit),len(res),thr
rng2=np.random.default_rng(4242)
print(f"\nQ3 31 个结局(逐轮重算 max-stat 阈值):")
ha,na_,ta=npass(Ca,rng2,'对齐 内容分')
hu,nu_,tu=npass(Cu,rng2,'未对齐 内容分(R210)')
Cra=resid(Ca,S); Cru=resid(Cu,S)
mr=np.isfinite(Cra)&np.isfinite(K)&ok; mru=np.isfinite(Cru)&np.isfinite(K)&ok
cKra=float(np.corrcoef(Cra[mr],K[mr])[0,1]); cKru=float(np.corrcoef(Cru[mru],K[mru])[0,1])
print(f"   ⚠ 扣掉 S 之后 corr(Cres, 勾选数):**对齐 {cKra:+.4f}** · 未对齐(`#189` 用的) {cKru:+.4f}")
hra,_,tra=npass(Cra,rng2,'对齐 Cres(扣掉 S)')
hru,_,tru=npass(Cru,rng2,'未对齐 Cres(#189)')

T=pd.DataFrame([dict(which='aligned',rel=ra[0],rel_sd=ra[1],corr_picks=cKa,corr_S=cSa,corr_picks_res=cKra,
                     n_pass=ha,n_pass_res=hra,n_out=na_,flipped=int((sg<0).sum())),
                dict(which='unaligned',rel=ru[0],rel_sd=ru[1],corr_picks=cKu,corr_S=np.nan,corr_picks_res=cKru,
                     n_pass=hu,n_pass_res=hru,n_out=nu_,flipped=0)])
check_columns(T,'R272'); T.to_csv(pathlib.Path(__file__).parent/'results'/'aligned_vs_not.csv',index=False)

g=Gate('把块符号对齐,判它伤到了谁')
g.asserted('正对照:种入内容信号+打乱符号 -> 对齐法救回来,未对齐法救不回来',
           (pa[0]>pu[0]+0.10) and (pu[0]<0.30), f"对齐 {pa[0]:.4f} vs 未对齐 {pu[0]:.4f}")
g.offset_control('Q1 对齐 vs 未对齐(真实数据分半信度)',ra[0],ru[0],float(np.hypot(ra[1],ru[1])),
                 null_kind='同一条管道不对齐时的信度 —— 不是零假设,是"若符号本来就齐,对齐应当什么都不改变"')
g.asserted('Q2 kill:|corr(内容分, 勾选数)| > 0.5 -> 这个维度的正名是 breadth',
           abs(cKa)>0.5, f"corr = {cKa:+.4f}(未对齐时 {cKu:+.4f})")
g.asserted('Q2b 扣掉 S 之后仍 > 0.5 -> 连 `#189` 用的 Cres 也是 breadth,不是内容',
           abs(cKra)>0.5, f"corr(Cres_对齐, 勾选数) = {cKra:+.4f}(未对齐 {cKru:+.4f})")
g.asserted('Q3 kill:对齐后越阈值结局数不增反降 -> 对齐不是修复',
           ha>=hu, f"对齐 {ha}/{na_} vs 未对齐 {hu}/{nu_};扣掉 S 后 {hra} vs {hru}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
