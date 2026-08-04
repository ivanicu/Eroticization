import os,sys,pathlib,time
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A94 R342 -- 用一个能通过正对照的结局侧代理,再问一次

`#296a`:答题风格代理**收回 −18.5%**,抓不住已知的信度差,所以它的「没抓到」什么也不说明。
它失败的原因是结构性的:**答题风格是人层的三个数,按定义压不动逐结局的信度。**

**这一轮换成逐结局的量:人层自助的 R² 噪声底 `sd_k`** ——
对每个结局,在**每组内部**重抽人 B 次,记该格六坐标 R² 的自助标准差。
那**直接就是那一格的分辨率下限**,而它随该组在该结局上的 n · 全距 · 信度而变。

ESTIMAND        逐结局 `sd_k` 曲线(29 维,每组一条);用它偏出六坐标 R² 曲线后重测
                「按性别劈 vs 随机劈」。
⛔ GATE         **正对照必须先过,否则第二臂不许读**(`#296a` 的教训):
                ① 给一组的 6 个结局加噪 -> `sd_k` **必须在那 6 格上升**(单独报,它决定代理合不合格);
                ② 偏出后人为差**收回 > +40%**;③ guard 17 先挡收回率的定义域。
KILL            **若偏出后差塌掉 -> `#286`–`#295` 测的是结局侧的测量差;
                若差还在(且正对照过了)-> 结构差是真的,那十轮升到 D6。**
NEGATIVE CTRL   由 offset 承担(随机劈,**同样偏出**)。
IMPOSSIBLE      自助固定了六坐标(只重抽人),所以 `sd_k` 是**给定坐标下**那格 R² 的抽样噪声,
                不含坐标本身的估计噪声 —— 它是结局侧的下限,不是全部噪声。
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
NB=len(MB); cov=np.zeros(NN); pos=np.zeros(NN)
for M,ppl in MB:
    rr=-np.log(np.clip(M.mean(0),1e-4,1.)); n=M.sum(1)
    cov[ppl]+=1; pos[ppl]+=np.where(n>0,(M@rr)/np.maximum(n,1),0.0)
ok=cov>=8; S=np.where(ok,pos/np.maximum(cov,1),np.nan)
SEX=pd.to_numeric(d.get('biomale'),errors='coerce').values.astype(float)
RATE=np.array([M.mean() for M,_ in MB]); o_=np.argsort(-RATE); ORD,TRG=o_[:NB//2],o_[NB//2:]
rg=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    o=rg.permutation(M.shape[1]); k=M.shape[1]//2
    A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
def coords(rows):
    m=np.zeros(NN,bool); m[rows]=True
    def prof_(X):
        F=np.isfinite(X)&m[None,:]; Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
        for b in range(NB):
            lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
            R[b]=np.where(F[b],X[b]-lo,np.nan); R[b]=R[b]-np.nanmean(R[b])
        return R
    Ra,Rb=prof_(A),prof_(B)
    st=np.full(NN,np.nan); G=np.isfinite(Ra)&np.isfinite(Rb)
    for i in rows:
        mm=G[:,i]
        if mm.sum()<8: continue
        x,y_=Ra[mm,i],Rb[mm,i]
        if x.std()>1e-9 and y_.std()>1e-9: st[i]=float(np.corrcoef(x,y_)[0,1])
    def z(v):
        mk=np.isfinite(v)&m; w=np.full(NN,np.nan); w[mk]=(v[mk]-v[mk].mean())/v[mk].std(); return w
    def hs(R,cols):
        sub=R[cols]; F2=np.isfinite(sub)
        return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
    D=(z(hs(Ra,TRG))-z(hs(Ra,ORD))+z(hs(Rb,TRG))-z(hs(Rb,ORD)))/2
    C=np.full((NB,NB),np.nan)
    for i in range(NB):
        for j in range(NB):
            mm=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
            if mm.sum()>200: C[i,j]=np.corrcoef(Ra[i][mm],Rb[j][mm])[0,1]
    C=np.where(np.isfinite(C),C,0.0); C=(C+C.T)/2; Vv=np.linalg.eigh(C)[1][:,::-1]
    R=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
    R=np.where(np.isfinite(Ra)|np.isfinite(Rb),R,np.nan)
    Fm=np.isfinite(R); Zm=np.where(Fm,R,0.0); csx=[]
    for k in range(3):
        num=(Vv[:,k][:,None]*Zm).sum(0); den=(Fm*np.abs(Vv[:,k])[:,None]).sum(0)
        csx.append(np.where(den>1e-9,num/np.maximum(den,1e-9),np.nan))
    return [S,D]+csx+[st]
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
lik=[c for c in lik if c!='biomale']
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
GEN=[c for c in [x for x,_ in OUT] if any(k in str(c) for k in
     ('opposite gender','biological *female*','biological *male*',
      'biological female','biological male'))]
SEEN=[c for c in [x for x,_ in OUT] if str(c).startswith('"I sometimes find people') or c=='animated']
NONGEN=[c for c,_ in OUT if c not in GEN]
HOLD=[c for c in NONGEN if c not in SEEN]
print(f"⚠ 分组已在跑之前冻结提交(见 PREREGISTRATION.md):")
print(f"   性别指涉 **{len(GEN)}** 道 · 非性别指涉 **{len(NONGEN)}** 道 · "
      f"**保留集(剔除 `#286a` 已看过的 2 格)= {len(HOLD)} 道**")

NBOOT=int(os.environ.get('NBOOT','40'))   # ⚠ 不能叫 B —— head 里 B 是选项劈分的后半矩阵
_C={}
def coords6(rows):
    k=rows.tobytes()
    if k not in _C: _C[k]=coords(rows)
    return _C[k]
def r2vec(rows,Q,outs):
    out=[]
    for nm,y in outs:
        m=np.isfinite(y[rows])
        for q_ in Q: m&=np.isfinite(q_[rows])
        r=rows[m]
        if len(r)<250: out.append(np.nan); continue
        X=np.column_stack([np.ones(len(r))]+[(q_[r]-q_[r].mean())/max(q_[r].std(),1e-12) for q_ in Q])
        yy=y[r]; yy=(yy-yy.mean())/max(yy.std(),1e-12)
        b=np.linalg.lstsq(X,yy,rcond=None)[0]
        out.append(float(1-np.var(yy-X@b)/max(np.var(yy),1e-12)))
    return np.array(out)
def curve(rows,outs=None): return r2vec(rows,coords6(rows),outs if outs is not None else OUT)
def sdcurve(rows,outs=None,seed=0):
    """★ 逐结局的人层自助噪声底 —— 那一格 R² 的分辨率下限,是结局侧的量。"""
    Q=coords6(rows); oo=outs if outs is not None else OUT
    rg=np.random.default_rng(seed); acc=[]
    for _ in range(NBOOT):
        bs=rows[rg.integers(0,len(rows),len(rows))]
        acc.append(r2vec(bs,Q,oo))
    return np.nanstd(np.array(acc),0)
def sim(a,b):
    m=np.isfinite(a)&np.isfinite(b); return float(np.corrcoef(a[m],b[m])[0,1]) if m.sum()>=10 else np.nan
def sim_p(a,b,ra,rb):
    m=np.isfinite(a)&np.isfinite(b)&np.isfinite(ra)&np.isfinite(rb)
    if m.sum()<10: return np.nan
    def rs(y,x):
        x=(x[m]-x[m].mean())/max(x[m].std(),1e-12); return y[m]-np.polyval(np.polyfit(x,y[m],1),x)
    return float(np.corrcoef(rs(a,ra),rs(b,rb))[0,1])
have=ok&np.isfinite(SEX)
g0=np.flatnonzero(have&(SEX==0)); g1=np.flatnonzero(have&(SEX==1))
rngS=np.random.default_rng(20260804); RND=[]
for _ in range(4):
    p=rngS.permutation(np.flatnonzero(have)); RND.append((p[:len(g0)].copy(),p[len(g0):len(g0)+len(g1)].copy()))
t0=time.time(); _=sdcurve(g0,seed=1); print(f"NBOOT={NBOOT} 单条 sd 曲线 {time.time()-t0:.1f}s;两组 n={len(g0):,}/{len(g1):,}")

# ---- ⛔ 门:正对照先跑,不过就不读第二臂 ----
HIT=set(range(0,len(OUT),5)); rgN=np.random.default_rng(4242)
def noisy(rows,lvl):
    m=np.zeros(NN,bool); m[rows]=True; o=[]
    for k,(nm,y) in enumerate(OUT):
        yy=y.copy()
        if k in HIT and lvl>0:
            z=lvl*np.nanstd(y)*np.random.default_rng(700+k).standard_normal(NN)
            yy=np.where(m,y+z,y)
        o.append((nm,yy))
    return o
print(f"\n⛔ 门 ①:给**组1**的 {len(HIT)} 个结局加噪 -> `sd_k` 必须**在那 6 格上升**")
hb=np.array([k in HIT for k in range(len(OUT))])
sd_lo=sdcurve(g1,noisy(g1,0.0),seed=11)
SDS=[]
for lvl in (0.0,0.35,0.7,1.4,2.0):     # ⚠ 用 guard 13 的**扫描**验证代理,不是两点(#297b)
    sd_=sdcurve(g1,noisy(g1,lvl),seed=11)
    rh=float(np.nanmean(sd_[hb])/np.nanmean(sd_lo[hb]))
    r2_=float(np.nanmean(curve(g1,noisy(g1,lvl))[hb]))
    SDS.append((lvl,rh))
    print(f"   lvl={lvl:.2f}: 加噪格 sd 比 **{rh:.3f}×**  (那 6 格的平均 R² = {r2_:.4f})")
sd_hi=sdcurve(g1,noisy(g1,2.0),seed=11)
rat_hit=float(np.nanmean(sd_hi[hb])/np.nanmean(sd_lo[hb]))
rat_oth=float(np.nanmean(sd_hi[~hb])/np.nanmean(sd_lo[~hb]))
print(f"   -> lvl=2 时 加噪格 **{rat_hit:.3f}×** · 其余格 **{rat_oth:.3f}×** -> "
      f"选择性 **{rat_hit/max(rat_oth,1e-9):.2f}×**")
print(f"\n⛔ 门 ②:偏出后人为差必须收回 > +40%")
for lvl in (0.0,2.0):
    oo0,oo1=noisy(g1,0.0),noisy(g1,lvl)
    a,b=curve(g0,oo0),curve(g1,oo1); ra,rb=sdcurve(g0,oo0,seed=21),sdcurve(g1,oo1,seed=22)
    raw_,par_=sim(a,b),sim_p(a,b,ra,rb)
    print(f"   lvl={lvl:.1f}: 原样 **{raw_:+.4f}** · 偏出后 **{par_:+.4f}**")
    if lvl==0.0: b_raw,b_par=raw_,par_
    else: h_raw,h_par=raw_,par_
rec=abs(h_par-b_par)/max(abs(h_raw-b_raw),1e-9); RECOV=1-rec
print(f"   原样掉 **{h_raw-b_raw:+.4f}** · 偏出后掉 **{h_par-b_par:+.4f}** -> 收回 **{100*RECOV:+.1f}%**")
GATE_OK=(rat_hit/max(rat_oth,1e-9)>1.20) and (RECOV>0.40) and (-2<=RECOV<=2)
print(f"\n⛔ 门:{'**过** —— 允许读第二臂' if GATE_OK else '**没过** —— 第二臂不读,整轮 UNVERIFIED'}")
RES={}
if GATE_OK:
    def arm(fn,tag):
        s=fn(g0,g1); rr=[fn(a,b) for a,b in RND]
        rm,rd=float(np.mean(rr)),float(np.std(rr))
        print(f"  {tag:<22} 性别劈 **{s:+.4f}** vs 随机劈 **{rm:+.4f} ± {rd:.4f}** -> "
              f"差 **{s-rm:+.4f}**({abs(s-rm)/max(2*rd,1e-9):.1f}× 的 2×展布)")
        return s,rm,rd
    print()
    RES['raw']=arm(lambda a,b: sim(curve(a),curve(b)),'① 原样(`#286a`)')
    RES['par']=arm(lambda a,b: sim_p(curve(a),curve(b),sdcurve(a,seed=31),sdcurve(b,seed=32)),
                   '② **偏出噪声底后**')
    d0=RES['raw'][0]-RES['raw'][1]; d1=RES['par'][0]-RES['par'][1]
    print(f"\n偏出后差从 **{d0:+.4f}** 变成 **{d1:+.4f}** —— 保留 **{100*d1/d0:.1f}%**")
    print(f"⚠ 展布:**{RES['raw'][2]:.4f}** -> **{RES['par'][2]:.4f}**({100*RES['par'][2]/RES['raw'][2]:.0f}%)")
T=pd.DataFrame([dict(gate='①sd 选择性',v=rat_hit/max(rat_oth,1e-9)),dict(gate='②收回率',v=RECOV)]+
               [dict(gate=k,v=v[0]-v[1]) for k,v in RES.items()])
check_columns(T,'R342'); T.to_csv(pathlib.Path(__file__).parent/'results'/'noise_floor.csv',index=False)
gg=Gate('用一个能通过正对照的结局侧代理,再问一次')
gg.bounded_statistic_out_of_range('⚠ guard 17:收回率先过自己的定义域',RECOV,-2,2,'收回率')
gg.plant_direction_from_sweep('⚠ 代理的单调性(guard 13 扫描,不是两点)',SDS,SDS[0][1])
gg.asserted('⛔ 门①:加噪格的 `sd_k` 必须选择性上升(>1.20×)',rat_hit/max(rat_oth,1e-9)>1.20,
            f"加噪格 {rat_hit:.3f}× · 其余 {rat_oth:.3f}× -> 选择性 {rat_hit/max(rat_oth,1e-9):.2f}×")
gg.asserted('⛔ 门②:偏出后人为差收回 > +40%',RECOV>0.40,
            f"原样掉 {h_raw-b_raw:+.4f} -> 偏出后掉 {h_par-b_par:+.4f}(收回 {100*RECOV:+.1f}%)")
if GATE_OK:
    for k,lab in (('raw','① 原样'),('par','★ ② **偏出噪声底后**')):
        s,rm,rd=RES[k]
        gg.offset_control(f'{lab}:按性别劈 vs 随机劈',s,rm,rd,
            null_kind='随机劈同样大小两组,**同样偏出** —— 不是零假设,是「若结构相同该落在哪」')
    gg.component_difference_is_not_mechanism('★ 注册的 kill:偏出结局侧噪声底后差是否塌掉',
        RES['raw'][0]-RES['raw'][1],RES['par'][0]-RES['par'][1],2*RES['par'][2],'结局侧噪声底')
else:
    gg.asserted('★ 注册的 kill',False,'门没过 -> 不评估(P16:能在坏仪器上开火的 kill 不是承诺)')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
