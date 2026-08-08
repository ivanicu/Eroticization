import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A95 R344 -- 把 `#298c` 的九成拿回来:最不一致的那些选项

`#298c`:`corr(S_global, S_local)` = **+0.8999**,所以 `#298` 只说得了残差的那一成。
**但重合不等于不可分** —— 两者的差别集中在**参照群体流行度与全样本流行度差得最多的那些选项**上。

ESTIMAND        逐选项算跨参照群体的流行度离散 `sd_cells(p)`,取每块内**最不一致的前 f**,
                只用它们重算 `S_g`/`S_l`,再问 Ⓐ/Ⓑ 谁带羞耻。
⚠ 劈半          **选项是看着数据挑的** -> 在**一半人**上挑,在**另一半人**上读(`#254a`),两向都报。
KILL            **若子集下 `corr(S_g,S_l)` 明显降(可分性真的提高)且 Ⓐ 仍显著大于 Ⓑ ->
                内容污名的读法扩到更大的一块;
                若可分性没提高 -> 这条路到此为止,`#298` 的一成就是这个设计的上限。**
POSITIVE CTRL   `S_local` 驱动的合成结局,**必须在这个子集上重跑**(子集变了仪器就变了 `#252c`)。
⚠ CONTROL       子集会让很多人在某块**一个都没选** -> 覆盖下降 -> guard 12 必须过(`#239a`)。
⚠ KNOB          f ∈ {0.10(注册的主口径), 0.25, 0.50} -> `CALIBER.md`。
IMPOSSIBLE      参照群体是观察到的;而且「最不一致」是**用同一份数据定义的**,劈半只能控住选择偏倚,
                控不住「这些选项恰好也是别的东西的标记」。
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

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
SXL=np.where(np.isfinite(SEX),np.where(SEX==1,'M','F'),'?')
CELL=np.array([f"{a}|{s}" for a,s in zip(d['age'].astype(str).values,SXL)],dtype=object)
LABS=[c for c in pd.unique(CELL) if isinstance(c,str)]

def divergence(sel_rows):
    """逐选项:跨参照群体的流行度离散。**只用 sel_rows 的人算** —— 这是「挑」的那一半。"""
    m=np.zeros(NN,bool); m[sel_rows]=True; out=[]
    for M,ppl in MB:
        keep=m[ppl]; P=[]
        for L in LABS:
            s2=keep&(CELL[ppl]==L)
            if s2.sum()>=150: P.append(M[s2].mean(0))
        out.append(np.std(np.array(P),0) if len(P)>=4 else np.zeros(M.shape[1]))
    return out

def S_pair(rows,pick,cells_on,est=None):
    """在 rows 上读,但**流行度在 est 上估**(默认全样本)。

    ⚠ 劈半只该管「挑哪些选项」这一步 —— 那是唯一看着数据做的选择。
    把**流行度**也劈了半会让 10 个参照格各自不足 200,覆盖塌到 17%,
    而那不是控制,是把仪器弄坏了(第一版就是这么坏的)。
    选择用的是流行度离散,**完全不碰结局**,所以这里没有结局泄漏可言。"""
    m=np.zeros(NN,bool); m[rows]=True
    me=np.ones(NN,bool) if est is None else np.zeros(NN,bool)
    if est is not None: me[est]=True
    cv=np.zeros(NN); ps=np.zeros(NN)
    for b,(M,ppl) in enumerate(MB):
        cols=pick[b]
        if cols.sum()<2: continue
        Ms_all=M[:,cols]
        for L in ([None] if not cells_on else LABS):
            est_s=me[ppl] if L is None else (me[ppl]&(CELL[ppl]==L))
            sel  =m[ppl]  if L is None else (m[ppl] &(CELL[ppl]==L))
            if est_s.sum()<200 or sel.sum()<50: continue
            rr=-np.log(np.clip(Ms_all[est_s].mean(0),1e-4,1.))    # ⚠ 流行度在 est 上估
            Ms=Ms_all[sel]; n=Ms.sum(1)
            v=np.where(n>0,(Ms@rr)/np.maximum(n,1),np.nan)
            ppl_s=ppl[sel]; g=np.isfinite(v)
            cv[ppl_s[g]]+=1; ps[ppl_s[g]]+=v[g]
    return np.where(cv>=6,ps/np.maximum(cv,1),np.nan),cv

def resid(a,b,m):
    out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
    out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
def cor(u,v):
    m=np.isfinite(u)&np.isfinite(v); return (float(np.corrcoef(u[m],v[m])[0,1]),int(m.sum())) if m.sum()>200 else (np.nan,0)

have=ok&np.isfinite(SEX)&np.isfinite(sh); ALLR=np.flatnonzero(have)
rgH=np.random.default_rng(31337); pp=rgH.permutation(ALLR)
HALF=[(pp[:len(pp)//2],pp[len(pp)//2:]),(pp[len(pp)//2:],pp[:len(pp)//2])]
print(f"人 n={len(ALLR):,};劈半 {len(HALF[0][0]):,} / {len(HALF[0][1]):,};参照群体 {len(LABS)} 个\n")
rows=[]
for f in (0.10,0.25,0.50,1.00):
    accs=[]
    for hi,(pick_on,read_on) in enumerate(HALF):
        DV=divergence(pick_on)
        pick=[]
        for b,(M,_) in enumerate(MB):
            k=max(2,int(round(f*M.shape[1]))) if f<1.0 else M.shape[1]
            idx=np.argsort(-DV[b])[:k]; msk=np.zeros(M.shape[1],bool); msk[idx]=True; pick.append(msk)
        Sg,cvg=S_pair(read_on,pick,False); Sl,_=S_pair(read_on,pick,True)
        mm=np.isfinite(Sg)&np.isfinite(Sl)
        rGL,_=cor(Sg,Sl); Rg,Rl=resid(Sg,Sl,mm),resid(Sl,Sg,mm)
        cg,ng=cor(Rg,sh); cl,_=cor(Rl,sh)
        accs.append((rGL,cg,cl,ng,int(mm.sum())))
    a=np.array([x[:3] for x in accs]); n_=int(np.mean([x[3] for x in accs]))
    nkeep=int(np.mean([x[4] for x in accs]))
    mde=1.96/np.sqrt(max(n_,1))
    ntot=int(np.mean([len(h[1]) for h in HALF]))
    print(f"【f={f:.2f}】 每块前 {int(100*f)}% 最不一致的选项 · 可读 n={n_:,}(覆盖 {100*nkeep/ntot:.0f}%)")
    print(f"   `corr(S_g,S_l)` = **{a[:,0].mean():+.4f}**(两向 {a[0,0]:+.4f} / {a[1,0]:+.4f})· MDE **{mde:.4f}**")
    print(f"   **Ⓐ 内容污名 {a[:,1].mean():+.4f}**(两向 {a[0,1]:+.4f} / {a[1,1]:+.4f})· "
          f"**Ⓑ 参照偏离 {a[:,2].mean():+.4f}**(两向 {a[0,2]:+.4f} / {a[1,2]:+.4f})")
    rows.append(dict(frac=f,r_gl=float(a[:,0].mean()),v_a=float(a[:,1].mean()),v_b=float(a[:,2].mean()),
                     mde=mde,n_read=n_,cover=nkeep/ntot,
                     a_spread=float(abs(a[0,1]-a[1,1])/2),b_spread=float(abs(a[0,2]-a[1,2])/2)))
T=pd.DataFrame(rows); check_columns(T,'R344')
T.to_csv(pathlib.Path(__file__).parent/'results'/'divergent_options.csv',index=False)
# ---- 正对照:在**这个子集**上重跑 S_local 驱动的合成结局(#252c) ----
pick_on,read_on=HALF[0]; DV=divergence(pick_on); pick=[]
for b,(M,_) in enumerate(MB):
    k=max(2,int(round(0.10*M.shape[1]))); idx=np.argsort(-DV[b])[:k]
    msk=np.zeros(M.shape[1],bool); msk[idx]=True; pick.append(msk)
Sg,_=S_pair(read_on,pick,False); Sl,_=S_pair(read_on,pick,True)
mm=np.isfinite(Sg)&np.isfinite(Sl); Rg,Rl=resid(Sg,Sl,mm),resid(Sl,Sg,mm)
zl=np.where(np.isfinite(Sl),(Sl-np.nanmean(Sl))/np.nanstd(Sl),np.nan)
SW=[]
for g in (0.0,0.3,0.6,1.0):
    y=g*zl+np.random.default_rng(5).standard_normal(NN)
    a_,_=cor(Rg,y); b_,_=cor(Rl,y); SW.append((g,b_))
    print(f"  正对照(**子集上重跑**)g={g:.1f}:全局残差 **{a_:+.4f}** · 局部残差 **{b_:+.4f}**")
# ⚠ guard 12 的交集重比:n 从 3,283 掉到 1,825,必须分清是**效应变了**还是**人变了**(#239a)
def AB_on(rows_read,f,common=None):
    accs=[]
    for pick_on,read_on in HALF:
        DV=divergence(pick_on); pick=[]
        for b,(M,_) in enumerate(MB):
            k=max(2,int(round(f*M.shape[1]))) if f<1.0 else M.shape[1]
            idx=np.argsort(-DV[b])[:k]; msk=np.zeros(M.shape[1],bool); msk[idx]=True; pick.append(msk)
        Sg,_=S_pair(read_on,pick,False); Sl,_=S_pair(read_on,pick,True)
        mm=np.isfinite(Sg)&np.isfinite(Sl)
        if common is not None: mm&=common
        Rg,Rl=resid(Sg,Sl,mm),resid(Sl,Sg,mm)
        accs.append((cor(Rg,sh)[0],cor(Rl,sh)[0],int(mm.sum())))
    a=np.array([x[:2] for x in accs])
    return float(a[:,0].mean()),float(a[:,1].mean()),int(np.mean([x[2] for x in accs]))
COM=np.zeros(NN,bool)
for pick_on,read_on in HALF:
    DV=divergence(pick_on); pick=[]
    for b,(M,_) in enumerate(MB):
        k=max(2,int(round(0.10*M.shape[1]))); idx=np.argsort(-DV[b])[:k]
        msk=np.zeros(M.shape[1],bool); msk[idx]=True; pick.append(msk)
    Sg,_=S_pair(read_on,pick,False); Sl,_=S_pair(read_on,pick,True)
    COM|=(np.isfinite(Sg)&np.isfinite(Sl))
aF,bF,nF=AB_on(None,1.00,COM); a1,b1,n1=AB_on(None,0.10,COM)
print(f"\n⚠ 交集重比(只看 f=0.10 也可读的那 {n1:,} 个人):")
print(f"   f=1.00 在交集上 Ⓐ **{aF:+.4f}** · Ⓑ **{bF:+.4f}**  (全样本上是 +0.0414 / +0.0147)")
print(f"   f=0.10 在交集上 Ⓐ **{a1:+.4f}** · Ⓑ **{b1:+.4f}**")
r0=T[T.frac==0.10].iloc[0]; rF=T[T.frac==1.00].iloc[0]
gg=Gate('最不一致的那些选项:可分性能不能提高')
gg.plant_direction_from_sweep('正对照:子集上 `S_local` 驱动的结局 -> 局部残差随强度上升',SW,SW[0][1])
gg.control_kept_the_sample('⚠ guard 12:子集让覆盖掉了多少(附交集重比)',
                           float(rF.v_a),float(r0.v_a),int(rF.n_read),int(r0.n_read),
                           before_common=aF,after_common=a1,n_common=n1)
gg.asserted('⚠ 规格曲线:Ⓐ−Ⓑ 在 f 上的符号是否稳定',
            len(set(np.sign(T.v_a-T.v_b)))==1,
            ' · '.join(f"f={r.frac:.2f} {r.v_a-r.v_b:+.4f}" for _,r in T.iterrows()))
gg.asserted('★ 注册的 kill ①:子集下可分性是否真的提高',r0.r_gl<0.8999-0.05,
            f"f=0.10 `corr(S_g,S_l)` = **{r0.r_gl:+.4f}** vs 全选项 `#298` 的 **+0.8999**")
gg.asserted('★ 注册的 kill ②:Ⓐ 是否仍明显大于 Ⓑ',
            (r0.v_a-r0.v_b)>2*r0.mde,
            f"Ⓐ **{r0.v_a:+.4f}** − Ⓑ **{r0.v_b:+.4f}** = {r0.v_a-r0.v_b:+.4f} vs 2×MDE {2*r0.mde:.4f}")
gg.asserted('⚠ 劈半两向的一致性(Ⓐ 的两向半差)',r0.a_spread<abs(r0.v_a),
            f"Ⓐ 两向半差 {r0.a_spread:.4f} vs |Ⓐ| {abs(r0.v_a):.4f}")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
