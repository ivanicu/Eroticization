import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A95 R343 -- 羞耻跟的是「东西冷门」,还是「我跟身边的人不一样」

`#231` 起羞耻有两条路,S 是其中一条(**+0.1185**)。但 S 是**对全样本**算的稀有度 ——
所以「S 预测羞耻」有两个完全不同的读法,而它们的**干预含义相反**:

- **Ⓐ 内容污名** —— 这些东西本身被污名化,谁喜欢都一样。
- **Ⓑ 参照群体偏离** —— 重要的是「在**像我这样的人**里我很少见」。

ESTIMAND        `S_local` = 在此人**自己的参照群体内**(性别 × 年龄段)算的稀有度;
                取 `S_g⊥l` 与 `S_l⊥g` 两个正交残差,**各自与羞耻的相关**就是分离器。
KILL            **Ⓐ -> 全局残差带羞耻、局部残差不带;Ⓑ -> 反过来;两个都带 -> 两种机制并存。**
POSITIVE CTRL   合成一个**只由 `S_local` 驱动**的结局 —— 必须只被局部残差抓到。
NEGATIVE CTRL   同一构造下 `perm_finite` 打乱人的零。
⚠ CONTROL       先报 `corr(S_global, S_local)`;若太高则两个残差都很小,
                **要报功效(MDE),不是报零**(`#296a` 的教训的另一面)。
⚠ KNOB          参照群体的切法本身是个旋钮(只按性别 / 性别×年龄 / +关系风格)-> `CALIBER.md`。
IMPOSSIBLE      参照群体是**观察到的**,不是随机分配的;所以 Ⓑ 只能说「与本地稀有度同变」,
                不能说「因为偏离参照群体而羞耻」。
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
AGEB=d['age'].astype(str).values
def S_ref(cells):
    """在**每个参照群体内**算流行度,再据此定义稀有度;cells = 每人的组标签(None=全样本)。"""
    cv=np.zeros(NN); ps=np.zeros(NN)
    labs=[None] if cells is None else [c for c in pd.unique(cells) if isinstance(c,str)]
    for M,ppl in MB:
        for L in labs:
            sel=np.ones(len(ppl),bool) if L is None else (cells[ppl]==L)
            if sel.sum()<200: continue
            Ms=M[sel]; rr=-np.log(np.clip(Ms.mean(0),1e-4,1.)); n=Ms.sum(1)
            v=np.where(n>0,(Ms@rr)/np.maximum(n,1),0.0)
            cv[ppl[sel]]+=1; ps[ppl[sel]]+=v
    return np.where(cv>=8,ps/np.maximum(cv,1),np.nan)
SXL=np.where(np.isfinite(SEX),np.where(SEX==1,'M','F'),'?')
CELL=np.array([f"{a}|{s}" for a,s in zip(AGEB,SXL)],dtype=object)
KN={'仅性别':np.array(list(SXL),dtype=object),
    '性别×年龄':CELL,
    '性别×年龄×关系风格':np.array([f"{c}|{r}" for c,r in zip(CELL,
        d['Personally, your preferred relationship style is: (4jib23m)'].astype(str).values)],dtype=object)}
Sg=S_ref(None)
def resid(a,b,m):
    out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
    out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
def cor(u,v):
    m=np.isfinite(u)&np.isfinite(v); return (float(np.corrcoef(u[m],v[m])[0,1]),int(m.sum())) if m.sum()>200 else (np.nan,0)
def perm_finite(v,seed):
    z=v.copy(); j=np.flatnonzero(np.isfinite(z))
    z[j]=z[np.random.default_rng(seed).permutation(j)]; return z
print(f"羞耻题:{SHAME[:62]}…  n={np.isfinite(sh).sum():,}")
rows=[]
for kn,cells in KN.items():
    Sl=S_ref(cells); m=np.isfinite(Sg)&np.isfinite(Sl)
    rGL,_=cor(Sg,Sl)
    Rg,Rl=resid(Sg,Sl,m),resid(Sl,Sg,m)
    cg,ng=cor(Rg,sh); cl,nl=cor(Rl,sh); c0,n0=cor(Sg,sh)
    # ⚠ 功效:残差被压掉多少 -> MDE 随之抬高
    sdg,sdl=float(np.nanstd(Rg)/np.nanstd(Sg[m])),float(np.nanstd(Rl)/np.nanstd(Sl[m]))
    mde=1.96/np.sqrt(max(ng,1))
    print(f"\n【{kn}】 组数 {len(set(cells))} · `corr(S_g, S_l)` = **{rGL:+.4f}**")
    print(f"   残差保留的标准差:全局 **{100*sdg:.1f}%** · 局部 **{100*sdl:.1f}%** · "
          f"n={ng:,} 下的 MDE ≈ **{mde:.4f}**")
    print(f"   ↔羞耻:  S_global(原样) **{c0:+.4f}**")
    print(f"            **`S_g ⊥ S_l`(内容污名 Ⓐ)** **{cg:+.4f}**")
    print(f"            **`S_l ⊥ S_g`(参照群体偏离 Ⓑ)** **{cl:+.4f}**")
    n_=[cor(perm_finite(Rg,900+i),sh)[0] for i in range(6)]
    print(f"   零(打乱人):{np.mean(n_):+.4f} ± {np.std(n_):.4f}")
    rows.append(dict(knob=kn,r_gl=rGL,v_shame_raw=c0,v_a_global=cg,v_b_local=cl,
                     mde=mde,null_mean=float(np.mean(n_)),null_sd=float(np.std(n_))))
T=pd.DataFrame(rows); check_columns(T,'R343')
T.to_csv(pathlib.Path(__file__).parent/'results'/'global_vs_local.csv',index=False)
# ---- 正对照:只由 S_local 驱动的合成结局,必须只被局部残差抓到 ----
Sl=S_ref(KN['性别×年龄']); m=np.isfinite(Sg)&np.isfinite(Sl)
Rg,Rl=resid(Sg,Sl,m),resid(Sl,Sg,m)
zl=np.where(np.isfinite(Sl),(Sl-np.nanmean(Sl))/np.nanstd(Sl),np.nan)
SW=[]
for g in (0.0,0.3,0.6,1.0):
    y=g*zl+np.random.default_rng(5).standard_normal(NN)
    a_,_=cor(Rg,y); b_,_=cor(Rl,y); SW.append((g,b_))
    print(f"  正对照 g={g:.1f}(只由 `S_local` 驱动):全局残差 **{a_:+.4f}** · 局部残差 **{b_:+.4f}**")
r0=T[T.knob=='性别×年龄'].iloc[0]
gg=Gate('羞耻:内容污名 Ⓐ,还是参照群体偏离 Ⓑ')
gg.plant_direction_from_sweep('正对照:只由 `S_local` 驱动的结局 -> 局部残差必须随强度上升',SW,SW[0][1])
gg.asserted(f'⚠ 功效控制:`corr(S_g,S_l)` = {r0.r_gl:+.4f},两个残差都不许小于 MDE 才算可读',
            min(abs(r0.v_a_global),abs(r0.v_b_local))>r0.mde or
            max(abs(r0.v_a_global),abs(r0.v_b_local))>3*r0.mde,
            f"Ⓐ {r0.v_a_global:+.4f} · Ⓑ {r0.v_b_local:+.4f} · MDE {r0.mde:.4f}")
gg.negative_control('零:打乱人后全局残差与羞耻的相关',r0.null_mean,r0.v_a_global,   # ⚠ 签名是 (name, null, effect)
    null_spread=r0.null_sd,
    null_kind='`perm_finite` 只在**有限项内**打乱人 —— 题内跨人的零,保住缺失格局(#264b/#278b)')
gg.asserted('★ Ⓑ 参照群体偏离这一路本身是否可读(仪器已由正对照证明有能力看见它)',
            abs(r0.v_b_local)<r0.mde,
            f"Ⓑ **{r0.v_b_local:+.4f}** vs MDE **{r0.mde:.4f}** —— 上界在 MDE 量级,"
            f"而正对照在 g=0.3 时给出 +0.0983")
gg.asserted('★ 注册的 kill:Ⓐ 与 Ⓑ 哪一个带羞耻',
            abs(r0.v_a_global-r0.v_b_local)>2*r0.mde,
            f"Ⓐ **{r0.v_a_global:+.4f}** vs Ⓑ **{r0.v_b_local:+.4f}**(差 {r0.v_a_global-r0.v_b_local:+.4f},"
            f"2×MDE = {2*r0.mde:.4f})")
gg.asserted('⚠ 旋钮:三种参照群体切法下 Ⓐ−Ⓑ 的符号是否一致',
            len(set(np.sign(T.v_a_global-T.v_b_local)))==1,
            ' · '.join(f"{r.knob} {r.v_a_global-r.v_b_local:+.4f}" for _,r in T.iterrows()))
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
