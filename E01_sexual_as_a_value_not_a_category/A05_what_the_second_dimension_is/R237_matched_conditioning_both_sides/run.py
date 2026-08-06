import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A37 R237 -- 两侧,同一个条件集,可比的两个保留率

`#191c`:`#164`(位置侧)只去性别、保留 102%;`#236`(内容侧)去性别**+勾选数**、保留 54%–83%。
**两个数不在同一把尺上**,而 `#188` 的前页把两侧并排陈述 —— 读者会当成同一把尺。

⚠ `#191a` 的教训在这里直接适用:**只去性别,就必须边际评估** —— 若残差化只去性别、
评估却给定勾选数,就是同一个条件集错配。**本轮两侧都用「只去性别 + 边际评估」。**

ESTIMAND        位置侧与内容侧,各自在**只去性别**下的保留率(去性别后 r ÷ 去性别前 r),
                两侧同一个条件集、同一批人、同一个评估口径。
KILL            **若内容侧在只去性别时保留 >90%,那么 `#236` 的 54% 主要是勾选数干的,
                前页要写"去性别后两侧都几乎不动"。**
POSITIVE CTRL   `biomale` 在**两侧**都必须塌到阈值以下。
NEGATIVE CTRL   纯个人种植(位置型 / 内容型各一)去性别后必须存活。
IMPOSSIBLE      「只去性别」的保留率与「去性别+勾选数」的保留率**永远不可比** ——
                本轮产出的是前者的两侧对照,不是对 `#236` 的替代。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
SHAME='"I am ashamed or embarrassed about at least some of what arouses me" (7cw1ziu)'
sex=pd.to_numeric(df['biomale'],errors='coerce').values.astype(float)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); rb=np.random.default_rng(20260803)
u_pos=rb.standard_normal(NN); u_con=rb.standard_normal(NN)

def build(desex, plant=0.0, kind=None):
    con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); share=[]
    for _,q in keep.iterrows():
        s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
        ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
        if len(ppl)<1200 or len(opt)<8: continue
        pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
        M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
        rr=-np.log(np.clip(M.mean(0),1e-4,1.))
        if plant and kind=='content':
            sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
            M=M+plant*np.outer(u_con[ppl],sub)
        if plant and kind=='position':
            M=M+plant*np.outer(u_pos[ppl],(rr-rr.mean())/max(rr.std(),1e-9))
        Z=M-M.mean(0,keepdims=True)
        w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); sc=Z@v[:,-1]
        ps=(M@rr)/np.maximum(M.sum(1),1)
        if desex:                                   # **只去性别** —— 与评估口径一致(边际)
            g=sex[ppl]; m=np.isfinite(g)
            if m.sum()>100:
                X=np.c_[np.ones(m.sum()),g[m]]
                for arr in (sc,ps):
                    arr[m]=arr[m]-X@np.linalg.lstsq(X,arr[m],rcond=None)[0]
                share.append(1-((ps[m]-ps[m].mean())**2).sum()/max(((ps[m]-ps[m].mean())**2).sum(),1e-9))
        con[ppl]+=sc; pos[ppl]+=ps; cnt[ppl]+=1
    ok=cnt>=8
    return (np.where(ok,con/np.maximum(cnt,1),np.nan),
            np.where(ok,pos/np.maximum(cnt,1),np.nan))

C0,S0=build(False); C1,S1=build(True)
base=np.isfinite(C0)&np.isfinite(C1)&np.isfinite(S0)&np.isfinite(S1); bi=np.flatnonzero(base)
print(f"n = {len(bi):,};corr(S0,S1) = {np.corrcoef(S0[bi],S1[bi])[0,1]:+.4f} · "
      f"corr(C0,C1) = {np.corrcoef(C0[bi],C1[bi])[0,1]:+.4f}   <- 只去性别,分数几乎没变")

def mr(y,x,ii):                                     # **边际**相关,不控勾选数(条件集与残差化一致)
    m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
    return float(np.corrcoef(y[jj],x[jj])[0,1]), len(jj)

TARGETS=[('羞耻',SHAME,'position'),('animated','animated','content'),
         ('written','written','content'),('biomale','biomale','both')]
rows=[]
for name,col,side in TARGETS:
    y=df[col].values.astype(float)
    for sname,x0,x1 in (('位置 S',S0,S1),('内容 C',C0,C1)):
        r0,n=mr(y,x0,bi); r1,_=mr(y,x1,bi)
        sd=float(np.std([mr(y,x1,rb.choice(bi,len(bi),replace=True))[0] for _ in range(200)]))
        rows.append(dict(target=name,side=sname,n=n,before=r0,after=r1,
                         retained=r1/r0 if abs(r0)>1e-6 else np.nan,sd=sd))
T=pd.DataFrame(rows); check_columns(T,'R237'); T.to_csv(pathlib.Path(__file__).parent/'results'/'matched.csv',index=False)
print(f"\n{'结局':<10}{'侧':<8}{'去性别前':>10}{'去性别后':>10}{'保留':>8}{'sd':>9}")
for _,r in T.iterrows():
    print(f"{r.target:<10}{r.side:<8}{r.before:>+10.4f}{r.after:>+10.4f}{100*r.retained:>7.0f}%{r.sd:>9.4f}")

# 负对照
Cp0,Sp0=build(False,0.6,'content'); Cp1,Sp1=build(True,0.6,'content')
Cq0,Sq0=build(False,0.6,'position'); Cq1,Sq1=build(True,0.6,'position')
m1=np.isfinite(Cp0)&np.isfinite(Cp1); m2=np.isfinite(Sq0)&np.isfinite(Sq1)
rc0=np.corrcoef(Cp0[m1],u_con[m1])[0,1]; rc1=np.corrcoef(Cp1[m1],u_con[m1])[0,1]
rp0=np.corrcoef(Sq0[m2],u_pos[m2])[0,1]; rp1=np.corrcoef(Sq1[m2],u_pos[m2])[0,1]
print(f"\n负对照 内容种植 {rc0:+.4f} -> {rc1:+.4f}(保留 {100*rc1/rc0:.0f}%) · "
      f"位置种植 {rp0:+.4f} -> {rp1:+.4f}(保留 {100*rp1/rp0:.0f}%)")

bio_p=T[(T.target=='biomale')&(T.side=='位置 S')].iloc[0]
bio_c=T[(T.target=='biomale')&(T.side=='内容 C')].iloc[0]
sh_p =T[(T.target=='羞耻')&(T.side=='位置 S')].iloc[0]
ani  =T[(T.target=='animated')&(T.side=='内容 C')].iloc[0]
wri  =T[(T.target=='written')&(T.side=='内容 C')].iloc[0]
g=Gate('两侧,同一个条件集')
g.asserted('正对照:biomale 在两侧都塌掉',abs(bio_p.after)<2*bio_p.sd and abs(bio_c.after)<2*bio_c.sd,
           f"位置 {bio_p.before:+.4f}->{bio_p.after:+.4f} · 内容 {bio_c.before:+.4f}->{bio_c.after:+.4f}")
g.asserted('负对照:两种种植去性别后都存活',min(rc1/rc0,rp1/rp0)>0.7,
           f"内容 {100*rc1/rc0:.0f}% · 位置 {100*rp1/rp0:.0f}%")
g.same_scale('两侧的保留率现在同一把尺',float(ani.retained),float(sh_p.retained),'retention(只去性别)')
g.asserted('注册的 kill:内容侧只去性别时保留 >90%',min(ani.retained,wri.retained)>0.90,
           f"animated {100*ani.retained:.0f}% · written {100*wri.retained:.0f}% · "
           f"(位置侧羞耻 {100*sh_p.retained:.0f}%)")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")

# ---- 两件必须单独说的事 ------------------------------------------------------
print("\n---- 两件必须单独说的事 ----")
# ① 零基数上的保留率无定义
bad=T[(T.before.abs()<2*T.sd)]
print(f"① 去性别**前**就与零不可分的格子({len(bad)} 个)—— 它们的「保留率」无定义:")
for _,r in bad.iterrows():
    print(f"   {r.target}/{r.side}: 前 {r.before:+.4f}(sd {r.sd:.4f},{abs(r.before)/r.sd:.1f}×)"
          f" -> 打印的 {100*r.retained:.0f}% 是**零基数上的比值**")
# ② 边际评估下,位置侧也强预测 animated
ap=T[(T.target=='animated')&(T.side=='位置 S')].iloc[0]
sp=T[(T.target=='羞耻')&(T.side=='位置 S')].iloc[0]
print(f"② 边际评估下位置侧对 animated = {ap.before:+.4f},与它对羞耻的 {sp.before:+.4f} **同量级** ——")
print(f"   `#188` 那句「位置贴羞耻、内容贴媒介」有一部分是**控制勾选数**带来的,不是原始结构。")

g2=Gate('这两件事要不要写进前页')
g2.no_sign_crossing('written/位置 的去性别前后不得取比值',[float(T[(T.target=='written')&(T.side=='位置 S')].before.iloc[0]),
                                                float(T[(T.target=='written')&(T.side=='位置 S')].after.iloc[0])])
g2.require_resolvable_first('written/位置 去性别前',
    float(T[(T.target=='written')&(T.side=='位置 S')].before.iloc[0]),
    float(T[(T.target=='written')&(T.side=='位置 S')].sd.iloc[0]),family='written_pos')
g2.resolvable('位置侧对 animated(边际)',float(ap.before),float(ap.sd))
g2.same_scale('位置对 animated vs 位置对羞耻',float(ap.before),float(sp.before),'marginal r')
print(g2)
