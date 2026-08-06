import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R21 -- 留出块重算 S。R20 的结果必须先过这一关。

R20 得到:稀有亲和 S -> 选中选项的冷门程度 +0.1822(7.9x,10/10 同向)。
**但我跑之前没写下最强的那个混淆**,跑完才查:强制单选的选项与多选块的选项重叠
——10 个题里 7 个重叠 89-100%。对那 7 个,S 与结局共享同一批 item,是恒等式不是测量。

3 个题重叠为 0(youfeelmost、otherfeel1most、pregnancy2most),它们的效应是
+0.136 / +0.123 / +0.205 —— 干净且为正。但 n=3。

正确修法:**对每个强制单选题,把与它重叠的多选块从 S 里剔掉再算**(留一块 S)。
这样 10 个题全部干净,而不是只剩 3 个。

ESTIMAND        选中选项冷门程度 对 留出块 S 的偏相关(控制广度·性别·年龄·答题数)。
KILL            threshold-free;gate 顺序按 #120d;negative_control 传 null_spread(#125)。
POSITIVE CTRL   按 S 排序把最高的一部分人的选择换成最冷门项;g=0 用
                degenerate_matches_reference 精确复现(#124f)。
NEGATIVE CTRL   S 与广度在(性别 x 年龄)分层内**联合**打乱。
额外报告        逐题重叠率,以及"重叠为 0 的三题"与"其余七题"分开的结果 —— 若两组一致,
                重叠不是驱动因素;若只有重叠组有效应,那才是恒等式。
IMPOSSIBLE      "真的更爱冷门"与"在问卷上更愿意挑冷门"分不开。
"""
import pandas as pd, numpy as np, warnings, hashlib
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns, check_coverage
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
P=(df[rate].apply(pd.to_numeric,errors='coerce')>0).astype(float).fillna(0.).values
breadth=P.sum(1)
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
BLK={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    BLK[q.qi]=dict(ppl=ppl,M=M,ref=-np.log(np.clip(M.mean(0),1e-4,1.)),opts=set(opt))
TOT=np.zeros(len(df)); CNT=np.zeros(len(df))
per={}
for k,b in BLK.items():
    t=b['M']@b['ref']; c=b['M'].sum(1)
    per[k]=(b['ppl'],t,c); TOT[b['ppl']]+=t; CNT[b['ppl']]+=c
def S_without(drop):
    tot=TOT.copy(); cnt=CNT.copy()
    for k in drop:
        p,t,c=per[k]; tot[p]-=t; cnt[p]-=c
    out=np.full(len(df),np.nan); ok=cnt>=15
    out[ok]=tot[ok]/cnt[ok]; return out
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
FC=inv[inv['kind']=='FORCED_CHOICE_MOST']['col'].tolist()
nfc=df[FC].notna().sum(axis=1).values
zs=lambda v:(v-np.nanmean(v))/(np.nanstd(v)+1e-9)
def pcorr(y,x,Z):
    D=np.c_[np.ones(len(y)),Z]
    ry=y-D@np.linalg.lstsq(D,y,rcond=None)[0]; rx=x-D@np.linalg.lstsq(D,x,rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])
rows=[]
for c in FC:
    v0=df[c].dropna().astype(str); o=set(v0.unique())
    drop=[k for k,b in BLK.items() if len(o&b['opts'])>=0.5*len(o)]     # 与该题重叠过半的块全剔
    ovl=max([len(o&b['opts'])/len(o) for b in BLK.values()],default=0.)
    Sv=S_without(drop)
    s=df[c]; m=s.notna()&np.isfinite(Sv)&np.isfinite(breadth); idx=np.flatnonzero(m)
    if len(idx)<800: continue
    v=s.iloc[idx]; pop=v.map(v.value_counts()/len(v))
    y=-np.log(np.clip(pop.values.astype(float),1e-4,1.))
    b_=breadth[idx]; sa=Sv[idx]
    Zs=np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(b_)]
    rs=pcorr(y,zs(sa),Zs)
    ps={}
    for g in [0.0,0.05,0.15]:
        yp=y.copy()
        if g>0: yp[np.argsort(-sa)[:int(g*len(idx))]]=y.max()
        ps[g]=pcorr(yp,zs(sa),Zs)
    rp=np.random.default_rng(31); st=(male[idx]>0).astype(int)*5+agev[idx].astype(int)
    bp=b_.copy(); sp=sa.copy()
    for q in np.unique(st):
        w=np.flatnonzero(st==q)
        if len(w)>1:
            pm=rp.permutation(len(w)); bp[w]=b_[w][pm]; sp[w]=sa[w][pm]
    ns=pcorr(y,zs(sp),np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(bp)])
    rows.append(dict(v_col=c[:22],n=len(idx),ovl=ovl,ndrop=len(drop),r_S=rs,n_S=ns,
                     p0=ps[0.0],p05=ps[0.05],p15=ps[0.15]))
    print(f"  {c[:20]:<20} n={len(idx):5,} 重叠{ovl:4.0%} 剔块{len(drop)}  S {rs:+.4f}  零 {ns:+.4f}",flush=True)
check_coverage(len(rows),len(FC),'A11R21',tol=0.25)
D=check_columns(pd.DataFrame(rows),'A11R21')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
E=D.r_S.values; N=D.n_S.values
print(f"\n=== 留出块 S -> 选中选项冷门程度(R20 未留出时是 +0.1822) ===")
print(D.round(4).to_string(index=False))
print(f"\n  效应 {E.mean():+.4f} ± {se(E):.4f} ({abs(E.mean())/se(E):.1f}x)  零 {N.mean():+.4f}  "
      f"逐题同向 {(E>0).sum()}/{len(E)}")
lo=D[D.ovl<0.5]; hi=D[D.ovl>=0.5]
print(f"\n  重叠<50% 的 {len(lo)} 题:{lo.r_S.mean():+.4f}     重叠>=50% 的 {len(hi)} 题:{hi.r_S.mean():+.4f}")
print(f"  -> 若两组一致,重叠不是驱动因素")
g=Gate("留出块之后,稀有亲和还能预测强制单选吗?")
g.degenerate_matches_reference("g=0 精确复现 (#124f)", degenerate=D.p0.mean(), reference=E.mean())
g.no_sign_crossing("种植阶梯", [D.p05.mean()-D.p0.mean(), D.p15.mean()-D.p0.mean()])
g.positive_control("种植 0.15 被测出", planted=D.p15.mean()-D.p0.mean(), floor=0.0, spread=se(E))
g.require_resolvable_first("留出块 S 的效应", effect=E.mean(), spread=se(E))
g.negative_control("S 与广度联合打乱", null=N.mean(), effect=E.mean(), null_spread=se(N))
g.asserted("重叠组与非重叠组一致(重叠不是驱动因素)",
           abs(lo.r_S.mean()-hi.r_S.mean())<2*max(se(lo.r_S.values),se(hi.r_S.values)),
           f"{lo.r_S.mean():+.4f} vs {hi.r_S.mean():+.4f}")
print(); print(g)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
