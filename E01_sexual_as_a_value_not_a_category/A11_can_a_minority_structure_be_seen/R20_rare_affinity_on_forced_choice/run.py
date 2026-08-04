import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A11 R20 -- 稀有偏好特质,搬到一个「什么都说是」无处施力的仪器上。

#125 的 NEXT 指向跨题变异,但设计时看到一个更锋利、也更直接可解释的版本,而且它同时回答
本项目最发达那个声明的最大弱点:

  #95/#99/#100 立起一条「稀有选项亲和」的人格维度(信度 +0.432,地板 -0.022)。
  它全部建立在 0-5 多选的勾选上 —— 而那正是「什么都说是」有施力点的地方。
  #100c 论证过它不是勾选数(移除后存活 67%),但那是**同一族仪器内部**的论证。

强制单选按构造消除作答水平(#123b)。于是:

  问 (a) 口味广的人,被迫只选一个时,选的是不是更冷门的那个?
  问 (b) #95 的稀有亲和特质,能不能预测他在**强制单选**里选多冷门的?
         -> 若能,这条特质第一次拿到**跨仪器**验证,而且是在作答风格无处施力的仪器上
         -> 若不能,它可能一直是同一族仪器内部的东西

⚠ 跑之前写下的混淆(#96a:必须断言在代码里,不能写在散文里):
  (a) 性别·年龄 -> 协变量
  (b) 答了多少个强制单选题 -> 协变量
  (c) **S 与广度相关**(#100c 测得 +0.608)-> 两个问题必须**互相控制**,否则 (b) 只是 (a) 的回声
  (d) 每题的选项流行度用**该题自己的**分布算,不跨题

ESTIMAND        选中选项的意外度 -log(流行度) 对 (a) 广度 (b) 稀有亲和 S 的偏相关,
                两者互相控制,逐题算再合并。
KILL            threshold-free;零的种类必须命名;gate 顺序按 #120d;
                negative_control 传 null_spread(#125)。
POSITIVE CTRL   把一部分人的选择替换成「该题最冷门的选项」,强度按其 S 排序 -> 必须单调,
                g=0 用 degenerate_matches_reference 精确复现真实臂(#124f)。
NEGATIVE CTRL   S 与广度在(性别 x 年龄)分层内**联合**打乱(保住两者的相关)。
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
# #95 的稀有亲和 S:每人勾选项的平均 -log(基率),用多选块算
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
tot=np.zeros(len(df)); cnt=np.zeros(len(df))
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    ref=-np.log(np.clip(M.mean(0),1e-4,1.))
    tot[ppl]+=M@ref; cnt[ppl]+=M.sum(1)
ok=cnt>=15
S=np.full(len(df),np.nan); S[ok]=tot[ok]/cnt[ok]
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
FC=inv[inv['kind']=='FORCED_CHOICE_MOST']['col'].tolist()
nfc=df[FC].notna().sum(axis=1).values
zs=lambda v:(v-np.nanmean(v))/(np.nanstd(v)+1e-9)
print(f"稀有亲和 S 可算的人 {int(ok.sum()):,}   corr(S, 广度) = "
      f"{np.corrcoef(S[ok],breadth[ok])[0,1]:+.3f}  (#100c 测得 +0.608)",flush=True)
def pcorr(y,x,Z):
    """x 对 y 的偏相关,控制 Z。"""
    D=np.c_[np.ones(len(y)),Z]
    ry=y-D@np.linalg.lstsq(D,y,rcond=None)[0]
    rx=x-D@np.linalg.lstsq(D,x,rcond=None)[0]
    return float(np.corrcoef(ry,rx)[0,1])
rows=[]
for c in FC:
    s=df[c]; m=s.notna()&np.isfinite(S)&np.isfinite(breadth)
    idx=np.flatnonzero(m)
    if len(idx)<800: continue
    v=s.iloc[idx]
    pop=v.map(v.value_counts()/len(v))          # (d) 用该题自己的分布
    y=-np.log(np.clip(pop.values.astype(float),1e-4,1.))
    b=breadth[idx]; sa=S[idx]
    Zb=np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(sa)]   # (c) 互相控制
    Zs=np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(b)]
    rb=pcorr(y,zs(b),Zb); rs=pcorr(y,zs(sa),Zs)
    # 正对照:按 S 排序,把最高的一部分人的选择换成该题最冷门的选项
    rare_y=y.max()
    ps={}
    for g in [0.0,0.05,0.15]:
        yp=y.copy()
        if g>0:
            hi=np.argsort(-sa)[:int(g*len(idx))]
            yp[hi]=rare_y
        ps[g]=pcorr(yp,zs(sa),Zs)
    # 负对照:S 与广度**联合**打乱(保住两者相关)
    rp=np.random.default_rng(31); st=(male[idx]>0).astype(int)*5+agev[idx].astype(int)
    bp=b.copy(); sp=sa.copy()
    for q in np.unique(st):
        w=np.flatnonzero(st==q)
        if len(w)>1:
            pm=rp.permutation(len(w)); bp[w]=b[w][pm]; sp[w]=sa[w][pm]
    nb=pcorr(y,zs(bp),np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(sp)])
    ns=pcorr(y,zs(sp),np.c_[zs(male[idx]),zs(agev[idx]),zs(nfc[idx]),zs(bp)])
    rows.append(dict(v_col=c[:22],n=len(idx),k=v.nunique(),r_breadth=rb,r_S=rs,
                     n_breadth=nb,n_S=ns,p0=ps[0.0],p05=ps[0.05],p15=ps[0.15]))
    print(f"  {c[:20]:<20} n={len(idx):5,} 选项{v.nunique():3d}  广度 {rb:+.4f}  S {rs:+.4f}  "
          f"零 {nb:+.4f}/{ns:+.4f}  种植 {ps[0.05]:+.4f}/{ps[0.15]:+.4f}",flush=True)
check_coverage(len(rows),len(FC),'A11R20',tol=0.25)
print(f"  纳入 {len(rows)}/{len(FC)}")
D=check_columns(pd.DataFrame(rows),'A11R20')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
print(f"\n=== 选中选项的冷门程度,对广度 / 对稀有亲和 S(互相控制) ===")
print(D.round(4).to_string(index=False))
for nm,e,n in [('广度',D.r_breadth.values,D.n_breadth.values),('稀有亲和 S',D.r_S.values,D.n_S.values)]:
    print(f"\n  {nm:10s} 效应 {e.mean():+.4f} ± {se(e):.4f} ({abs(e.mean())/se(e):.1f}x)   "
          f"零 {n.mean():+.4f} ± {se(n):.4f}   逐题同向 {(np.sign(e)==np.sign(e.mean())).sum()}/{len(e)}")
g=Gate("稀有亲和特质,在一个作答风格无处施力的仪器上还在吗?")
g.degenerate_matches_reference("g=0 精确复现真实臂 (#124f)", degenerate=D.p0.mean(), reference=D.r_S.mean())
g.no_sign_crossing("种植阶梯(相对 g=0)", [D.p05.mean()-D.p0.mean(), D.p15.mean()-D.p0.mean()])
g.positive_control("种植 0.15 被测出", planted=D.p15.mean()-D.p0.mean(), floor=0.0, spread=se(D.r_S.values))
g.require_resolvable_first("S -> 选中冷门程度", effect=D.r_S.mean(), spread=se(D.r_S.values))
g.negative_control("S 与广度联合打乱", null=D.n_S.mean(), effect=D.r_S.mean(), null_spread=se(D.n_S.values))
print(); print(g)
gb=Gate("口味广的人,被迫只选一个时选得更冷门吗?")
gb.require_resolvable_first("广度 -> 选中冷门程度", effect=D.r_breadth.mean(), spread=se(D.r_breadth.values))
gb.negative_control("联合打乱", null=D.n_breadth.mean(), effect=D.r_breadth.mean(),
                    null_spread=se(D.n_breadth.values))
print(); print(gb)
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
