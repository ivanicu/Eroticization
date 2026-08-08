import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A65 R286 -- 把守卫 12 回溯扫一遍全账本,并当场补算最贵的那一格

**类型:CLOSURE**(如实标注)。产物是一张待办表 + 一格补算,不是一个新发现。

`#240b`:12 个守卫里只有 1 个防**假否定**。这一族的缺口可能不止 `#210` 一处。

⚠ **而在写扫描之前,我已经用推理找到一条,必须先说明它不是扫描的功劳**:
`#231`(`R276`)的头条 —— 普通半 −0.0618(n=5,322)vs 越轨半 +0.1007(n=6,113),差 +0.1625 ——
**两条臂的 n 不同**,因为各自要求 ≥4 块覆盖。**这正是 `#239a` 的形状,而它两轮前刚上了公开页。**

ESTIMAND        ① 扫描:账本里所有含撤回/降级语言**且**同条目出现 ≥2 个不同样本量的条目;
                ② 对最 load-bearing 的候选(`#231`)做**交集样本重比**。
KILL(逐项)      扫描:**`#210` 必须被扫出来**(已知阳性);
                含撤回语言但只有一个 n 的条目**必须不被**扫出来(已知阴性)。
                `#231`:**若交集样本上差值仍 > 2×展布 -> 结论不变,只是加一行作用域;
                若塌进展布 -> 公开页那条必须改写。**
IMPOSSIBLE      扫描是**词面**的:它找的是「撤回语言 + 两个 n」,不是「控制是否真的改变纳入」。
                所以它的产物是**候选**,每条仍需人读。假阴性(用了同一个 n 却换了人)扫不出来。
"""
import numpy as np, pandas as pd, warnings, re, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

txt=(ROOT/'RETRACTIONS.md').read_text()
ents=re.split(r'\n## Entry ',txt)[1:]
RET=('撤回','降级','UNVERIFIED','收窄')
rows=[]
for e in ents:
    num=re.match(r'(\d+)',e)
    if not num: continue
    n=int(num.group(1))
    # ⚠ 修:第一版只认 `n = 1,234`,而 `#210` 的样本量写在**表格单元**里(`| 9,944 |`)。
    #    正对照当场把它抓了出来 —— 那一版报「只有 1 条候选」,是 P5★ 的假无罪形状。
    ns=set()
    for m in re.finditer(r'n\s*=\s*\**([\d,]{3,})',e):
        v=int(m.group(1).replace(',',''))
        if v>=100: ns.add(v)
    for m in re.finditer(r'(?<![\d.])(\d{1,3},\d{3})(?![\d.])',e):
        v=int(m.group(1).replace(',',''))
        if 100<=v<=100000: ns.add(v)
    has=any(w in e for w in RET)
    rows.append(dict(entry=n,has_retraction=has,n_distinct=len(ns),
                     sizes=','.join(str(x) for x in sorted(ns)[:6])))
T=pd.DataFrame(rows)
cand=T[(T.has_retraction)&(T.n_distinct>=2)].copy()
readme=(ROOT/'README.md').read_text()
cand['on_readme']=[f"`#{int(r.entry)}`" in readme for _,r in cand.iterrows()]
print(f"账本条目 {len(T)};含撤回/降级语言 {int(T.has_retraction.sum())};"
      f"**候选(撤回语言 + ≥2 个不同 n)= {len(cand)}**;其中上了公开页的 **{int(cand.on_readme.sum())}**")
print("\n上了公开页的候选(按条目号倒序,越新越可能没被这一族检查过):")
for _,r in cand[cand.on_readme].sort_values('entry',ascending=False).head(12).iterrows():
    print(f"  #{int(r.entry):<4} n = {r.sizes}")
neg=T[(T.has_retraction)&(T.n_distinct<=1)]
print(f"\n已知阴性(有撤回语言但只有 ≤1 个 n):{len(neg)} 条,例如 "
      + ' · '.join(f"#{int(x)}" for x in neg.entry.tail(5)))
check_columns(cand,'R286'); cand.to_csv(pathlib.Path(__file__).parent/'results'/'candidates.csv',index=False)

# ---------- 当场补算 `#231` 的交集样本 ----------
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
rgm=np.random.default_rng(500); A=np.full((NB,NN),np.nan); B=np.full((NB,NN),np.nan)
for b,(M,ppl) in enumerate(MB):
    o=rgm.permutation(M.shape[1]); k=M.shape[1]//2
    A[b,ppl]=M[:,o[:k]].mean(1); B[b,ppl]=M[:,o[k:2*k]].mean(1)
def prof(X):
    F=np.isfinite(X); Z=np.where(F,X,0.0); tot=Z.sum(0); ct=F.sum(0); R=np.full_like(X,np.nan)
    for b in range(NB):
        lo=np.where(ct-F[b]>=6,(tot-Z[b])/np.maximum(ct-F[b],1),np.nan)
        R[b]=X[b]-lo; R[b]=R[b]-np.nanmean(R[b])
    return R
Ra,Rb=prof(A),prof(B)
def hs(R,cols):
    sub=R[cols]; F2=np.isfinite(sub)
    return np.where(F2.sum(0)>=4,np.nansum(np.where(F2,sub,0.0),0)/np.maximum(F2.sum(0),1),np.nan)
SC={nm:(hs(Ra,c)+hs(Rb,c))/2 for nm,c in (('普通',ORD),('越轨',TRG))}
SHN=[c for c in d.columns if str(c).lower().startswith('"i am ashamed')][0]
y=pd.to_numeric(d[SHN],errors='coerce').values.astype(float)
rng=np.random.default_rng(20260804)
def cr(x,mask=None):
    m=np.isfinite(x)&np.isfinite(y)&ok
    if mask is not None: m&=mask
    r=float(np.corrcoef(x[m],y[m])[0,1])
    sd=float(np.std([np.corrcoef(x[i],y[i])[0,1] for i in
        (rng.choice(np.flatnonzero(m),int(m.sum()),True) for _ in range(200))]))
    return r,sd,int(m.sum())
own={nm:cr(SC[nm]) for nm in SC}
common=np.isfinite(SC['普通'])&np.isfinite(SC['越轨'])&np.isfinite(y)&ok
com={nm:cr(SC[nm],common) for nm in SC}
g_own=own['越轨'][0]-own['普通'][0]; s_own=float(np.hypot(own['越轨'][1],own['普通'][1]))
g_com=com['越轨'][0]-com['普通'][0]; s_com=float(np.hypot(com['越轨'][1],com['普通'][1]))
print(f"\n`#231` 补算:")
print(f"  各自样本:普通 {own['普通'][0]:+.4f}(n={own['普通'][2]:,})· 越轨 {own['越轨'][0]:+.4f}"
      f"(n={own['越轨'][2]:,})· 差 **{g_own:+.4f}** vs 2×展布 {2*s_own:.4f}")
print(f"  **交集样本(n={com['普通'][2]:,}):普通 {com['普通'][0]:+.4f} · 越轨 {com['越轨'][0]:+.4f}"
      f" · 差 {g_com:+.4f} vs 2×展布 {2*s_com:.4f};保留 {100*g_com/g_own:.1f}%**")

g=Gate('回溯扫描 + `#231` 的交集重比')
g.asserted('⚠ 类型:CLOSURE —— 保护现存结论,产物是待办表 + 一格补算',True,'§0 三类动作')
g.asserted('扫描的已知阳性:`#210` 必须被扫出来',
           210 in set(cand.entry), f"候选 {len(cand)} 条;#210 {'在' if 210 in set(cand.entry) else '不在'}")
g.asserted('扫描的已知阴性:存在含撤回语言但只有 ≤1 个 n 的条目,且它们没被扫出来',
           len(neg)>0, f"{len(neg)} 条,例如 "+' · '.join(f"#{int(x)}" for x in neg.entry.tail(3)))
g.control_kept_the_sample('★ `#231` 的两条臂',before=own['普通'][0],after=own['越轨'][0],
                          n_before=own['普通'][2],n_after=own['越轨'][2],
                          before_common=com['普通'][0],after_common=com['越轨'][0],n_common=com['普通'][2])
g.asserted('★ 注册的 kill:`#231` 在交集样本上差值仍 > 2×展布 -> 结论不变,只加作用域',
           g_com>2*s_com, f"交集差 {g_com:+.4f} vs 2×展布 {2*s_com:.4f};保留 {100*g_com/g_own:.1f}%")
print(g)
print(f"\nsha1 {hashlib.sha1(cand.to_csv(index=False).encode()).hexdigest()[:12]}")
