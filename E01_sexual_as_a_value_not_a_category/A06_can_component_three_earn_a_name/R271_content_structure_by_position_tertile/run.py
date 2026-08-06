import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A55 R271 -- C(递归)有没有一个正面证据

**换方向**(`#225` 的 NEXT:连着六轮方法轮,`§0.2` 的信号已亮)。
`#179` 杀掉了 B(位置是一个不带立场的坐标),而 **C 从没被正面检验过**。
**C 有一个具体的、A 与 B 都不预测的预测**:

    C 递归:赋值回流并**重塑表征** -> 在"位置"这一维上越强的人,**内容维度越有结构**
    A 专用内容系统:内容的结构由内容本身决定,与一个人在位置维上的位置无关
    B 对普通表征的赋值:位置是一个坐标,不回流 -> 同样无交互

ESTIMAND        按 `S`(位置分)分三层,各层各算一次**内容分半信度**(`#165` 同款:块劈半,
                各半算内容主成分分,跨人相关 + Spearman-Brown)。
KILL            **若信度随 S 单调上升且最高层明显高于最低层(差 > 2× 展布)-> C 得到第一个正面证据;
                若三层相同 -> C 在这份数据上没有支持。**
⚠ 最强混杂(跑之前写下)
                **`S` 与勾选数相关 +0.608**(`#100`)—— **勾得多的人信度天然更高**,
                那会单独造出 C 的图样。
控制(同一迭代内)
                ① 三层内**按勾选数卡钳匹配**(每层抽出勾选数分布相同的子样本),重算;
                ② 并排报**未匹配**与**匹配后**两条曲线。
NEGATIVE CTRL   每层内跨人置换 -> 该层自己的零信度。
POSITIVE CTRL   种一个"高 S 人群里更强"的内容信号 -> 梯度必须被测到。
IMPOSSIBLE      分层是按 `S` 切的,而 `S` 本身有测量误差 -> 层间有渗漏,梯度是**下界**。
"""
import numpy as np, pandas as pd, warnings, hashlib, zlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df); BL=[]
pos=np.zeros(NN); cnt=np.zeros(NN); K=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    BL.append(dict(M=M,ppl=ppl,rar=rr))
    pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); K[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
S=np.where(ok,pos/np.maximum(cnt,1),np.nan); K=np.where(ok,K,np.nan)
base=np.flatnonzero(np.isfinite(S)&np.isfinite(K))
print(f"n = {len(base):,};corr(S, 勾选数) = {np.corrcoef(S[base],K[base])[0,1]:+.4f}(`#100` 报 +0.608)")

def content_reliability(rows, half_seed, perm=False, plant=0.0, u=None):
    """块劈半,各半算内容主成分分,跨人相关 + Spearman-Brown。"""
    rg=np.random.default_rng(half_seed); p=rg.permutation(len(BL)); h=len(BL)//2
    sc=[np.zeros(NN),np.zeros(NN)]; ct=[np.zeros(NN),np.zeros(NN)]
    for side,cols in ((0,p[:h]),(1,p[h:2*h])):
        for j in cols:
            b=BL[j]; M=b['M'].copy(); ppl=b['ppl']
            sel=np.isin(ppl,rows)
            if sel.sum()<200: continue
            if perm:
                rg2=np.random.default_rng(zlib.crc32(f'{j}{half_seed}'.encode())%(1<<30))
                M=M[rg2.permutation(len(M))]
            if plant and u is not None:
                sub=(np.arange(M.shape[1])<max(2,M.shape[1]//3)).astype(float)
                M=M+plant*np.outer(u[ppl],sub)
            Msub=M[sel]; Z=Msub-Msub.mean(0,keepdims=True)
            w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
            sc[side][ppl[sel]]+=Z@v[:,-1]; ct[side][ppl[sel]]+=1
    a=np.where(ct[0]>=3,sc[0]/np.maximum(ct[0],1),np.nan)
    b_=np.where(ct[1]>=3,sc[1]/np.maximum(ct[1],1),np.nan)
    m=np.isfinite(a)&np.isfinite(b_)&np.isin(np.arange(NN),rows)
    if m.sum()<300: return np.nan,int(m.sum())
    r=abs(float(np.corrcoef(a[m],b_[m])[0,1]))
    return (2*r/(1+r) if r<0.999 else np.nan), int(m.sum())

q1,q2=np.nanpercentile(S[base],[33.3,66.7])
TERT={'低 S':base[S[base]<=q1],'中 S':base[(S[base]>q1)&(S[base]<=q2)],'高 S':base[S[base]>q2]}
print(f"三层 n = {[len(v) for v in TERT.values()]};各层勾选数中位 = "
      f"{[int(np.median(K[v])) for v in TERT.values()]}")

rng=np.random.default_rng(20260803)
def curve(tert, perm=False, plant=0.0, u=None, seeds=5):
    out={}
    for name,rows in tert.items():
        vs=[content_reliability(rows,900+s,perm,plant,u)[0] for s in range(seeds)]
        out[name]=(float(np.nanmean(vs)),float(np.nanstd(vs)))
    return out
raw=curve(TERT); nul=curve(TERT,perm=True)
print(f"\n{'层':<6}{'内容分半信度':>14}{'展布':>9}{'置换零':>10}{'净':>10}")
rows=[]
for name in TERT:
    net=raw[name][0]-nul[name][0]
    rows.append(dict(tertile=name,n=len(TERT[name]),rel=raw[name][0],sd=raw[name][1],
                     null=nul[name][0],net=net))
    print(f"{name:<6}{raw[name][0]:>14.4f}{raw[name][1]:>9.4f}{nul[name][0]:>10.4f}{net:>10.4f}")

# 勾选数卡钳匹配
tgt=np.percentile(K[TERT['低 S']],[10,30,50,70,90])
MT={}
for name,rows_ in TERT.items():
    sel=[]
    for lo,hi in zip([0]+list(tgt),list(tgt)+[1e9]):
        pool=rows_[(K[rows_]>=lo)&(K[rows_]<hi)]
        n_take=min(len(pool),int(0.2*min(len(v) for v in TERT.values())))
        if n_take>0: sel.append(rng.choice(pool,n_take,replace=False))
    MT[name]=np.concatenate(sel) if sel else rows_
print(f"\n勾选数匹配后 n = {[len(v) for v in MT.values()]};各层勾选数中位 = "
      f"{[int(np.median(K[v])) for v in MT.values()]}")
mraw=curve(MT); mnul=curve(MT,perm=True)
print(f"{'层':<6}{'匹配后信度':>14}{'置换零':>10}{'净':>10}")
for name in MT:
    rows[[r['tertile'] for r in rows].index(name)]['rel_matched']=mraw[name][0]
    rows[[r['tertile'] for r in rows].index(name)]['net_matched']=mraw[name][0]-mnul[name][0]
    print(f"{name:<6}{mraw[name][0]:>14.4f}{mnul[name][0]:>10.4f}{mraw[name][0]-mnul[name][0]:>10.4f}")

T=pd.DataFrame(rows); check_columns(T,'R271'); T.to_csv(pathlib.Path(__file__).parent/'results'/'tertiles.csv',index=False)
# 正对照:种一个"高 S 更强"的内容信号
u=np.zeros(NN); u[TERT['高 S']]=rng.standard_normal(len(TERT['高 S']))*1.0
u[TERT['中 S']]=rng.standard_normal(len(TERT['中 S']))*0.4
pl=curve(TERT,plant=0.6,u=u,seeds=3)
print(f"\n正对照(高 S 层种更强的内容信号):"
      + ' · '.join(f"{k} {v[0]:.4f}" for k,v in pl.items()))

lo_,hi_=T[T.tertile=='低 S'].iloc[0],T[T.tertile=='高 S'].iloc[0]
sd_pool=float(np.hypot(lo_.sd,hi_.sd))
mono=T.net_matched.iloc[0]<=T.net_matched.iloc[1]<=T.net_matched.iloc[2]
g=Gate('C 有没有正面证据')
g.asserted('正对照:高 S 层种更强内容信号 -> 梯度必须被测到',
           pl['高 S'][0]>pl['低 S'][0]+0.05,f"低 {pl['低 S'][0]:.4f} -> 高 {pl['高 S'][0]:.4f}")
g.asserted('⚠ 最强混杂已控制:勾选数卡钳匹配后三层勾选数中位相同',
           len(set(int(np.median(K[v])) for v in MT.values()))<=2,
           f"匹配前 {[int(np.median(K[v])) for v in TERT.values()]} -> 匹配后 {[int(np.median(K[v])) for v in MT.values()]}")
g.negative_control('各层置换零(高 S 层)',abs(float(hi_.null)),abs(float(hi_.rel)))
g.offset_control('高 S 层 vs 低 S 层(匹配后净信度)',float(hi_.net_matched),float(lo_.net_matched),sd_pool,
                 null_kind='同一条管道在低 S 层上的净信度 —— 不是零假设,是"若没有交互,高 S 层该落在哪"')
g.asserted('注册的 kill:信度随 S 单调上升且最高层明显高于最低层 -> C 得到正面证据',
           mono and (hi_.net_matched-lo_.net_matched>2*sd_pool),
           f"匹配后净信度 {T.net_matched.round(4).tolist()};单调={mono};"
           f"高−低 = {hi_.net_matched-lo_.net_matched:+.4f} vs 2×sd {2*sd_pool:.4f}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
