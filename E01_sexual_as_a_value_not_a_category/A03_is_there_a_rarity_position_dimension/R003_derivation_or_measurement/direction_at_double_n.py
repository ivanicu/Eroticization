import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R14 -- 加强 #117c(2.1x),并补上它根本没有的正对照。

#117 的 NEXT:方向命题只有 2.1x,最便宜的杠杆是 n —— 配对门槛从 >=400 降到 >=250,对数约翻倍,
而统计量是逐对均值,直接受益。**预先承诺(#111c):双倍 n 下若仍不过 3x,就按 2.1x 引用然后走开,
不追第三轮。**

而且复查 #117c 时发现一个更严重的问题:**它根本没有正对照。** 它只过了负对照(分层置换)和
可分辨性。P5 星号规则:一个从未返回过非零的仪器给出的零(或任何读数)是沉默,不是测量。
这一轮补上,而且是分级的。

统计量(与 #117 相同,以便可比):
  w = 「更喜欢 A 的人」减「更喜欢 B 的人」的残差化偏好轮廓差,单位化
  位移 = A 先组在 w 上的投影均值 − B 先组的
  残差化掉:这两项各自的评分、它们的差、以及全部人层协变量(含人均评分)

ESTIMAND        逐对位移的均值,配对门槛 >=250。
KILL            threshold-free。预先承诺:>=3x 才算加强,否则按 2.1x 引用并停止。
POSITIVE CTRL   分级种植:把 A 先组沿 w 真的推 g 个标准差(g = 0, 0.05, 0.15)。
                必须单调,且 g=0 必须不开火。
NEGATIVE CTRL   顺序标签在协变量分层内打乱。
IMPOSSIBLE      方向不判别(#113c 已查)—— 本轮仍标 DESCRIPTION,只加强强度不改变性质。
"""
import pandas as pd, numpy as np, warnings, hashlib, itertools, re
sys.path.insert(0,str(ROOT))
from lib.gates import Gate, check_columns
warnings.filterwarnings('ignore')
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False); inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
V=pd.DataFrame({c:df[c].map(BIN) for c in ons}).values
CONCRETE={1,2,4,13,14,16,17,18,20,26,29,30}; RELATIONAL={3,5,6,8,9,10,11,12,15,21,22,23,24,27}
KIND={**{i:'C' for i in CONCRETE},**{i:'R' for i in RELATIONAL}}
rate=[c for c in inv[inv['kind']=='RATING_0_5']['col'] if c in df.columns]
R=df[rate].apply(pd.to_numeric,errors='coerce').values
P=(np.nan_to_num(R)>0).astype(float); breadth=P.sum(1)
meanrating=np.nanmean(np.where(np.isfinite(R),R,np.nan),axis=1)
meanrating=np.where(np.isfinite(meanrating),meanrating,np.nanmean(meanrating))
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}
male=pd.to_numeric(df.get('biomale'),errors='coerce').fillna(-1).values
agev=df['age'].map(AGEMAP).fillna(2).values
prec=np.nanmean(V,axis=1); prec=np.where(np.isfinite(prec),prec,np.nanmean(prec))
mc=np.nanmean(V[:,sorted(CONCRETE)],axis=1); mr=np.nanmean(V[:,sorted(RELATIONAL)],axis=1)
mc=np.where(np.isfinite(mc),mc,np.nanmean(mc)); mr=np.where(np.isfinite(mr),mr,np.nanmean(mr))
zs=lambda X:(X-X.mean(0))/(X.std(0)+1e-9)
BASE=zs(np.c_[male,agev,breadth,prec,mc,mr,mc-mr,meanrating])
def norm(s): return re.sub(r'[^a-z]',' ',s.lower())
best={}
for j,c in enumerate(ons):
    m=re.search(r'interest in ([a-z /-]+)',norm(c))
    if not m: continue
    ws=set(w for w in m.group(1).split() if len(w)>4)
    if not ws: continue
    sc=[(len(ws&set(norm(rc).split())),i) for i,rc in enumerate(rate)]
    s,i=max(sc)
    if s>=1: best[j]=i
short={j:re.sub(r'\s+',' ',re.search(r'interest in ([a-z /-]+)',norm(ons[j])).group(1)).strip()
       for j in best}
def one_pair(a,b,MINN,plant=0.0,shuffle=False,seed=1):
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<MINN: return None
    y=(V[idx,a]<V[idx,b]).astype(float)
    ra=np.nan_to_num(R[idx,best[a]]); rb=np.nan_to_num(R[idx,best[b]])
    C=np.c_[BASE[idx],zs(np.c_[ra,rb,ra-rb])]
    Xp=np.delete(P[idx],[best[a],best[b]],axis=1)
    D0=np.c_[np.ones(len(idx)),C]
    Xr=Xp-D0@np.linalg.lstsq(D0,Xp,rcond=None)[0]
    pref=np.sign(ra-rb)
    if (pref>0).sum()<60 or (pref<0).sum()<60: return None
    w=Xr[pref>0].mean(0)-Xr[pref<0].mean(0); nw=np.linalg.norm(w)
    if nw<1e-9: return None
    w=w/nw
    proj=Xr@w; proj=(proj-proj.mean())/(proj.std()+1e-9)
    lab=y.copy()
    if shuffle:
        st=(male[idx]>0).astype(int)*3+np.digitize(breadth[idx],np.quantile(breadth[idx],[.33,.66]))
        rp=np.random.default_rng(seed+5)
        for s in np.unique(st):
            wq=np.flatnonzero(st==s)
            if len(wq)>1: lab[wq]=y[wq][rp.permutation(len(wq))]
    if plant>0: proj=proj+plant*lab                      # 真的把 A 先组沿 w 推 g 个标准差
    return dict(a=a,b=b,n=len(idx),
                v_shift=float(proj[lab==1].mean()-proj[lab==0].mean()),
                na=short[a],nb=short[b])
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)
       if KIND[a]==KIND[b] and a in best and b in best]
print(f"同类且能对上评分列的对: {len(pairs)}",flush=True)
rows=[]
for MINN,tag in [(400,'n>=400'),(250,'n>=250')]:
    for arm,pl,sh in [('real',0.,False),('shuf',0.,True),
                      ('plant0',0.,False),('plant05',0.05,False),('plant15',0.15,False)]:
        got=[]
        for a,b in pairs:
            r=one_pair(a,b,MINN,plant=pl,shuffle=sh)
            if r: got.append(dict(thresh=tag,arm=arm,**r))
        rows+=got
    print(f"  {tag}: {len([r for r in rows if r['thresh']==tag and r['arm']=='real'])} 对",flush=True)
D=check_columns(pd.DataFrame(rows),'R14')
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False)
se=lambda v: np.std(v)/np.sqrt(len(v))
S=D.groupby(['thresh','arm']).v_shift.agg(['size','mean',lambda s:np.std(s)/np.sqrt(len(s))])
S.columns=['k','shift','se']; S['ratio']=S['shift']/S['se']
print("\n=== 方向位移,两个门槛 ===")
print(S.reindex(pd.MultiIndex.from_product([['n>=400','n>=250'],
      ['real','shuf','plant0','plant05','plant15']])).round(4).to_string())
for tag in ['n>=400','n>=250']:
    T=S.loc[tag]
    g=Gate(f"#117c 加强 @ {tag}:先来的把其余偏好拉向自己吗?")
    g.negative_control("顺序标签分层内打乱", null=T.loc['shuf','shift'], effect=T.loc['real','shift'])
    g.positive_control("种植 0.15 个标准差被测出", planted=T.loc['plant15','shift'],
                       floor=T.loc['plant0','shift'], spread=T.loc['plant15','se'])
    g.no_sign_crossing("种植阶梯(相对 g=0)",
        [T.loc['plant05','shift']-T.loc['plant0','shift'],
         T.loc['plant15','shift']-T.loc['plant0','shift']])
    g.asserted("g=0 不开火(等于真实臂)", abs(T.loc['plant0','shift']-T.loc['real','shift'])<1e-9,
               f"{T.loc['plant0','shift']:+.4f} vs {T.loc['real','shift']:+.4f}")
    g.resolvable("方向位移", effect=T.loc['real','shift'], spread=T.loc['real','se'])
    print(); print(g)
    print(f"   位移为正的对: {(D[(D.thresh==tag)&(D.arm=='real')].v_shift>0).sum()}"
          f"/{int(T.loc['real','k'])}   强度 {T.loc['real','ratio']:.1f}x")
r400=S.loc[('n>=400','real'),'ratio']; r250=S.loc[('n>=250','real'),'ratio']
print(f"\n  预先承诺的判定:{r400:.1f}x -> {r250:.1f}x")
print(f"  -> {'加强成功,按新强度引用' if r250>=3 else '未达 3x,按 2.1x 引用并停止(不追第三轮 #111c)'}")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
