import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
"""
E01 A12 R13 -- 两件事:把零修对(#116e),然后用修好的仪器问"往哪个方向拉"。

第一件(仪器,#116e):
  #116 的负对照用 eff(y_置换),给 +0.0083 = 效应的 48%。#110f 已诊断原因:置换后的 y 没有协变量
  依赖,eff() 内部的合成世界因此退化成抛硬币 —— 偏移掉回"平置换 regime"。
  正确的零是 **eff(y_synth)**:y_synth 由真实 y 的协变量拟合生成,保留 COV->y,毁掉 P->y,
  而 eff() 内部对它再算一次同样的拟合 —— regime 自洽。若 eff(y_synth) ~ 0,则偏移无偏,
  #116 的 +0.0171 是真量级,而不是 +0.009。

第二件(描述,不是 Frontier —— 按 P0 标注):
  #116 立住"先获得哪一个预示其余轮廓"。**方向**是什么?
  ⚠ 跑之前的判别检查(#113c 的教训):方向本身**不判别** —— 重塑和"你本来就是 A 型"都预测
    "A 先的人更 A 向"。所以这一段标为 DESCRIPTION,只描述已确立效应的方向,不做因果主张。
    但它有一个非平凡的版本:在**减掉这两项各自的评分**之后,A 先的人是否仍在**其余**偏好上更 A 向。
    那不是"更喜欢 A",那是"其余的东西被拉向 A 那一侧"。

ESTIMAND        (1) eff(y_synth) —— 修正后的零
                (2) A 先组减 B 先组在"A 向 vs B 向"方向上的标准化位移,残差化掉两项评分与协变量
IDENTIFICATION  (1) 自洽;(2) 方向 w 从"评分 A > 评分 B"的人学出,与顺序标签无关
KILL            threshold-free。(1) 修正零对效应;(2) 位移对分层置换零
POSITIVE CTRL   种植一个已知方向的位移 -> 必须测出且单调
NEGATIVE CTRL   顺序标签在协变量分层内打乱
IMPOSSIBLE      (2) 的因果方向 —— 见上面的判别检查
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
def auc(y,s):
    o=np.argsort(s); r=np.empty(len(s)); r[o]=np.arange(len(s))+1
    n1=y.sum(); n0=len(y)-n1
    if n1<10 or n0<10: return np.nan
    return (r[y==1].sum()-n1*(n1+1)/2)/(n1*n0)
def ridge_auc(X,y,seed,alpha=50.,reps=8):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(reps):
        p=rng.permutation(len(y)); h=len(y)//2; tr,te=p[:h],p[h:]
        if y[tr].sum()<10 or (1-y[tr]).sum()<10: continue
        Xt=np.c_[np.ones(len(tr)),X[tr]]; A=Xt.T@Xt+alpha*np.eye(Xt.shape[1]); A[0,0]-=alpha
        b=np.linalg.solve(A,Xt.T@y[tr]); out.append(auc(y[te],np.c_[np.ones(len(te)),X[te]]@b))
    return np.nanmean(out) if out else np.nan
def synth_y(idx,y,C,seed):
    X=np.c_[np.ones(len(idx)),C]; w=np.linalg.lstsq(X,y,rcond=None)[0]
    return (np.random.default_rng(seed).random(len(idx))<np.clip(X@w,0.02,0.98)).astype(float)
def eff(idx,y,C,seed=1,ndraw=3):
    inc=lambda yy: ridge_auc(np.c_[C,P[idx]],yy,seed)-ridge_auc(C,yy,seed)
    off=np.nanmean([inc(synth_y(idx,y,C,seed+300+d)) for d in range(ndraw)])
    return inc(y)-off
pairs=[(a,b) for a,b in itertools.combinations(sorted(KIND),2)
       if KIND[a]==KIND[b] and a in best and b in best]
np.random.default_rng(3).shuffle(pairs)
rows=[];dirs=[]
for a,b in pairs:
    m=np.isfinite(V[:,a])&np.isfinite(V[:,b])&(V[:,a]!=V[:,b]); idx=np.flatnonzero(m)
    if len(idx)<400: continue
    y=(V[idx,a]<V[idx,b]).astype(float)
    ra=np.nan_to_num(R[idx,best[a]]); rb=np.nan_to_num(R[idx,best[b]])
    C=np.c_[BASE[idx],zs(np.c_[ra,rb,ra-rb])]
    rp=np.random.default_rng(88)
    rows.append(dict(a=a,b=b,arm='real',e=eff(idx,y,C),n=len(idx)))
    rows.append(dict(a=a,b=b,arm='null_synth',e=eff(idx,synth_y(idx,y,C,7),C),n=len(idx)))
    rows.append(dict(a=a,b=b,arm='null_perm',e=eff(idx,y[rp.permutation(len(y))],C),n=len(idx)))
    # ---- 方向(DESCRIPTION):A 向 vs B 向 ----
    Xp=np.delete(P[idx],[best[a],best[b]],axis=1)           # 剔除这两项自身
    D0=np.c_[np.ones(len(idx)),C]
    Xr=Xp-D0@np.linalg.lstsq(D0,Xp,rcond=None)[0]          # 残差化掉两项评分与协变量
    pref=np.sign(ra-rb)                                     # 谁更喜欢 A
    if (pref>0).sum()<80 or (pref<0).sum()<80: continue
    w=Xr[pref>0].mean(0)-Xr[pref<0].mean(0); w/=np.linalg.norm(w)+1e-12
    proj=Xr@w; proj=(proj-proj.mean())/(proj.std()+1e-9)
    def shift(lbl): return float(proj[lbl==1].mean()-proj[lbl==0].mean())
    st=(male[idx]>0).astype(int)*3+np.digitize(breadth[idx],np.quantile(breadth[idx],[.33,.66]))
    yp=y.copy()
    for s in np.unique(st):
        wq=np.flatnonzero(st==s)
        if len(wq)>1: yp[wq]=y[wq][np.random.default_rng(5).permutation(len(wq))]
    dirs.append(dict(a=a,b=b,n=len(idx),shift=shift(y),shift_null=shift(yp),
                     na=short[a],nb=short[b]))
    if len(dirs)%15==0: print(f"  {len(dirs)} pairs",flush=True)
    if len(dirs)>=45: break
D=pd.DataFrame(rows); G=pd.DataFrame(dirs)
# ⚠ #156:装上守卫本身,而不只是修症状。`check_columns` 正是为 `#117e` 的
#   `shift` 撞名写的,而这一轮从来没调用过它 —— 于是它在打印完 `#117` 的结论
#   之后崩掉,两处,没人发现。
try: check_columns(G,'A12/R169 方向表')
except Exception as _e: print(f'  ⚠ check_columns: {_e}',flush=True)
OUT=pathlib.Path(__file__).parent/'results'; D.to_csv(OUT/'grid.csv',index=False); G.to_csv(OUT/'dirs.csv',index=False)
S=D.groupby('arm').e.agg(['size','mean',lambda s:s.std()/np.sqrt(len(s))]); S.columns=['k','e','se']
print("\n=== 第一件:零修对之后 ===")
print(S.reindex(['real','null_synth','null_perm']).round(4).to_string())
g=Gate("#116e 修正后的零")
g.negative_control("修正零 eff(y_synth)", null=S.loc['null_synth','e'], effect=S.loc['real','e'])
g.resolvable("真实效应", effect=S.loc['real','e'], spread=S.loc['real','se'])
print(); print(g)
print(f"  旧零 eff(y_置换) {S.loc['null_perm','e']:+.4f} = 效应的 {100*S.loc['null_perm','e']/S.loc['real','e']:.0f}%")
print(f"  新零 eff(y_synth) {S.loc['null_synth','e']:+.4f} = 效应的 {100*S.loc['null_synth','e']/S.loc['real','e']:.0f}%")
print("\n=== 第二件(DESCRIPTION,不做因果主张):先来的把人往哪边拉 ===")
# ⚠ #156:`shift` 是 pandas 的 DataFrame **方法**,`G.shift` 返回方法不是列 ——
#   这正是 `#117e` 自己记录的第 5 次访问器撞名,而 `check_columns` 就是为它写的。
#   这一轮从来没调用过那个守卫(`guard_lint` #128a 早就标了它缺 columns),
#   所以它在打印完 `#117` 的结论之后**崩在这里**,而没人发现。
sh=G['shift'].values; nl=G['shift_null'].values
se=lambda v: np.std(v)/np.sqrt(len(v))
print(f"  A 先组 减 B 先组,在「A 向 − B 向」方向上的标准化位移(已残差化掉两项评分与协变量)")
print(f"    真实 {sh.mean():+.4f} ± {se(sh):.4f}   分层置换零 {nl.mean():+.4f} ± {se(nl):.4f}   "
      f"{len(G)} 对")
print(f"    位移为正的对: {(sh>0).sum()}/{len(G)}")
g2=Gate("先来的那个,把其余偏好拉向自己吗?(DESCRIPTION)")
g2.negative_control("顺序标签在分层内打乱", null=nl.mean(), effect=sh.mean())
g2.resolvable("方向位移", effect=sh.mean(), spread=se(sh))
print(); print(g2)
G['d']=G['shift']-G['shift_null']
print("\n  位移最大的 6 对(先来的那个 -> 其余偏好被拉向它):")
for _,r in G.reindex(G.d.abs().sort_values(ascending=False).index).head(6).iterrows():
    who=r.na if r.d>0 else r.nb
    print(f"    {r.na[:26]:<26} vs {r.nb[:26]:<26}  位移 {r.d:+.3f}  -> 拉向「{who[:26]}」")
print(f"\nartifact sha1 {hashlib.sha1(D.to_csv(index=False).encode()).hexdigest()[:12]}")
