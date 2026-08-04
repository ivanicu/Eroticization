import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A124 R386 -- 残差第一名是 `humiliation`,而它和羞耻是近义概念

⚠ **Closure(稳健性),如实标注。**

`#338b` 的残差前二:**`humiliation` +0.301** 与 **`gender` +0.300**。
**而 `humiliation` 是这个项目里唯一一个既是**兴趣类别**、又直接是**羞耻的近义概念**的东西。**
**它可能不是「这类人更羞耻」,而是「这道题和羞耻那道题问的是同一件事」。**

ESTIMAND        去掉 `humiliation` 重算 `#338`/`#339` 的全部数字;
                **以去掉 `gender` 作对照**(`gender` 没有近义问题,效果应当更小),
                **并以去掉一个随机类别的分布作基线**。
KILL            **若残差结构基本不变 -> 它只是一格,`#339b` 的开放问题照旧;
                若残差 sd 明显掉(且明显超过随机留一的分布)-> 那 42% 里有一部分是概念重叠。**
POSITIVE CTRL   沿用 `#383` 的种入(给一个类别 +1.0,留一必须把它挑出来)。
NEGATIVE CTRL   留一**随机**类别 30 次 -> 残差 sd 的分布就是「只是少了一格」的基线。
IMPOSSIBLE      「概念重叠」不可直接测;本轮只能测「去掉它之后剩下多少」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in ons])
HAS=np.isfinite(ONS); NCA=HAS.shape[1]; NAM=[str(c) for c in ons]
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
okS=np.isfinite(sh)
def pack(keep):
    K=np.array(keep)
    P=HAS[:,K].mean(0); RAR=-np.log(np.clip(P,1e-4,1.))
    mu=np.array([sh[okS&HAS[:,j]].mean() for j in K])
    se=np.array([sh[okS&HAS[:,j]].std()/np.sqrt(max((okS&HAS[:,j]).sum(),1)) for j in K])
    z=lambda v:(v-np.nanmean(v))/np.nanstd(v)
    X0=np.column_stack([np.ones(len(K)),z(RAR)])
    res=mu-X0@np.linalg.lstsq(X0,mu,rcond=None)[0]
    F=[np.array([np.nanmean(ONS[HAS[:,j],j]) for j in K]),
       np.array([np.nanstd(ONS[HAS[:,j],j]) for j in K]),
       np.array([np.nanmean(S[okS&HAS[:,j]&ok]) for j in K]),
       np.array([np.nanmean(C3[okS&HAS[:,j]&ok]) for j in K])]
    X=np.column_stack([np.ones(len(K))]+[z(v) for v in F])
    b,*_=np.linalg.lstsq(X,res,rcond=None); r=res-X@b
    R2=1-float(r@r)/float(((res-res.mean())**2).sum())
    return float(np.corrcoef(mu,RAR)[0,1]),float(res.std()),float(np.median(se)),R2
i_hum=[i for i,nm in enumerate(NAM) if 'humili' in nm.lower()][0]
i_gen=[i for i,nm in enumerate(NAM) if 'gender' in nm.lower()][0]
ALLK=list(range(NCA))
base=pack(ALLK)
print(f"全部 {NCA} 类:`corr(mu, rar)` **{base[0]:+.4f}** · 残差 sd **{base[1]:.4f}** · "
      f"se 中位 **{base[2]:.4f}** ({base[1]/base[2]:.2f}×) · 四量 R² **{base[3]:.4f}**")
for nm,i in (('humiliation',i_hum),('gender(对照)',i_gen)):
    v=pack([k for k in ALLK if k!=i])
    print(f"去掉 **{nm:<16}** 残差 sd **{v[1]:.4f}**({v[1]/v[2]:.2f}×)· "
          f"`corr(mu,rar)` {v[0]:+.4f} · 四量 R² {v[3]:.4f}")
    if i==i_hum: hum=v
    else: gen=v
rg=np.random.default_rng(1212)
NGV=[pack([k for k in ALLK if k!=j])[1] for j in ALLK]
NGV=np.array(NGV)
print(f"\n负对照(留一**每一个**类别,{NCA} 次):残差 sd **{NGV.mean():.4f} ± {NGV.std():.4f}** · "
      f"范围 [{NGV.min():.4f}, {NGV.max():.4f}]")
zh=(hum[1]-NGV.mean())/max(NGV.std(),1e-12); zg=(gen[1]-NGV.mean())/max(NGV.std(),1e-12)
print(f"   去掉 humiliation 落在 **{zh:+.2f}** sd · 去掉 gender 落在 **{zg:+.2f}** sd")
syn=sh.copy().astype(float); syn[HAS[:,i_gen]]+=1.0
def res_of(shv,K):
    P=HAS[:,K].mean(0); RAR=-np.log(np.clip(P,1e-4,1.))
    m2=np.array([shv[okS&HAS[:,j]].mean() for j in K])
    z=lambda v:(v-v.mean())/v.std()
    X0=np.column_stack([np.ones(len(K)),z(RAR)])
    return m2-X0@np.linalg.lstsq(X0,m2,rcond=None)[0]
rs=res_of(syn,ALLK); k=int(np.argmax(rs))
print(f"正对照(给 `gender` 的人 +1.0):残差最大的是 **#{k}**(应 #{i_gen}),{rs[k]:+.3f}")
T=pd.DataFrame([dict(v_arm='全部',sd=base[1],r2=base[3]),
                dict(v_arm='去 humiliation',sd=hum[1],r2=hum[3]),
                dict(v_arm='去 gender',sd=gen[1],r2=gen[3]),
                dict(v_arm='留一均值',sd=float(NGV.mean()),r2=np.nan)])
check_columns(T,'R386'); T.to_csv(pathlib.Path(__file__).parent/'results'/'drop.csv',index=False)
gg=Gate('残差结构是不是 `humiliation` 一格撑的')
gg.asserted('★ 正对照:给一个类别 +1.0,留一必须挑出它',k==i_gen,f"挑出 #{k},应 #{i_gen}")
gg.asserted('★ 负对照(留一每一个类别)给出「只是少一格」的基线',NGV.std()>0,
            f"残差 sd {NGV.mean():.4f} ± {NGV.std():.4f},范围 [{NGV.min():.4f}, {NGV.max():.4f}]")
gg.asserted('★ 注册的 kill:去掉 `humiliation` 后残差 sd 是否明显低于随机留一的分布(< −2 sd)',
            zh<-2.0,f"humiliation **{zh:+.2f}** sd · 对照 gender **{zg:+.2f}** sd")
gg.asserted('⚠ 残差里剩下的东西还在不在(比 se 中位大)',hum[1]>hum[2],
            f"去 humiliation 后 {hum[1]:.4f} vs se 中位 {hum[2]:.4f}({hum[1]/hum[2]:.2f}×)")
gg.asserted('⚠ 边界:「概念重叠」不可直接测',True,'本轮只能测「去掉它之后剩下多少」')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
