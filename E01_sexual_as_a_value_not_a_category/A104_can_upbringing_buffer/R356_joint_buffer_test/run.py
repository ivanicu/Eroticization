import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A104 R356 -- 四个「本该缓冲」的变量,合起来有没有在缓冲

`#310c`:`#309b` 有功效而 `#310a` 没有,差别在**聚合**。**那就把缓冲问题也做成聚合量。**

⚠ **先写下,免得结果一出来就去挑最大的那个**(`#309c` 神经质那格的同一个陷阱):
**联合检验只能说「有没有任何一个在缓冲」,说不出是哪个。**
若联合检验开火,下一轮才是定位,而定位需要它自己的多重性处理。

ESTIMAND        四个「本该缓冲」的变量 —— **成长期性开放度 · 关系风格 · 年龄 · 开放性** ——
                与两条路的**全部 8 个交互项**,联合为零吗?
                统计量 = 去掉这 8 项后的 **R² 下降**;零 = `perm_finite` 打乱人。
KILL            **若 R² 下降明显超过置换零 -> 至少有一个在缓冲,下一轮定位;
                若没有 -> 报**联合 MDE**:只在一个变量上植入多大的缓冲才抓得到。**
POSITIVE CTRL   只在**一个**变量上植入缓冲,强度扫描 -> 联合检验必须抓到
                (这同时证明**聚合确实换来了功效**,而不只是换了个数)。
NEGATIVE CTRL   `perm_finite` 题内跨人打乱。
⚠ 报           **置换零的分布**,不报解析 p 值 —— 设计里 8 个交互项彼此共线。
IMPOSSIBLE      联合检验对「所有变量都缓冲一点点」最敏感,对「只有一个缓冲很多」最不敏感;
                它与 `#310a` 的单系数检验**互补**,不是它的升级。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A97_is_the_headline_optimistic/R347_nested_cv/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('rgF=np.random.default_rng')[0])

SHAME=next(c for c in d.columns if 'ashamed' in str(c))
sh=pd.to_numeric(d[SHAME],errors='coerce').values.astype(float)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
MODS={'成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map(
        {'Repressed':-1.,'Neutral':0.,'Liberated':1.}).values.astype(float),
      '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map(
        {'Monogamous':0.,'Not monogamous':1.}).values.astype(float),
      '年龄':d['age'].map(AGE).values.astype(float),
      '开放性':pd.to_numeric(d['opennessvariable'],errors='coerce').values.astype(float)}
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&ok
for v in MODS.values(): m0&=np.isfinite(v)
def resid(a,b,m):
    out=np.full(NN,np.nan); x=b[m]; x=(x-x.mean())/x.std()
    out[m]=a[m]-np.polyval(np.polyfit(x,a[m],1),x); return out
RS=resid(S,C3,m0); RC=resid(C3,S,m0); n=int(m0.sum())
z=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zr,zc=z(RS),z(RC); ZM={k:z(v) for k,v in MODS.items()}
BASE=np.column_stack([np.ones(n),zr,zc]+[ZM[k] for k in MODS])
INT=np.column_stack([zr*ZM[k] for k in MODS]+[zc*ZM[k] for k in MODS])
print(f"n={n:,} · 4 个变量 × 2 条路 = **{INT.shape[1]} 个交互项**")
print("⚠ 先写下:**联合检验只能说「有没有任何一个在缓冲」,说不出是哪个。**")
def r2(X,y): 
    b,*_=np.linalg.lstsq(X,y,rcond=None); r=y-X@b; return 1-float(r@r)/float(((y-y.mean())**2).sum())
def stat(y):
    yy=(y-y.mean())/y.std()
    return r2(np.column_stack([BASE,INT]),yy)-r2(BASE,yy)
def perm_finite_local(v,seed):
    zz=v.copy(); j=np.flatnonzero(np.isfinite(zz))
    zz[j]=zz[np.random.default_rng(seed).permutation(j)]; return zz
obs=stat(sh[m0])
nul=np.array([stat(perm_finite_local(sh,500+i)[m0]) for i in range(60)])
print(f"\n★ 观测 R² 下降 **{100*obs:.4f}pp** · 置换零 **{100*nul.mean():.4f} ± {100*nul.std():.4f}pp** "
      f"-> **{(obs-nul.mean())/max(2*nul.std(),1e-12):.2f}× 的 2×展布** · "
      f"零里 ≥ 观测的比例 **{100*float((nul>=obs).mean()):.1f}%**")
rg=np.random.default_rng(77)
print(f"\n正对照:只在**一个**变量(成长期性开放度)上植入缓冲,强度扫描")
U=ZM['成长期性开放度']; SW=[]
for buf in (0.0,0.10,0.20,0.30,0.50):
    y=0.12*zr*(1-buf*U)+0.12*zc*(1-buf*U)+rg.standard_normal(n)
    s_=stat(y); SW.append((buf,s_))
    print(f"   缓冲 {100*buf:>3.0f}%: R² 下降 **{100*s_:.4f}pp** "
          f"({(s_-nul.mean())/max(2*nul.std(),1e-12):>5.2f}× 的 2×展布)")
det=[b for b,s_ in SW if (s_-nul.mean())>2*nul.std()]
MDE=min(det) if det else None
print(f"   -> **联合 MDE ≈ {100*MDE:.0f}%**" if MDE is not None
      else "   -> **扫描到 50% 都没抓到 -> 聚合也没换来功效**")
print(f"   对比 `#310a` 的单系数 MDE:**54.8%**")
# ⚠ 上面的 MDE 只在**成长期性开放度**上测过。四个变量的方差与分布不同 -> 逐个测。
print(f"\n⚠ 逐变量的联合 MDE(每次只在那一个变量上植入):")
MDEV={}
for kk_ in MODS:
    Uk=ZM[kk_]; got=None
    for buf in (0.10,0.15,0.20,0.30,0.50):
        y=0.12*zr*(1-buf*Uk)+0.12*zc*(1-buf*Uk)+np.random.default_rng(91).standard_normal(n)
        if (stat(y)-nul.mean())>2*nul.std(): got=buf; break
    MDEV[kk_]=got
    print(f"   {kk_:<16} 联合 MDE ≈ **{'%.0f%%'%(100*got) if got else '> 50%'}**")
WORST=max([v for v in MDEV.values() if v] or [None]) if any(MDEV.values()) else None
print(f"   -> **最保守的一个:{'%.0f%%'%(100*WORST) if WORST else '> 50%'}**"
      f"{'(有变量到 50% 都抓不到)' if any(v is None for v in MDEV.values()) else ''}")
# ⚠ #300a:换一个**聚合方式**当旋钮 —— R² 下降 vs 8 个交互项里的 max|t|。
def stat2(y):
    yy=(y-y.mean())/y.std(); X=np.column_stack([BASE,INT])
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(len(yy)-X.shape[1]); se_=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    t=np.abs(b[BASE.shape[1]:]/np.maximum(se_[BASE.shape[1]:],1e-12)); return float(t.max())
o2=stat2(sh[m0]); n2=np.array([stat2(perm_finite_local(sh,500+i)[m0]) for i in range(60)])
print(f"\n发明的旋钮 · 换聚合方式:max|t| 观测 **{o2:.3f}** vs 置换零 **{n2.mean():.3f} ± {n2.std():.3f}** "
      f"-> **{(o2-n2.mean())/max(2*n2.std(),1e-12):.2f}× 的 2×展布** · "
      f"零里 ≥ 观测的 **{100*float((n2>=o2).mean()):.1f}%**")

T=pd.DataFrame([dict(v_kind='观测',v_val=obs)]+[dict(v_kind=f'植入{int(100*b)}%',v_val=s_) for b,s_ in SW]
               +[dict(v_kind='零均值',v_val=float(nul.mean())),dict(v_kind='零sd',v_val=float(nul.std()))])
check_columns(T,'R356'); T.to_csv(pathlib.Path(__file__).parent/'results'/'joint.csv',index=False)
gg=Gate('四个本该缓冲的变量,合起来有没有在缓冲')
gg.plant_direction_from_sweep('★ 正对照:单变量缓冲强度扫描 -> 联合统计量必须随强度上升',
                              SW,SW[0][1],baseline_spread=2*float(nul.std()))
gg.negative_control('★ 负对照:`perm_finite` 打乱人',float(nul.mean()),obs,
    null_spread=float(nul.std()),null_kind='`perm_finite` 题内跨人打乱 —— 保住缺失格局(#264b)')
gg.asserted('★ 发明的旋钮:换成 max|t| 的聚合方式,结论一不一样',
            ((o2-n2.mean())<2*n2.std())==((obs-nul.mean())<2*nul.std()),
            f"R² 下降 {(obs-nul.mean())/max(2*nul.std(),1e-12):+.2f}× · "
            f"max|t| {(o2-n2.mean())/max(2*n2.std(),1e-12):+.2f}× —— 两种聚合同一结论")
gg.asserted('★ 注册的 kill:联合 R² 下降是否明显超过置换零',
            (obs-nul.mean())>2*nul.std(),
            f"观测 {100*obs:.4f}pp vs 零 {100*nul.mean():.4f} ± {100*nul.std():.4f}pp")
gg.asserted('★ 逐变量 MDE:四个都测过了吗(不能只用一个变量的 MDE 代表全体)',
            all(v is not None for v in MDEV.values()),
            ' · '.join(f"{k} {'%.0f%%'%(100*v) if v else '>50%'}" for k,v in MDEV.items()))
gg.asserted('★ 聚合有没有换来功效(联合 MDE vs `#310a` 的单系数 MDE 54.8%)',
            MDE is not None and 100*MDE<54.8,
            f"联合 MDE ≈ {100*MDE:.0f}%" if MDE is not None else "扫描到 50% 都没抓到")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
