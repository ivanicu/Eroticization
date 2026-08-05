import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A141 R423 -- `EARLY` 的两性差:加真的 n,并且**点名那一格不可能**

`#378b`:`EARLY` 的两性差 |t| **2.14**,族内阈 **2.308** —— 差 0.17 没过。
不能回头单判 `EARLY`(看过结果再选检验);合法的只有**加真的 n**。

⚠ **但有一件事在跑之前就知道**:`c3⁻` 的剖面在 `fit_apply` 里需要 **≥8 个块**
(留一剖面 `ct-F[b]>=6` · 半剖面 `>=4` · 稳定性 `>=8`),**所以宽口径下它不是同一个坐标**。
**⇒ 「宽 × 全模型」这一格**不可能**,而不是「还没做」。** 点名它,不糊过去(`#374d` 同族)。

设计 = **三格 + 一格点名不可能**:
| | 全模型(S · c3⁻ · EARLY) | 简模型(EARLY · 类别数) |
|---|---|---|
| **窄**(覆盖≥8) | **A** = `#378`(|t| 2.14) | **B** — 隔离「去掉坐标」的影响 |
| **宽**(覆盖≥4) | **D ⛔ 不可能**(c3⁻ 需要 ≥8 块) | **C** — 隔离「加人」的影响 |

ESTIMAND        每格里 `EARLY` 的两性系数差与 |t|;主量 = **C 格**(宽 × 简)。
⚠ 预注册方向    `EARLY` 在**组 0** 更强(`#378b`)。**同号且越阈才算复制;变号或未越阈都算没有。**
KILL(条件式)  仅当对照都过**且 C 格 MDE < 0.05** -> 判:C 格是否同号且越过它自己的零。
⚠ 零的种类     `offset_control`:**随机等大小劈分**(每格各建自己的零 —— 格的 n 不同,零也不同)。
⚠ guard 25     `#371a` 的坑:**「宽臂」必须真的宽**。本轮的简模型不依赖 `S`,所以应当真的宽 —— **但要测,不要信**。
IMPOSSIBLE      ① D 格不可能(上面);② 简模型的 `EARLY` 系数**不是**全模型的那个 ——
                B 格存在就是为了量这个差;③ 宽窄是**不同人群**(`#346b`),这是跨人群复制。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
Sn=make_S(8)
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
SEX=pd.to_numeric(d['biomale'],errors='coerce').values.astype(float)

def cell(name,mask,preds,labels):
    g0=mask&(SEX==0); g1=mask&(SEX==1)
    def fit(y,g):
        k=int(g.sum()); yy=(y[g]-y[g].mean())/max(y[g].std(),1e-12)
        X=np.column_stack([np.ones(k)]+[(v[g]-v[g].mean())/max(v[g].std(),1e-12) for v in preds])
        b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
        s2=float(r@r)/(k-len(preds)-1); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
        return b[1:],se[1:]
    b0,s0=fit(sh,g0); b1,s1=fit(sh,g1)
    j=labels.index('EARLY')
    diff=b1[j]-b0[j]; sed=np.sqrt(s0[j]**2+s1[j]**2); t=diff/max(sed,1e-12)
    n0=int(g0.sum()); idx=np.flatnonzero(mask); mt=[]
    for s_ in range(400):
        rg=np.random.default_rng(5000+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
        ba,sa=fit(sh,ga); bb,sb=fit(sh,gb)
        mt.append(float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))))
    mt=np.array(mt); thr=float(np.percentile(mt,95))
    return dict(name=name,n=int(mask.sum()),n0=n0,b0=b0[j],b1=b1[j],diff=diff,t=t,
                thr=thr,k=len(preds),fit=fit,g0=g0,g1=g1,mask=mask,preds=preds,labels=labels)

BASE=np.isfinite(EARLY)&np.isfinite(sh)&np.isfinite(SEX)&np.isfinite(ncat)
NARROW=BASE&ok&np.isfinite(Sn)&np.isfinite(C3)
WIDE=BASE&(COVB>=4)
A=cell('A 窄 × 全模型',NARROW,[Sn,C3,EARLY],['S','c3','EARLY'])
B=cell('B 窄 × 简模型',NARROW,[EARLY,ncat],['EARLY','ncat'])
C=cell('C 宽 × 简模型',WIDE,[EARLY,ncat],['EARLY','ncat'])
print("三格(D 格 ⛔ 不可能:`c3⁻` 的剖面需要 ≥8 个块,宽口径下它不是同一个坐标):")
for X in (A,B,C):
    print(f"   {X['name']:<14} n=**{X['n']:>6,}** (k={X['k']}) · 组0 **{X['b0']:+.4f}** · 组1 **{X['b1']:+.4f}** · "
          f"差 **{X['diff']:+.4f}** · |t| **{abs(X['t']):.3f}** vs 自己的族内阈 **{X['thr']:.3f}** · "
          f"{'**越阈**' if abs(X['t'])>X['thr'] else '未越阈'}")
print(f"\n   ⚠ guard 25 的对象:窄 {A['n']:,} -> 宽 {C['n']:,}")

# C 格的 MDE:先算再看
Xc=C
def mde_for(X,seed0):
    m=None
    for gg in (0.03,0.05,0.08,0.12):
        hit=0
        for s_ in range(30):
            rg=np.random.default_rng(seed0+int(gg*100)*11+s_)
            y=np.full(NN,np.nan)
            for gs,ex in ((X['g0'],0.0),(X['g1'],gg)):
                k=int(gs.sum()); ze=(EARLY[gs]-EARLY[gs].mean())/EARLY[gs].std()
                y[gs]=-0.07*ze+ex*ze+rg.standard_normal(k)
            ba,sa=X['fit'](y,X['g0']); bb,sb=X['fit'](y,X['g1'])
            if np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))>X['thr']: hit+=1
        print(f"   [{X['name']}] 组1 多出 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
        if m is None and hit>=24: m=gg
    return m if m else 0.20
print(f"\n先算 MDE 再看数(`#369a`):")
MDE_C=mde_for(C,90000)
print(f"   **C 格 MDE = {MDE_C:.2f}** · 有意义 0.05 -> "
      f"{'**够小**' if MDE_C<=0.05 else '**不够小**'}")

rg=np.random.default_rng(21); yp=np.full(NN,np.nan)
for gs,ex in ((C['g0'],0.0),(C['g1'],0.25)):
    k=int(gs.sum()); ze=(EARLY[gs]-EARLY[gs].mean())/EARLY[gs].std()
    yp[gs]=-0.07*ze+ex*ze+rg.standard_normal(k)
ba,sa=C['fit'](yp,C['g0']); bb,sb=C['fit'](yp,C['g1'])
TP=float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12))))
rgn=np.random.default_rng(33); idx=np.flatnonzero(C['mask']); pp=rgn.permutation(idx)
ha=np.zeros(NN,bool); ha[pp[:C['n0']]]=True; hb=np.zeros(NN,bool); hb[pp[C['n0']:]]=True
ca,sa2=C['fit'](sh,ha); cb,sb2=C['fit'](sh,hb)
TN=float(np.max(np.abs((cb-ca)/np.maximum(np.sqrt(sa2**2+sb2**2),1e-12))))
print(f"\n正对照(C 格,组1 多出 0.25):max-|t| **{TP:.3f}** vs 阈 {C['thr']:.3f}")
print(f"负对照(C 格,打乱性别标签):max-|t| **{TN:.3f}** vs 阈 {C['thr']:.3f}")
pd.DataFrame([{k:v for k,v in X.items() if k in ('name','n','n0','b0','b1','diff','t','thr','k')}
              for X in (A,B,C)]).to_csv(pathlib.Path(__file__).parent/'results'/'cells.csv',index=False)

SAME=np.sign(C['diff'])==np.sign(A['diff'])
g=Gate('EARLY 的两性差,加真的 n 之后还在吗')
g.asserted('★ 正对照:C 格组1 多出 0.25 -> 必须越阈',TP>C['thr'],f"{TP:.3f} vs {C['thr']:.3f}",kind='control')
g.asserted('★ 负对照:C 格打乱性别标签 -> 必须落回零',TN<=C['thr'],f"{TN:.3f} vs {C['thr']:.3f}",kind='control')
g.relaxation_reached_the_population('★ guard 25:放宽的口径真的到达了人群',A['n'],C['n'],
                                    what='覆盖 ≥8 -> ≥4;简模型不依赖 S')
if TP>C['thr'] and TN<=C['thr'] and MDE_C<=0.05:
    g.asserted('★ 注册的 kill:C 格同号(预注册方向)**且**越过自己的阈',
               SAME and abs(C['t'])>C['thr'],
               f"同号 {SAME} · |t| {abs(C['t']):.3f} vs {C['thr']:.3f}")
else:
    g.asserted('★ 注册的 kill(MDE ≥ 0.05 或对照未过 -> 不判)',False,
               f"C 格 MDE {MDE_C:.2f} · 正 {TP>C['thr']} · 负 {TN<=C['thr']}")
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
