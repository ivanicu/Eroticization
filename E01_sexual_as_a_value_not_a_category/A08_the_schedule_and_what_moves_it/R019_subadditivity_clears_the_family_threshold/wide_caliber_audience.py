import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A159 R455 -- 那个「近两倍」是不是噪声:宽口径 + 简模型,带窄口径参照臂

`#410b`:差 −0.0361、|t| 1.294,MDE 0.08 > 有意义 0.05 —— **既不能说有,也不能说没有**。
`#410c`:**加 n 修不了混淆**(没有性伴的人既没观众也没机会)——
**但它能修「是不是噪声」,而那正是现在挡路的那一个。**

三格 + 一格**点名不可能**(`#379b` 的做法):
| | 全模型(含 `S`·`c3⁻`) | 简模型(年龄 + 类别数) |
|---|---|---|
| **窄**(覆盖≥8) | **A** = `#410a`(|t| 1.294) | **B** — 隔离「去掉坐标」 |
| **宽**(覆盖≥4) | **D ⛔ 不可能**(`S`/`c3⁻` 需 ≥8 块) | **C** — 隔离「加人」 |

判据(**先标支**,`#379c`)
                【两支】**guard 25**(宽臂必须真的宽,`#371a` 的坑)· 负对照用**越阈率** ·
                        guard 26 **显式传 branch**。
                【非零支】**预注册方向:差为负**(有过性伴的人更强)——
                          **C 格同号且越阈才算;变号或未越阈都算没有。**
                【零支】未越阈时启用 MDE。
⚠ 零的种类     `offset_control`:**随机等大小分层**(每格各建自己的零)。
⚠ 结论的携带   **无论结果如何,`#410c` 的混淆原样带进结论** —— 本轮只解决「是不是噪声」。
IMPOSSIBLE      ① D 格不可能;② 简模型的系数**不是**全模型的那个(B 格就是为量这个差存在的);
                ③ 宽窄是**不同人群**(`#346b`)-> **跨人群**,不是加功率。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
AC=next(c for c in d.columns if '41kpfir' in str(c)); ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
SXMAP={'0':0.,'1-2':1.5,'3-7':5.,'8-20':14.,'21+':25.}
SX=d['sexcount'].map(SXMAP).values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
BASE=np.isfinite(sh)&np.isfinite(ACTED)&np.isfinite(SX)&np.isfinite(AG)&np.isfinite(ncat)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def cell(name,mask,preds):
    g0=mask&(SX==0); g1=mask&(SX>0)
    def fit(y,g):
        k=int(g.sum())
        X=np.column_stack([np.ones(k)]+[z(v,g) for v in [sh]+preds])
        yy=z(y,g); b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
        s2=float(r@r)/(k-len(preds)-2); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
        return float(b[1]),float(se[1])
    b0,s0=fit(ACTED,g0); b1,s1=fit(ACTED,g1); diff=b1-b0
    idx=np.flatnonzero(mask); n0=int(g0.sum()); nul=[]
    for s_ in range(400):
        rg=np.random.default_rng(6900+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
        nul.append(fit(ACTED,gb)[0]-fit(ACTED,ga)[0])
    nul=np.array(nul); thr=float(np.percentile(np.abs(nul),95))
    return dict(name=name,n=int(mask.sum()),n0=n0,n1=int(g1.sum()),b0=b0,b1=b1,
                diff=diff,thr=thr,fit=fit,g0=g0,g1=g1,mask=mask,nulsd=float(nul.std()))
NARROW=BASE&ok&np.isfinite(S)&np.isfinite(C3)
WIDE=BASE&(COVB>=4)
A=cell('A 窄 × 全模型',NARROW,[S,C3,ncat,AG])
B=cell('B 窄 × 简模型',NARROW,[ncat,AG])
C=cell('C 宽 × 简模型',WIDE,[ncat,AG])
print("三格(D 格 ⛔ 不可能:`S`/`c3⁻` 需 ≥8 块):")
for X in (A,B,C):
    print(f"   {X['name']:<14} n=**{X['n']:>6,}**(没有过 {X['n0']:,} · 有过 {X['n1']:,})· "
          f"没有过 **{X['b0']:+.4f}** · 有过 **{X['b1']:+.4f}** · 差 **{X['diff']:+.4f}** vs 阈 **{X['thr']:.4f}** · "
          f"{'**越阈**' if abs(X['diff'])>X['thr'] else '未越阈'}")
negs=[]
idx=np.flatnonzero(C['mask']); n0=C['n0']
for s_ in range(200):
    rg=np.random.default_rng(99500+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    negs.append(C['fit'](ACTED,gb)[0]-C['fit'](ACTED,ga)[0])
negs=np.array(negs); rate=float((np.abs(negs)>C['thr']).mean())
print(f"\n负对照(**越阈率**,C 格随机分层 200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(C 格),每级 30 次:")
MDE=None
# ⚠ 第一版网格 (0.03, 0.05, 0.08, 0.12):0.05 给 66.7%、0.08 给 96.7%,**80% 点落在两者之间**,
# MDE 被报成 0.08 > 实测 0.0618,于是 guard 26 FAIL。
# **`#403b` 的规矩**:网格分辨率会把 MDE 系统性报**高**,而报高的方向恰好是让门失败的方向;
# **加密网格不是搬门柱** —— 改的是**测量的分辨率**,判据(实测 0.0618)一个字没动。
for gg in (0.03,0.05,0.055,0.06,0.065,0.07,0.08):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(40+int(gg*1000)*113+s_)
        y=np.full(NN,np.nan)
        for gs,ex in ((C['g0'],0.0),(C['g1'],gg)):
            k=int(gs.sum()); zs=z(sh,gs); y[gs]=-0.08*zs-ex*zs+rg.standard_normal(k)
        if abs(C['fit'](y,C['g1'])[0]-C['fit'](y,C['g0'])[0])>C['thr']: hit+=1
    print(f"   「有过」层多出 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.15
NONNULL=abs(C['diff'])>C['thr']
CONT=abs(C['diff']) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 **{CONT:.4f}**")
pd.DataFrame([{k:v for k,v in X.items() if k in ('name','n','n0','n1','b0','b1','diff','thr')}
              for X in (A,B,C)]).to_csv(pathlib.Path(__file__).parent/'results'/'cells.csv',index=False)
SAME=np.sign(C['diff'])<0
g=Gate('那个「近两倍」是不是噪声')
g.relaxation_reached_the_population('★【两支】guard 25:放宽的口径真的到达了人群',
                                    A['n'],C['n'],what='覆盖 ≥8 -> ≥4;简模型不依赖 S/c3⁻')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='MDE 扫描 80% 检出',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】C 格**同号(负,预注册)且越阈** -> 不是噪声',SAME and NONNULL,
               f"C 差 {C['diff']:+.4f} vs 阈 {C['thr']:.4f} · 同号 {SAME}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **`#410c` 原样带进结论**:没有性伴的人**既没观众也没机会** ——"
      f"**本轮只解决「是不是噪声」,不解决「是不是观众」。**")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
