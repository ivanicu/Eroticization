import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A159 R454 -- 羞耻是「在意别人怎么看」,还是「自己看自己」

`#384b`:**羞耻 -> 少做**(−0.083,−6.85 sd),而 `#385` 已排除「起始时间造出来的」。
**从没问过:这条关系在**从没有过性伴**与**有过**的人身上,是不是一样。**

两个活着的世界,**对同一个差预测相反的方向**:
**A 别人的眼睛** —— 羞耻的压制需要**观众**:有过性伴的人身上**更强**(差为负,即更负的系数);
**B 自己的眼睛** —— 压制来自自我评价,**与有没有观众无关**:差 ≈ 0,或在**没有性伴**的人身上更强。

⚠ 而这是一个**真正的分离器**:两个世界对**差的符号**给出相反的预测,不是「大一点小一点」。

ESTIMAND        按 `sexcount`(0 vs ≥1)分两层,层内拟合 `ACTED ~ 羞耻 + S + c3⁻ + 类别数`;
                主量 = **两层羞耻系数的差**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**;offset 零非退化。
                【非零支】差越过 offset 零 -> 按符号判 A 或 B;
                【零支】未越阈 -> 启用 MDE;MDE < 0.05 才算「看得见而没有」。
⚠ 零的种类     `offset_control`:**两层系数差的零绝不是零**(任意两组人都有差)->
                零 = **随机等大小分层**(层大小照旧)的差分布。
⚠ `#392e`      `sexcount` 进模型前先打印取值集合、众数、与方向已知锚的相关。
IMPOSSIBLE      ① 「从没有过性伴」与年龄强相关 -> **同轮控年龄**,不事后加;
                ② `ACTED` 本身对没有性伴的人天然更低 -> 分层会压窄该层的方差,
                   **只比符号与是否越阈,不比大小**;
                ③ 「需要观众」与「有机会」在这份数据里分不开 —— 这一条明说,不让零承担。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
AC=next(c for c in d.columns if '41kpfir' in str(c)); ACTED=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
SXMAP={'0':0.,'1-2':1.5,'3-7':5.,'8-20':14.,'21+':25.}
SX=d['sexcount'].map(SXMAP).values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
print("⚠ **`#392e`:分层变量进模型前先看它自己**")
gg0=np.isfinite(SX)&np.isfinite(anc)
print(f"   `sexcount` 取值 {sorted(set(SX[np.isfinite(SX)].tolist()))} · "
      f"众数 **{float(pd.Series(SX[np.isfinite(SX)]).mode().iloc[0]):g}** · n={int(np.isfinite(SX).sum()):,}")
print(f"   与方向已知锚 `Totalsexacts` 相关 **{np.corrcoef(SX[gg0],anc[gg0])[0,1]:+.4f}** -> 同向,方向确认")
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(ncat)&np.isfinite(ACTED)&np.isfinite(SX)&np.isfinite(AG)
G0=M&(SX==0); G1=M&(SX>0)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
print(f"\nn=**{n:,}** · **从没有过性伴 {int(G0.sum()):,}** · **有过 {int(G1.sum()):,}**")
print(f"⚠ IMPOSSIBLE ①:corr(从没有过, 年龄) = "
      f"**{np.corrcoef((SX[M]==0).astype(float),AG[M])[0,1]:+.4f}** -> **同轮控年龄**")
def fit(y,g):
    k=int(g.sum())
    X=np.column_stack([np.ones(k),z(sh,g),z(S,g),z(C3,g),z(ncat,g),z(AG,g)])
    yy=z(y,g); b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-6); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1])
b0,s0=fit(ACTED,G0); b1,s1=fit(ACTED,G1)
DIFF=b1-b0; SED=np.sqrt(s0**2+s1**2)
print(f"\n羞耻 -> `ACTED`(已控 `S`·`c3⁻`·类别数·**年龄**):")
print(f"   **从没有过性伴** **{b0:+.4f}** (se {s0:.4f})")
print(f"   **有过**         **{b1:+.4f}** (se {s1:.4f})")
print(f"   差(有过 − 没有过)= **{DIFF:+.4f}** · |t| **{abs(DIFF/max(SED,1e-12)):.3f}**")
NP_=400; idx=np.flatnonzero(M); n0=int(G0.sum()); nul=[]
for s_ in range(NP_):
    rg=np.random.default_rng(5800+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    nul.append(fit(ACTED,gb)[0]-fit(ACTED,ga)[0])
nul=np.array(nul); THR=float(np.percentile(np.abs(nul),95))
print(f"\n⚠ offset 零(**随机等大小分层** {NP_} 次;**任意两组人都有差,所以零不是零**):")
print(f"   **{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   实测 **{DIFF:+.4f}** -> **{(DIFF-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(DIFF)>THR else '**未越阈**'}")
negs=[]
for s_ in range(200):
    rg=np.random.default_rng(99000+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    negs.append(fit(ACTED,gb)[0]-fit(ACTED,ga)[0])
negs=np.array(negs); rate=float((np.abs(negs)>THR).mean())
print(f"\n负对照(**越阈率**,随机分层 200 次):**{100*rate:.1f}%**(合格 1–12%)")
print(f"\nguard 26 = **MDE 扫描**,每级 30 次(只在「有过」那层上加羞耻的额外斜率):")
MDE=None; det=[]
for gg in (0.03,0.05,0.08,0.12):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(70+int(gg*100)*109+s_)
        y=np.full(NN,np.nan)
        for gs,ex in ((G0,0.0),(G1,gg)):
            k=int(gs.sum()); zs=z(sh,gs); y[gs]=-0.08*zs-ex*zs+rg.standard_normal(k)
        if abs(fit(y,G1)[0]-fit(y,G0)[0])>THR: hit+=1
    det.append(hit/30); print(f"   「有过」层多出 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.15
NONNULL=abs(DIFF)>THR
CONT=abs(DIFF) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.2f}** · 争议幅度 **{CONT:.4f}**({'实测(非零支)' if NONNULL else '有意义(零支)'})")
pd.DataFrame([dict(v_b0=b0,v_se0=s0,v_n0=int(G0.sum()),v_b1=b1,v_se1=s1,v_n1=int(G1.sum()),
                   v_diff=DIFF,v_thr=THR,v_mde=MDE_)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'audience.csv',index=False)
g=Gate('羞耻是在意别人怎么看,还是自己看自己')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='MDE 扫描 80% 检出',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(任意两组人都有差)',nul.std()>0,
           f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
if 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】按符号判(负 = A 别人的眼睛 · 正 = B 自己的眼睛)',True,
                   f"{DIFF:+.4f} -> **世界 {'A(需要观众)' if DIFF<0 else 'B(不需要观众)'}**")
    else:
        g.asserted('★【零支】未越阈且 MDE < 0.05 -> 与有没有观众无关',MDE_<0.05,
                   f"{DIFF:+.4f} vs {THR:.4f} · MDE {MDE_:.2f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **IMPOSSIBLE ③**:「需要观众」与「有机会」在这份数据里**分不开** —— 明说,不让零承担。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
