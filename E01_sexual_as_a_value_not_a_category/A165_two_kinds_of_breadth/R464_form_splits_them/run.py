import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A165 R464 -- 「对画的/写的东西有反应」是不是真的把两种广度劈开了

`#419d②`:`form` 与**性行为计数**负相关(−0.076 / −0.062),与**恋物类别计数**正相关(+0.095 / +0.035)。
⚠ **而那只是两个未经检验的相关** —— `#419f③` 明说它还没有自己的零。
**它读起来像一句关于人的话,那正是最容易被当成结论的形态。**

两个活着的世界:
**A 真的分裂** -> 控制 `S`·`c3⁻`·羞耻·年龄后,两个偏系数**仍反号且各自越阈**;
**B 只是两种计数的共线残留** -> 至少一个**塌掉**。
⚠ 两个计数彼此相关 **+0.3731**,所以 B 不是稻草人。

ESTIMAND        `Totalsexacts ~ form_i + S + c3⁻ + 羞耻 + 年龄` 与
                `totalfetishcategory ~ form_i + S + c3⁻ + 羞耻 + 年龄`(**结局不进自己的控制集**);
                主量 = **两个偏系数的符号是否相反,且各自越阈**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**,**网格一开始就加密**。
                【非零支】**反号且两个都越阈** -> 世界 A;
                【零支】至少一个未越阈 -> 世界 B,启用 MDE。
⚠ 零的种类     `offset_control`:**两个系数的零绝不是零** ——
                `form` 与两个计数都相关,而两个计数彼此也相关;
                零 = **`lib.nulls.perm_in` 打乱 `form`**(掩码内)后重算两个系数的分布。
⚠ 方向已查     `#419b` 已用家族标定确认 `LIKERT_PM3` 的「值越高 = 越同意」。
IMPOSSIBLE      ① 两个计数量纲不同 -> 都标准化,**只比符号与越阈**;
                ② 「劈开两种广度」是描述,不定因果;③ 全自报,同源方差。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
def num(c): return pd.to_numeric(d[c],errors='coerce').values.astype(float)
AN=num('animated'); WR=num('written')
FORM=np.where(np.isfinite(AN)&np.isfinite(WR),(AN+WR)/2,np.nan)
A1=num('Totalsexacts'); A2=num('totalfetishcategory')
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
M=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(sh)&np.isfinite(FORM)&np.isfinite(A1)&np.isfinite(A2)&np.isfinite(AG)
n=int(M.sum()); z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
print(f"⚠ **`#392e` / `#419b`:方向已查** —— `LIKERT_PM3` 的「值越高 = 越同意」已由家族标定确认。")
print(f"n=**{n:,}** · corr(两个计数) = **{np.corrcoef(A1[M],A2[M])[0,1]:+.4f}** -> **世界 B 不是稻草人**")
CTRL=[S,C3,sh,AG]
def fit(y,f=None):
    f=FORM if f is None else f
    X=np.column_stack([np.ones(n),z(f,M)]+[z(v,M) for v in CTRL]); yy=z(y,M)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(n-6); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[1]),float(se[1])
b1,s1=fit(A1); b2,s2_=fit(A2)
print(f"\n`form_i` 的偏系数(控制 `S`·`c3⁻`·羞耻·年龄):")
print(f"   -> **性行为计数** **{b1:+.4f}** (se {s1:.4f})")
print(f"   -> **恋物类别计数** **{b2:+.4f}** (se {s2_:.4f})")
print(f"   **反号:{np.sign(b1)!=np.sign(b2)}**")
NP_=400
n1=np.array([fit(A1,perm_in(FORM,M,9700+s))[0] for s in range(NP_)])
n2=np.array([fit(A2,perm_in(FORM,M,9700+s))[0] for s in range(NP_)])
T1=float(np.percentile(np.abs(n1),95)); T2=float(np.percentile(np.abs(n2),95))
print(f"\n⚠ offset 零(**`lib.nulls.perm_in` 打乱 `form`**;`form` 与两计数都相关、两计数彼此也相关"
      f" -> **零不是零**):")
print(f"   性行为计数阈 **{T1:.4f}** · 恋物类别计数阈 **{T2:.4f}**")
O1=abs(b1)>T1; O2=abs(b2)>T2
print(f"   -> 性行为 {'**越阈**' if O1 else '**塌掉**'} · 恋物类别 {'**越阈**' if O2 else '**塌掉**'}")
negs=np.array([fit(A1,perm_in(FORM,M,99930+s))[0] for s in range(200)])
rate=float((np.abs(negs)>T1).mean())
print(f"\n负对照(**越阈率**,200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ 网格一开始就加密),每级 30 次:")
MDE=None
for gg in (0.015,0.020,0.025,0.030,0.040,0.060):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(6+int(gg*1000)*163+s_)
        y=np.full(NN,np.nan); y[M]=gg*z(FORM,M)+rg.standard_normal(n)
        if abs(fit(y)[0])>T1: hit+=1
    print(f"   种植 **{gg:.3f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
NONNULL=O1 and O2 and (np.sign(b1)!=np.sign(b2))
CONT=min(abs(b1),abs(b2)) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.3f}** · 争议幅度 **{CONT:.4f}**({'两系数中较小的(实测)' if NONNULL else '有意义'})")
pd.DataFrame([dict(v_out='Totalsexacts',v_b=b1,v_se=s1,v_thr=T1,v_over=bool(O1)),
              dict(v_out='totalfetishcategory',v_b=b2,v_se=s2_,v_thr=T2,v_over=bool(O2))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'split.csv',index=False)
g=Gate('「对画的/写的东西有反应」是不是真的把两种广度劈开了')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='网格一开始就加密',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化',n1.std()>0 and n2.std()>0,
           f"{n1.std():.4f} / {n2.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】**反号且两个都越阈** -> 世界 A(真的分裂)',NONNULL,
               f"性行为 {b1:+.4f}/{T1:.4f} · 恋物类别 {b2:+.4f}/{T2:.4f} · "
               f"反号 {np.sign(b1)!=np.sign(b2)}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
