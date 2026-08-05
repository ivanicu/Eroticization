import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A141 R422 -- 三条羞耻的路,对男女是不是同一条

页面上三条路:`S` 0.491pp · `c3⁻` 0.949pp · `EARLY` 0.484pp。
`#309b`/`#311` 问过「有没有东西**缓冲**羞耻」(没有);
**没问过「有没有人身上这三条路本身就不成立」** —— 那是关于**人**的问题,不是关于变量的。

⚠ **先读 `#296` 的 IMPOSSIBLE 栏再设计**(`#376d`):那一轮停在 D5,因为**结局侧混淆**未控 ——
两组的**答题风格**可能不同。`#296` 是在 **29 个结局的曲线**上比,所以控起来很难。
**本轮只有一个结局(羞耻),所以控法是直接的:在**每组内部**标准化羞耻,再拟合。**
这拿掉了量表使用的**均值与方差**差异;**拿不掉形状差异**(如更极端作答)—— 这一条写进 IMPOSSIBLE,不让零承担。

ESTIMAND        `羞耻(组内标准化) ~ S + c3⁻ + EARLY`,男女各拟合一次;
                主量 = **三条路各自的两性系数差**。
KILL(条件式)  仅当对照都过**且 MDE < 0.05** -> 判:**是否有任何一条路的系数差越过族内阈**。
                越阈 -> 那条路对两性不是同一条,页面要写;未越阈 -> **三条路对两性是同一条**。
POSITIVE CTRL   合成一个「`c3⁻` 只在一组里起作用」的结局 -> 必须被抓到。
NEGATIVE CTRL   打乱性别标签 -> 三条差必须都落回零。
⚠ 零的种类     `offset_control`:**两性系数差的零绝不是零** —— 任意两个子样本都会有差。
                零 = **随机等大小劈分**(保住两组大小)后的差分布。
⚠ 多重性       3 条路 -> **族内 max-|t| 阈**(`#334` 的做法),不逐条判。
IMPOSSIBLE      ① 组内标准化拿不掉**形状**差异(极端作答倾向);
                ② `S`/`c3⁻` 的尺子是在**合并样本**上估的(`#293`:换尺子不解决问题,但换尺子会换掉被测的东西)
                   -> 本轮固定用合并尺子,差异因此是「同一把尺子下两组的斜率差」,不是「两组各自的路」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
MINCOV=8; S=make_S(MINCOV)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
EARLY=np.where(np.isfinite(O).sum(1)>0,np.nanmean(np.where(np.isfinite(O),O,np.nan),1),np.nan)
SEX=pd.to_numeric(d['biomale'],errors='coerce').values.astype(float)
mA=ok&np.isfinite(S)&np.isfinite(C3)&np.isfinite(EARLY)&np.isfinite(sh)&np.isfinite(SEX)
n=int(mA.sum()); g0=mA&(SEX==0); g1=mA&(SEX==1)
print(f"n={n:,} · 组 0(非生理男)**{int(g0.sum()):,}** · 组 1(生理男)**{int(g1.sum()):,}**\n")
LAB=['S 位置','c3⁻','EARLY 平均起始']
Xs=[S,C3,EARLY]
def fit(y,g):
    """组内标准化 y(⚠ #296 的结局侧控制,同轮做,不事后加),返回三条系数。"""
    k=int(g.sum()); yy=(y[g]-y[g].mean())/max(y[g].std(),1e-12)
    X=np.column_stack([np.ones(k)]+[(v[g]-v[g].mean())/max(v[g].std(),1e-12) for v in Xs])
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(k-4); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return b[1:4],se[1:4]
b0,s0=fit(sh,g0); b1,s1=fit(sh,g1)
DIFF=b1-b0; SED=np.sqrt(s0**2+s1**2); TT=DIFF/np.maximum(SED,1e-12)
print("三条路,两组各自的系数(羞耻已在**组内**标准化 —— `#296` 的结局侧控制):")
for i,l in enumerate(LAB):
    print(f"   {l:<12} 组0 **{b0[i]:+.4f}** (se {s0[i]:.4f}) · 组1 **{b1[i]:+.4f}** (se {s1[i]:.4f}) · "
          f"差 **{DIFF[i]:+.4f}** · |t| **{abs(TT[i]):.2f}**")

# ---- offset 零:随机等大小劈分 ----
NP_=400; n0=int(g0.sum()); idx=np.flatnonzero(mA)
maxt=[]
for s_ in range(NP_):
    rg=np.random.default_rng(4000+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    ba,sa=fit(sh,ga); bb,sb=fit(sh,gb)
    maxt.append(float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))))
maxt=np.array(maxt); THR=float(np.percentile(maxt,95))
print(f"\n⚠ offset 零(**随机等大小劈分**,{NP_} 次;**任意两个子样本都会有差,所以零不是零**):")
print(f"   族内 max-|t| 阈(95 分位)= **{THR:.3f}** · 零的 max-|t| 均值 {maxt.mean():.3f}")
print(f"   实测 max-|t| = **{np.max(np.abs(TT)):.3f}** ({LAB[int(np.argmax(np.abs(TT)))]}) -> "
      f"{'**越阈**' if np.max(np.abs(TT))>THR else '**未越阈:三条路对两性是同一条**'}")

# ---- MDE:先算再看 ----
print(f"\n先算 MDE 再看数(`#369a`),每级 30 次(只在组 1 上加 `c3⁻` 的额外斜率):")
MDE=None
for gg in (0.03,0.05,0.08,0.12,0.20):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(8000+int(gg*100)*11+s_)
        y=np.full(NN,np.nan)
        for gsel in (g0,g1):
            k=int(gsel.sum()); zc=(C3[gsel]-C3[gsel].mean())/C3[gsel].std()
            y[gsel]=0.12*zc+(gg*zc if gsel is g1 else 0)+rg.standard_normal(k)
        ba,sa=fit(y,g0); bb,sb=fit(y,g1)
        if np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))>THR: hit+=1
    print(f"   组 1 多出 **{gg:+.2f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.25
print(f"   **MDE = {MDE_:.2f}** · 有意义 0.05 -> "
      f"{'**MDE 够小**' if MDE_<=0.05 else '**MDE 不够小 —— 零的内容有限**'}")

# ---- 对照 ----
rg=np.random.default_rng(21); ypos=np.full(NN,np.nan)
for gsel,extra in ((g0,0.0),(g1,0.25)):
    k=int(gsel.sum()); zc=(C3[gsel]-C3[gsel].mean())/C3[gsel].std()
    ypos[gsel]=0.12*zc+extra*zc+rg.standard_normal(k)
ba,sa=fit(ypos,g0); bb,sb=fit(ypos,g1)
TP=float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12))))
rgn=np.random.default_rng(33); pp=rgn.permutation(idx)
ha=np.zeros(NN,bool); ha[pp[:n0]]=True; hb=np.zeros(NN,bool); hb[pp[n0:]]=True
ca,sa2=fit(sh,ha); cb,sb2=fit(sh,hb)
TN=float(np.max(np.abs((cb-ca)/np.maximum(np.sqrt(sa2**2+sb2**2),1e-12))))
print(f"\n正对照(组 1 多出 0.25 的 `c3⁻` 斜率):max-|t| **{TP:.3f}** vs 阈 {THR:.3f}")
print(f"负对照(打乱性别标签):max-|t| **{TN:.3f}** vs 阈 {THR:.3f}")
pd.DataFrame([dict(v_route=LAB[i],v_b0=b0[i],v_b1=b1[i],v_diff=DIFF[i],v_t=TT[i],
                   v_thr=THR,v_mde=MDE_) for i in range(3)]).to_csv(
    pathlib.Path(__file__).parent/'results'/'by_sex.csv',index=False)

g=Gate('三条羞耻的路对男女是不是同一条')
CP=TP>THR; CN=TN<=THR; CO=maxt.std()>0
g.asserted('★ 正对照:组 1 多出 0.25 斜率 -> 必须越阈',CP,f"{TP:.3f} vs {THR:.3f}",kind='control')
g.asserted('★ 负对照:打乱性别标签 -> 必须落回零',CN,f"{TN:.3f} vs {THR:.3f}",kind='control')
g.asserted('★ offset 零非退化(任意两个子样本都会有差)',CO,
           f"零 max-|t| {maxt.mean():.3f} ± {maxt.std():.3f}",kind='control')
if CP and CN and CO:
    g.asserted('★ 注册的 kill:某条路的两性系数差越过族内阈',float(np.max(np.abs(TT)))>THR,
               f"max-|t| {np.max(np.abs(TT)):.3f} vs {THR:.3f}")
    g.null_claim_uses_null_criteria('★ guard 21:这个零可发布吗','NULL',
        perm_quantile=float((maxt>np.max(np.abs(TT))).mean()),mde=MDE_,
        sensitivity_shown=True,meaningful=0.05)
else:
    g.asserted('★ 注册的 kill(对照未过 -> 不判)',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
