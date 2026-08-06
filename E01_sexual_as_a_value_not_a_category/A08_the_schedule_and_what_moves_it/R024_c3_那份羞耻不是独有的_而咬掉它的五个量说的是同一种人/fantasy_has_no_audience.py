import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A159 R456 -- 幻想没有观众,也不需要机会:羞耻的压制还在吗

`#411a`:羞耻对**行动**的压制,在从没有过性伴的人身上几乎消失(−0.015 vs −0.077)。
`#411c`:**但那分不清是「有人看着」还是「本来就有机会做」** —— 加多少 n 都改不了。

**⇒ 换一个结局,它同时把两个混淆都拿掉:**
`"Engaging with or **fantasizing** about what arouses me feels therapeutic or healing to me"`
—— **幻想既不需要观众,也不需要机会。**

三个世界,**而这一次它们的预测不同**:
**A 自我评价** —— 压制来自「我怎么看我自己」,幻想同样受它管
   -> `羞耻 → 治疗性` 的分层差**存在**(与 `#411a` 同号);
**B 观众** —— 幻想没有观众 -> 分层差**消失**;
**C 机会** —— 幻想不需要机会 -> 分层差**消失**。
⚠ **B 与 C 预测相同 -> 本轮**不分离**它们,只**排除 A 或排除 B+C**。这一条写在前面。**

ESTIMAND        按 `sexcount`(0 vs ≥1)分层,层内拟合 `治疗性 ~ 羞耻 + 类别数 + 年龄`;
                主量 = **两层羞耻系数的差**。
判据(**先标支**,`#379c`)
                【两支】**guard 25**(宽臂必须真的宽)· 负对照用**越阈率** ·
                        guard 26 **显式传 branch**,且**网格一开始就加密**(`#411b` 的教训)。
                【非零支】差越阈**且与 `#411a` 同号(负)** -> 世界 A(自我评价);
                【零支】未越阈**且 MDE < 0.05** -> 排除 A -> 世界 B 或 C(仍分不开)。
⚠ 零的种类     `offset_control`:**两层系数差的零绝不是零**(任意两组人都有差)->
                零 = **随机等大小分层**(层大小照旧)的差分布。
IMPOSSIBLE      ① B 与 C 本轮不分离(上面已明说);
                ② `corr(羞耻, 治疗性) = +0.0005`(`#401a`)-> **总体接近零**,
                   所以分层差若存在,意味着**零掩盖了两个相反的一半**,那本身是一个更强的主张;
                ③ 「从没有过性伴」与年龄强相关 -> **同轮控年龄**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
S=make_S(8); ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); C3=-Q[4]
THC=next(c for c in d.columns if 'vmq8jqw' in str(c)); TH=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
SXMAP={'0':0.,'1-2':1.5,'3-7':5.,'8-20':14.,'21+':25.}
SX=d['sexcount'].map(SXMAP).values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
anc=pd.to_numeric(d['Totalsexacts'],errors='coerce').values.astype(float)
print("⚠ **`#392e`:结局与分层变量各自先看清楚**")
g0_=np.isfinite(TH)&np.isfinite(anc)
print(f"   `治疗性` 取值 {np.unique(TH[np.isfinite(TH)]).tolist()} · "
      f"众数 **{float(pd.Series(TH[np.isfinite(TH)]).mode().iloc[0]):g}** · "
      f"与锚相关 **{np.corrcoef(TH[g0_],anc[g0_])[0,1]:+.4f}**")
g1_=np.isfinite(SX)&np.isfinite(anc)
print(f"   `sexcount` 众数 **{float(pd.Series(SX[np.isfinite(SX)]).mode().iloc[0]):g}** · "
      f"与锚相关 **{np.corrcoef(SX[g1_],anc[g1_])[0,1]:+.4f}**")
print(f"⚠ IMPOSSIBLE ②:`corr(羞耻, 治疗性)` 总体 = "
      f"**{np.corrcoef(sh[np.isfinite(sh)&np.isfinite(TH)],TH[np.isfinite(sh)&np.isfinite(TH)])[0,1]:+.4f}**"
      f"(`#401a`)-> **总体接近零**\n")
BASE=np.isfinite(sh)&np.isfinite(TH)&np.isfinite(SX)&np.isfinite(AG)&np.isfinite(ncat)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
def cell(name,mask,preds):
    g0=mask&(SX==0); g1=mask&(SX>0)
    def fit(y,g):
        k=int(g.sum())
        X=np.column_stack([np.ones(k)]+[z(v,g) for v in [sh]+preds]); yy=z(y,g)
        b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
        s2=float(r@r)/(k-len(preds)-2); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
        return float(b[1]),float(se[1])
    b0,s0=fit(TH,g0); b1,s1=fit(TH,g1); diff=b1-b0
    idx=np.flatnonzero(mask); n0=int(g0.sum()); nul=[]
    for s_ in range(400):
        rg=np.random.default_rng(7100+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
        nul.append(fit(TH,gb)[0]-fit(TH,ga)[0])
    nul=np.array(nul); thr=float(np.percentile(np.abs(nul),95))
    return dict(name=name,n=int(mask.sum()),n0=n0,n1=int(g1.sum()),b0=b0,s0=s0,b1=b1,s1=s1,
                diff=diff,thr=thr,fit=fit,g0=g0,g1=g1,mask=mask)
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
idx=np.flatnonzero(C['mask']); n0=C['n0']; negs=[]
for s_ in range(200):
    rg=np.random.default_rng(99700+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    negs.append(C['fit'](TH,gb)[0]-C['fit'](TH,ga)[0])
negs=np.array(negs); rate=float((np.abs(negs)>C['thr']).mean())
print(f"\n负对照(**越阈率**,C 格随机分层 200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(C 格,⚠ **网格一开始就加密**,`#411b` 的教训),每级 30 次:")
MDE=None
for gg in (0.020,0.030,0.040,0.045,0.050,0.060,0.080):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(20+int(gg*1000)*127+s_)
        y=np.full(NN,np.nan)
        for gs,ex in ((C['g0'],0.0),(C['g1'],gg)):
            k=int(gs.sum()); zs=z(sh,gs); y[gs]=-ex*zs+rg.standard_normal(k)
        if abs(C['fit'](y,C['g1'])[0]-C['fit'](y,C['g0'])[0])>C['thr']: hit+=1
    print(f"   「有过」层多出 **{gg:.3f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.10
NONNULL=abs(C['diff'])>C['thr']
CONT=abs(C['diff']) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.3f}** · 争议幅度 **{CONT:.4f}**({'实测' if NONNULL else '有意义'})")
pd.DataFrame([{k:v for k,v in X.items() if k in ('name','n','n0','n1','b0','s0','b1','s1','diff','thr')}
              for X in (A,B,C)]).to_csv(pathlib.Path(__file__).parent/'results'/'cells.csv',index=False)
SAMESIGN=C['diff']<0
g=Gate('幻想没有观众也不需要机会,羞耻的压制还在吗')
g.relaxation_reached_the_population('★【两支】guard 25:放宽的口径真的到达了人群',
                                    A['n'],C['n'],what='覆盖 ≥8 -> ≥4;简模型不依赖 S/c3⁻')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='网格一开始就加密',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
if 0.01<=rate<=0.12:
    if NONNULL:
        g.asserted('★【非零支】差越阈**且与 `#411a` 同号(负)** -> 世界 A(自我评价)',
                   SAMESIGN,f"C 差 {C['diff']:+.4f} vs 阈 {C['thr']:.4f} · 同号 {SAMESIGN}")
    else:
        g.asserted('★【零支】未越阈**且 MDE < 0.05** -> **排除 A**,剩 B 或 C(仍分不开)',
                   MDE_<0.05,f"C 差 {C['diff']:+.4f} vs 阈 {C['thr']:.4f} · MDE {MDE_:.3f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **IMPOSSIBLE ①**:B(观众)与 C(机会)在本轮**预测相同** -> 只排除,不分离。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
