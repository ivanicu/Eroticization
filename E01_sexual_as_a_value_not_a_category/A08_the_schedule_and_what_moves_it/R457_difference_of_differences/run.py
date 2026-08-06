import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A159 R457 -- 羞耻的压制是不是**特定于行动**:直接检验「差的差」

`#412c`:行动上的分层差 **−0.0618**(越阈)vs 幻想上 **−0.0185**(未越),同一设计、同一分层。
**但那两个数从来没有被相减过 —— 那个差没有自己的零。**
`#369` 的规矩:**换估计量,不换问题。**

主量 = **(行动的分层差) − (幻想的分层差)**,在**同一批人**上算。
**两个结局共享同一批人、同一个分层 -> 噪声大部分抵消,比「分别对零检验」有功率得多。**

预测矩阵:
**差的差显著为负** -> 羞耻的压制**特定于行动** -> **排除 A(自我评价)**,剩 B(观众)/ C(机会);
**差的差 ≈ 0** -> 两个结局上的压制一样 -> **支持 A**,而 B/C 被削弱。

ESTIMAND        同一掩码上,两个结局各自 `结局 ~ 羞耻 + 类别数 + 年龄`,层内拟合;
                `DiD = (b1_act − b0_act) − (b1_ther − b0_ther)`。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;guard 26 **显式传 branch**,**网格一开始就加密**;
                        **秩变换版必须同号**(`#384` 的做法,两个结局量纲不同)。
                【非零支】`DiD` 越阈**且为负** -> 压制特定于行动 -> **排除 A**;
                【零支】未越阈**且 MDE < 0.05** -> 两个结局一样 -> **支持 A**。
⚠ 零的种类     `offset_control`:**差的差的零绝不是零** ——
                两个结局彼此相关,**任意**分层都会给出一个非零的差的差。
                零 = **随机等大小分层**(层大小照旧)后重算 `DiD` 的分布(`lib.nulls`)。
IMPOSSIBLE      ① B(观众)与 C(机会)**仍然分不开**(`#412d`)—— 本轮只可能**排除 A**;
                ② 两个结局量纲不同 -> 都标准化,**只比符号与越阈**;
                ③ 同一批人 -> `DiD` 的零已含两结局的相关,但**不含**它们对羞耻的共同依赖的全部结构。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import controls as null_controls
from scipy.stats import rankdata
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
THC=next(c for c in d.columns if 'vmq8jqw' in str(c)); TH=pd.to_numeric(d[THC],errors='coerce').values.astype(float)
AC=next(c for c in d.columns if '41kpfir' in str(c)); ACT=pd.to_numeric(d[AC],errors='coerce').values.astype(float)
SXMAP={'0':0.,'1-2':1.5,'3-7':5.,'8-20':14.,'21+':25.}
SX=d['sexcount'].map(SXMAP).values.astype(float)
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
AG=d['age'].map(AGE).values.astype(float)
M=(COVB>=4)&np.isfinite(sh)&np.isfinite(TH)&np.isfinite(ACT)&np.isfinite(SX)&np.isfinite(AG)&np.isfinite(ncat)
n=int(M.sum()); G0=M&(SX==0); G1=M&(SX>0)
z=lambda v,g:(v[g]-v[g].mean())/max(v[g].std(),1e-12)
print(f"⚠ **同一批人、同一分层、两个结局**(这正是功率的来源)")
print(f"   n=**{n:,}** · 从没有过 **{int(G0.sum()):,}** · 有过 **{int(G1.sum()):,}**")
print(f"   corr(行动, 治疗性) = **{np.corrcoef(ACT[M],TH[M])[0,1]:+.4f}** -> "
      f"**两结局相关,所以差的差的零绝不是零**")
def fit(y,g):
    k=int(g.sum())
    X=np.column_stack([np.ones(k),z(sh,g),z(ncat,g),z(AG,g)]); yy=z(y,g)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); return float(b[1])
def DiD(ya,yt,g0,g1):
    return (fit(ya,g1)-fit(ya,g0))-(fit(yt,g1)-fit(yt,g0))
D_OBS=DiD(ACT,TH,G0,G1)
dA=fit(ACT,G1)-fit(ACT,G0); dT=fit(TH,G1)-fit(TH,G0)
print(f"\n行动的分层差 **{dA:+.4f}** · 幻想的分层差 **{dT:+.4f}** · **DiD = {D_OBS:+.4f}**")
RA=np.full(NN,np.nan); RT=np.full(NN,np.nan)
j=np.flatnonzero(M)
RA[j]=(rankdata(ACT[j])-rankdata(ACT[j]).mean())/rankdata(ACT[j]).std()
RT[j]=(rankdata(TH[j])-rankdata(TH[j]).mean())/rankdata(TH[j]).std()
D_RANK=DiD(RA,RT,G0,G1)
print(f"   秩变换版 **{D_RANK:+.4f}** -> **{'同号' if D_OBS*D_RANK>0 else '⚠ 变号'}**(`#384` 的做法)")
NP_=400; idx=np.flatnonzero(M); n0=int(G0.sum()); nul=[]
for s_ in range(NP_):
    rg=np.random.default_rng(7600+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    nul.append(DiD(ACT,TH,ga,gb))
nul=np.array(nul); THR=float(np.percentile(np.abs(nul),95))
print(f"\n⚠ offset 零(**随机等大小分层** {NP_} 次;"
      f"**两结局相关 -> 任意分层都会给出非零的差的差**):")
print(f"   **{nul.mean():+.5f} ± {nul.std():.5f}** · |值| 95 分位 **{THR:.5f}**")
print(f"   实测 **{D_OBS:+.4f}** -> **{(D_OBS-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈**' if abs(D_OBS)>THR else '**未越阈**'}")
negs=[]
for s_ in range(200):
    rg=np.random.default_rng(99900+s_); p=rg.permutation(idx)
    ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
    negs.append(DiD(ACT,TH,ga,gb))
negs=np.array(negs); rate=float((np.abs(negs)>THR).mean())
print(f"\n负对照(**越阈率**,随机分层 200 次):**{100*rate:.1f}%**")
print(f"\nguard 26 = **MDE 扫描**(⚠ 网格一开始就加密),每级 30 次:")
MDE=None
for gg in (0.015,0.020,0.025,0.030,0.040,0.060):
    hit=0
    for s_ in range(30):
        rg=np.random.default_rng(10+int(gg*1000)*131+s_)
        ya=np.full(NN,np.nan)
        for gs,ex in ((G0,0.0),(G1,gg)):
            k=int(gs.sum()); ya[gs]=-ex*z(sh,gs)+rg.standard_normal(k)
        yt=np.full(NN,np.nan)
        for gs in (G0,G1):
            k=int(gs.sum()); yt[gs]=rg.standard_normal(k)
        if abs(DiD(ya,yt,G0,G1))>THR: hit+=1
    print(f"   只在行动上、只在「有过」层多出 **{gg:.3f}** -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
    if MDE is None and hit>=24: MDE=gg
MDE_=MDE if MDE else 0.08
NONNULL=abs(D_OBS)>THR
CONT=abs(D_OBS) if NONNULL else 0.05
print(f"   **MDE = {MDE_:.3f}** · 争议幅度 **{CONT:.4f}**({'实测' if NONNULL else '有意义'})")
pd.DataFrame([dict(v_dA=dA,v_dT=dT,v_DiD=D_OBS,v_rank=D_RANK,v_thr=THR,v_mde=MDE_,
                   v_n=n,v_n0=int(G0.sum()),v_n1=int(G1.sum()))]).to_csv(
    pathlib.Path(__file__).parent/'results'/'did.csv',index=False)
g=Gate('羞耻的压制是不是特定于行动')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.asserted('★【两支】秩变换版必须同号(两结局量纲不同)',D_OBS*D_RANK>0,
           f"原始 {D_OBS:+.4f} · 秩 {D_RANK:+.4f}",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26:MDE 扫描 vs 争议幅度',MDE_,CONT,True,what='网格一开始就加密',
    branch='non_null' if NONNULL else 'null',main_quantity='continuous')
g.asserted('★【两支】offset 零非退化(两结局相关)',nul.std()>0,
           f"{nul.mean():+.5f} ± {nul.std():.5f}",kind='control')
if 0.01<=rate<=0.12 and D_OBS*D_RANK>0:
    if NONNULL:
        g.asserted('★【非零支】`DiD` 越阈**且为负** -> 压制特定于行动 -> **排除 A(自我评价)**',
                   D_OBS<0,f"{D_OBS:+.4f} vs 阈 {THR:.4f}")
    else:
        g.asserted('★【零支】未越阈**且 MDE < 0.05** -> 两结局一样 -> **支持 A**',MDE_<0.05,
                   f"{D_OBS:+.4f} vs {THR:.4f} · MDE {MDE_:.3f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **IMPOSSIBLE ①**:B(观众)与 C(机会)**仍然分不开** —— 本轮只可能**排除 A**。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
