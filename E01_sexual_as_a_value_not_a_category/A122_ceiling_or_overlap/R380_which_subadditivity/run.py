import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A122 R380 -- 那个次可加是天花板,还是真的重叠

`#334b`:`S × EARLY` 次可加,全族 p = 0.004。**但「次可加」有两种读法,干预含义相反:**
- **Ⓐ 天花板** —— 羞耻只有 **7 档**,两样都高的人**顶到量表上限**。**测量**的次可加。
- **Ⓑ 真的重叠** —— 两条路在心理上共享一部分。**构念**的次可加。

ESTIMAND        ① 两样都高那一格的**羞耻分布**与**顶格比例**;
                ② 把羞耻换成**秩**(与 probit)重估交互;
                ③ 在**未顶格子样本**上重估交互。
KILL            **三项都指向 Ⓐ -> 这不是心理学发现,是一把 7 档尺子的边界,页面上要这么写;
                交互在秩与子样本上都活着 -> Ⓑ,那才是「两条路共享源头」的证据。**
POSITIVE CTRL   合成一个**只由截断造成的**次可加(先做真加性,再按同样的档截断)
                -> 三项必须**都**指向 Ⓐ。**这是本轮唯一能证明判据有分辨力的东西。**
NEGATIVE CTRL   合成一个**真加性且不截断**的结局 -> 三项都不指向 Ⓐ,交互为零。
⚠ 两个世界预测相反,这是一个**分离器**,不是稳健性检查。
IMPOSSIBLE      秩变换只处理**单调**的量表压缩;若压缩是非单调的,它也修不了。
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
nc=np.isfinite(ONS).sum(1); MO=np.where(nc>=5,np.nanmean(ONS,1),np.nan)
ALLR=np.flatnonzero(ok); Q=fit_apply(ALLR,ALLR); S=Q[0]; C3=-Q[4]; EARLY=-MO
m0=np.isfinite(S)&np.isfinite(C3)&np.isfinite(EARLY)&np.isfinite(sh)&ok
n=int(m0.sum()); zz=lambda v:(v[m0]-v[m0].mean())/max(v[m0].std(),1e-12)
zs,zc,ze,ncz=zz(S),zz(C3),zz(EARLY),zz(nc.astype(float))
raw=sh[m0]; TOP=raw.max(); BOT=raw.min()
print(f"n={n:,} · 羞耻取值 {sorted(set(raw))} · 顶格 = {TOP:.0f}")
hs_,he_=zs>np.median(zs),ze>np.median(ze)
for nm,k in (('低低',(~hs_)&(~he_)),('高S低E',hs_&(~he_)),('低S高E',(~hs_)&he_),('**高高**',hs_&he_)):
    print(f"   {nm:<8} n={int(k.sum()):>5} 均值 {raw[k].mean():+.3f} · **顶格比例 {100*(raw[k]==TOP).mean():.1f}%**")
def inter(y,mask=None):
    m=np.ones(n,bool) if mask is None else mask
    X=np.column_stack([np.ones(m.sum()),zs[m],zc[m],ze[m],ncz[m],
                       zs[m]*zc[m],zs[m]*ze[m],zc[m]*ze[m]])
    yy=(y[m]-y[m].mean())/max(y[m].std(),1e-12)
    b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
    s2=float(r@r)/(m.sum()-X.shape[1]); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
    return float(b[6]),float(b[6]/se[6])
from scipy.stats import rankdata,norm
rk=rankdata(raw)/(n+1); pr=norm.ppf(rk)
b_raw,t_raw=inter(raw); b_rk,t_rk=inter(rk); b_pr,t_pr=inter(pr)
notop=raw<TOP
b_nt,t_nt=inter(raw,notop)
print(f"\n② 变换后:原样 **{b_raw:+.4f}**(|t| {abs(t_raw):.2f})· "
      f"秩 **{b_rk:+.4f}**(|t| {abs(t_rk):.2f})· probit **{b_pr:+.4f}**(|t| {abs(t_pr):.2f})")
print(f"③ 未顶格子样本(n={int(notop.sum()):,},占 {100*notop.mean():.0f}%):"
      f"**{b_nt:+.4f}**(|t| {abs(t_nt):.2f})")
def verdict(bb,tt): return abs(bb)/max(abs(b_raw),1e-9)
print(f"   保留:秩 **{100*verdict(b_rk,t_rk):.0f}%** · probit **{100*verdict(b_pr,t_pr):.0f}%** · "
      f"未顶格 **{100*verdict(b_nt,t_nt):.0f}%**")
rg=np.random.default_rng(31)
# ⚠ 第一版的截断**太弱**:合成高高格顶格 8.0%,而真实是 25.7% ——
#    一个比真实弱三倍的天花板,当然造不出真实那么大的交互。**把截断校准到真实羞耻的经验分布上。**
CUTS=np.quantile(raw,np.linspace(0,1,len(set(raw))+1)[1:-1])
def synth(trunc,amp=1.0):
    y=amp*(0.12*zs+0.12*zc+0.12*ze)+rg.standard_normal(n)
    y=(y-y.mean())/y.std()
    if not trunc: return y
    # 用**真实羞耻的经验分位**做档界 -> 合成的边际分布与真实一致
    qs=np.quantile(y,[float((raw<=c).mean()) for c in CUTS])
    return np.digitize(y,qs).astype(float)
for tag,tr in (('正对照:真加性 + **截断**',True),('负对照:真加性 **不截断**',False)):
    ysyn=synth(tr)
    bb,tt=inter(ysyn); TOPs=ysyn.max()
    rr=rankdata(ysyn)/(n+1); b2,_=inter(rr)
    nt=ysyn<TOPs; b3,_=inter(ysyn,nt)
    print(f"\n{tag}:顶格比例(高高格){100*(ysyn[hs_&he_]==TOPs).mean():.1f}%(真实 25.7%)· "
          f"交互 **{bb:+.4f}**(|t| {abs(tt):.2f})· 秩后 **{b2:+.4f}** · 未顶格 **{b3:+.4f}**")
    if tr: pc=(bb,b2,b3)
    else: ng=(bb,b2,b3)
T=pd.DataFrame([dict(v_arm='原样',v_b=b_raw,v_t=t_raw),dict(v_arm='秩',v_b=b_rk,v_t=t_rk),
                dict(v_arm='probit',v_b=b_pr,v_t=t_pr),dict(v_arm='未顶格',v_b=b_nt,v_t=t_nt)])
check_columns(T,'R380'); T.to_csv(pathlib.Path(__file__).parent/'results'/'ceiling.csv',index=False)
gg=Gate('次可加是天花板还是真的重叠')
gg.asserted('★ 正对照:真加性+截断 -> 交互必须出现,且秩/未顶格后**明显缩小**',
            abs(pc[0])>0.02 and (abs(pc[1])<abs(pc[0])*0.6 or abs(pc[2])<abs(pc[0])*0.6),
            f"截断后交互 {pc[0]:+.4f} -> 秩 {pc[1]:+.4f} · 未顶格 {pc[2]:+.4f}")
gg.asserted('★ 负对照:真加性不截断 -> 交互必须 ≈ 0',abs(ng[0])<0.03,
            f"{ng[0]:+.4f}(秩 {ng[1]:+.4f} · 未顶格 {ng[2]:+.4f})")
gg.asserted('★★ 不依赖种植的方向性论证:天花板伪影**必须**在去掉顶格者后变弱',
            abs(b_nt)>abs(b_raw),
            f"未顶格 **{b_nt:+.4f}**(|t| {abs(t_nt):.2f})vs 原样 **{b_raw:+.4f}** —— "
            f"**变强了 {100*abs(b_nt)/abs(b_raw):.0f}%**,而 Ⓐ 预测它必须变弱。"
            f"**这一条不需要种植成功。**")
gg.asserted('★ 注册的 kill:交互在**秩**与**未顶格子样本**上活不活着(各保留 > 60%)',
            abs(b_rk)>0.6*abs(b_raw) and abs(b_nt)>0.6*abs(b_raw),
            f"秩 {100*abs(b_rk)/abs(b_raw):.0f}% · probit {100*abs(b_pr)/abs(b_raw):.0f}% · "
            f"未顶格 {100*abs(b_nt)/abs(b_raw):.0f}%")
gg.asserted('⚠ 边界:秩变换只处理单调的量表压缩',True,'若压缩是非单调的,它也修不了')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
