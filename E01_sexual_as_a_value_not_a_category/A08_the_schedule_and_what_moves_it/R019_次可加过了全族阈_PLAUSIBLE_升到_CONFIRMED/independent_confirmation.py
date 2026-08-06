import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A141 R424 -- `#379` 的独立确认:换估计量,不换问题;而且**每个门先标它管哪一支**

`#379c`:我的预注册门写错在**开火的那一支**上(拿 MDE 这个「零的判据」去评一个非零结果),
而那是在记下 `#372c①` 的**下一轮**。**⇒ 惯例升级:每个门先标支。本轮是它第一次执行。**

⚠ **先读 `#378e` 的 IMPOSSIBLE 栏**(`#376d`):
① 组内标准化拿不掉**形状**差异(极端作答倾向)—— **本轮用组内秩变换,正是为它**;
② 尺子在合并样本上估 -> 与本轮无关(简模型不含 `S`/`c3⁻`)。

ESTIMAND        `羞耻(组内**秩**变换) ~ EARLY_k + 类别数`,男女各拟合,k ∈ {平均 · 中位 · 最早};
                主量 = 每个 k 的两性系数差与 |t|。
判据(**先标支**)
                【非零支】同号(组 0 更强,`#379a` 预注册方向)**且** |t| 越过自己的族内阈。
                【零支】  仅在结果**未越阈**时启用 MDE:MDE < 0.05 才算「看得见而没有」。
                【两支】  guard 25(口径真的到达人群)· 正/负对照。
POSITIVE CTRL   组 1 多出 0.25 斜率 -> 必须越阈。
NEGATIVE CTRL   打乱性别标签 -> 必须落回零。
⚠ 零的种类     `offset_control`:**随机等大小劈分**,每个 k 各建自己的零。
⚠ 不是三次独立检验  三个 `EARLY` 摘要**彼此高度相关** -> 报的是**一致性**,不是复制。跑前先量相关并印出来。
⚠ 只比符号与越阈  秩变换改了量纲(`#293`:换尺子会换掉被测的东西)-> **不比大小。**
IMPOSSIBLE      秩变换拿掉形状差异,但也拿掉了真实的**幅度**信息;
                「最早起始」对**报得少**的人系统性偏晚 -> 它与类别数的耦合比另外两个强(已同轮控类别数)。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from scipy.stats import rankdata
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first/R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('MINCOV,MINC=')[0])
SEX=pd.to_numeric(d['biomale'],errors='coerce').values.astype(float)
OM=np.where(np.isfinite(O),O,np.nan)
SUM={'平均起始':np.nanmean(OM,1),'中位起始':np.nanmedian(OM,1),'最早起始':np.nanmin(OM,1)}
for k in SUM: SUM[k]=np.where(np.isfinite(O).sum(1)>0,SUM[k],np.nan)
NARROW=ok&np.isfinite(sh)&np.isfinite(SEX)&np.isfinite(ncat)
WIDE=(COVB>=4)&np.isfinite(sh)&np.isfinite(SEX)&np.isfinite(ncat)
for k in SUM: WIDE&=np.isfinite(SUM[k]); NARROW&=np.isfinite(SUM[k])
print(f"n 窄 **{int(NARROW.sum()):,}** -> 宽 **{int(WIDE.sum()):,}**")
KS=list(SUM)
print(f"⚠ 三个摘要**彼此高度相关**(所以报一致性,不是复制):")
for i in range(3):
    for j in range(i+1,3):
        g=WIDE
        print(f"   corr({KS[i]}, {KS[j]}) = **{np.corrcoef(SUM[KS[i]][g],SUM[KS[j]][g])[0,1]:+.4f}**")

def rank_within(y,g):
    out=np.full(NN,np.nan); j=np.flatnonzero(g)
    r=rankdata(y[j]); out[j]=(r-r.mean())/r.std(); return out
def make_fit(pred,mask):
    g0=mask&(SEX==0); g1=mask&(SEX==1)
    yR=np.full(NN,np.nan)
    yR[np.flatnonzero(g0)]=rank_within(sh,g0)[np.flatnonzero(g0)]
    yR[np.flatnonzero(g1)]=rank_within(sh,g1)[np.flatnonzero(g1)]
    def fit(y,g):
        k=int(g.sum()); yy=y[g]; yy=(yy-yy.mean())/max(yy.std(),1e-12)
        X=np.column_stack([np.ones(k),
                           (pred[g]-pred[g].mean())/max(pred[g].std(),1e-12),
                           (ncat[g]-ncat[g].mean())/max(ncat[g].std(),1e-12)])
        b,*_=np.linalg.lstsq(X,yy,rcond=None); r=yy-X@b
        s2=float(r@r)/(k-3); se=np.sqrt(np.diag(s2*np.linalg.pinv(X.T@X)))
        return b[1:],se[1:]
    return fit,g0,g1,yR

print(f"\n三个摘要,宽口径,结局 = **组内秩变换**的羞耻:")
rows=[]; details={}
for k in KS:
    fit,g0,g1,yR=make_fit(SUM[k],WIDE)
    b0,s0=fit(yR,g0); b1,s1=fit(yR,g1)
    diff=b1[0]-b0[0]; t=diff/max(np.sqrt(s0[0]**2+s1[0]**2),1e-12)
    n0=int(g0.sum()); idx=np.flatnonzero(WIDE); mt=[]
    for s_ in range(300):
        rg=np.random.default_rng(6000+s_); p=rg.permutation(idx)
        ga=np.zeros(NN,bool); ga[p[:n0]]=True; gb=np.zeros(NN,bool); gb[p[n0:]]=True
        ba,sa=fit(yR,ga); bb,sb=fit(yR,gb)
        mt.append(float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12)))))
    mt=np.array(mt); thr=float(np.percentile(mt,95))
    rows.append(dict(v_k=k,v_b0=b0[0],v_b1=b1[0],v_diff=diff,v_t=t,v_thr=thr,v_n=int(WIDE.sum())))
    details[k]=(fit,g0,g1,yR,thr,idx,n0)
    print(f"   {k:<8} 组0 **{b0[0]:+.4f}** · 组1 **{b1[0]:+.4f}** · 差 **{diff:+.4f}** · "
          f"|t| **{abs(t):.3f}** vs 自己的阈 **{thr:.3f}** · "
          f"{'**越阈**' if abs(t)>thr else '未越阈'} · {'同号' if diff>0 else '⚠ 变号'}")
T=pd.DataFrame(rows); check_columns(T,'R424')
T.to_csv(pathlib.Path(__file__).parent/'results'/'summaries.csv',index=False)

kk=KS[0]; fit,g0,g1,yR,thr,idx,n0=details[kk]
rg=np.random.default_rng(21); yp=np.full(NN,np.nan)
for gs,ex in ((g0,0.0),(g1,0.25)):
    kx=int(gs.sum()); ze=(SUM[kk][gs]-SUM[kk][gs].mean())/SUM[kk][gs].std()
    yp[gs]=-0.10*ze+ex*ze+rg.standard_normal(kx)
ba,sa=fit(yp,g0); bb,sb=fit(yp,g1)
TP=float(np.max(np.abs((bb-ba)/np.maximum(np.sqrt(sa**2+sb**2),1e-12))))
pp=np.random.default_rng(33).permutation(idx)
ha=np.zeros(NN,bool); ha[pp[:n0]]=True; hb=np.zeros(NN,bool); hb[pp[n0:]]=True
ca,sa2=fit(yR,ha); cb,sb2=fit(yR,hb)
TN=float(np.max(np.abs((cb-ca)/np.maximum(np.sqrt(sa2**2+sb2**2),1e-12))))
print(f"\n正对照(组1 多出 0.25)max-|t| **{TP:.3f}** vs 阈 {thr:.3f} · "
      f"负对照(打乱性别标签)**{TN:.3f}**")
ALLSAME=bool((T.v_diff>0).all()); ALLOVER=bool((T.v_t.abs()>T.v_thr).all())
NOVER=int((T.v_t.abs()>T.v_thr).sum())

g=Gate('#379 的独立确认:换估计量,不换问题')
g.asserted('★【两支】正对照:组1 多出 0.25 -> 必须越阈',TP>thr,f"{TP:.3f} vs {thr:.3f}",kind='control')
g.asserted('★【两支】负对照:打乱性别标签 -> 必须落回零',TN<=thr,f"{TN:.3f} vs {thr:.3f}",kind='control')
g.relaxation_reached_the_population('★【两支】guard 25:口径真的到达人群',
                                    int(NARROW.sum()),int(WIDE.sum()),what='覆盖 ≥8 -> ≥4')
if TP>thr and TN<=thr:
    g.asserted('★【非零支】三个摘要全部同号(预注册方向:组 0 更强)',ALLSAME,
               f"差 {[round(x,4) for x in T.v_diff]}")
    g.asserted('★【非零支】三个摘要全部越过各自的阈',ALLOVER,
               f"越阈 {NOVER}/3 · |t| {[round(abs(x),3) for x in T.v_t]}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\n⚠ **不是三次独立检验**:三个摘要彼此高度相关 -> 上面报的是**一致性**,不是复制。")
print(f"SHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
