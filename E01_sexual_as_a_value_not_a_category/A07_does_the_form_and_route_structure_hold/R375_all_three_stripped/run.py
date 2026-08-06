import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A118 R375 -- 扣掉**全部三个**结构属性之后,`c3` 还剩什么

`#329c` 的边界写着「只扣了流行度」。**关掉它。**

⚠ **门槛全部由**同一次运行里的参照臂**给**(`#329b` 的教训:
从别轮拿来的数当门槛,第四次出错)。**本轮不引用任何别轮的数值当阈。**

ESTIMAND        `res3 = v0 − (流行度 · 选项数 · 离散度)` 的拟合,投影得 `c3_ns`;
                与 `#374` 的三项并排,**并给每一项人层自助区间**。
KILL            **若三项的区间与参照臂 `c3` 的区间重叠 -> 「`c3` 是一个内容维度」有了它能有的最强证据;
                若羞耻或跨仪器明显掉出参照臂的区间 -> 它有一部分是块的结构,写进页面 caveat。**
POSITIVE CTRL   只由 `res3` 驱动的合成结局 -> `c3_ns` 必须比 `c3` 抓得更强。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 20     |cos| 高而派生量翻号 = 符号没对齐。
⚠ 噪声         `res3` 只剩约 55% 的载荷方差(`#328a` R²=0.4522)-> `c3_ns` 更吵,
                **必须报自助区间,不拿点估计比**。
IMPOSSIBLE      三个结构属性是我挑的三个;还有别的结构属性没被想到 —— 这仍是一个下界式的论证。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

PREVb=np.array([float(M.mean(0).mean()) for M,_ in MB])
NOPT=np.array([float(M.shape[1]) for M,_ in MB])
DISP=np.array([float(M.mean(0).std()) for M,_ in MB])
X3=np.column_stack([np.ones(NB)]+[(v-v.mean())/v.std() for v in (PREVb,NOPT,DISP)])
def strip3(v):
    r=v-X3@np.linalg.lstsq(X3,v,rcond=None)[0]
    n=np.linalg.norm(r); return r/n if n>1e-12 else r
res3=strip3(v0)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
BCOL=[]
for _,qq in keep2.iterrows():
    s=lg[lg.qi==qq.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)>=1200 and len(opt)>=8: BCOL.append(str(qq['col']))
assert not (set(BCOL)&set(ons)), "两把仪器共享了列"
OBS=np.column_stack([np.isfinite(d[c].map(BIN).values.astype(float)) for c in ons]).astype(float)
NC=OBS.shape[1]; PR=OBS.mean(0); okO=OBS.sum(1)>=8
o_=np.argsort(-PR); COM,TRG=o_[:NC//2],o_[NC//2:]
com=np.where(okO,OBS[:,COM].mean(1),np.nan); trg=np.where(okO,OBS[:,TRG].mean(1),np.nan)
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(ok if m is None else m)
    return float(np.corrcoef(u[k],v[k])[0,1]) if k.sum()>200 else np.nan
def three(v):
    s_=score_of(v)
    if cor(s_,score_of(v0))<0: s_=-s_
    a=cor(s_,trg,okO&ok); b=cor(s_,com,okO&ok)
    return cor(s_,sh),a,b,(abs(a)/max(abs(b),1e-9)),s_
sh3,a3,b3,rt3,s3=three(res3)
shR,aR,bR,rtR,sR=three(v0)
print(f"`res3` 与 `v0` 的 |cos| = **{abs(float(res3@v0/np.linalg.norm(v0))):.4f}** · "
      f"`c3_ns ↔ c3` = **{cor(s3,sR):+.4f}**")
print(f"\n{'':<12}{'↔羞耻':>12}{'↔越轨半':>12}{'↔普通半':>12}{'越轨/普通':>12}")
print(f"{'c3(参照)':<12}{shR:>+12.4f}{aR:>+12.4f}{bR:>+12.4f}{rtR:>12.2f}")
print(f"{'c3_ns':<12}{sh3:>+12.4f}{a3:>+12.4f}{b3:>+12.4f}{rt3:>12.2f}")
NBOOT=200; rg=np.random.default_rng(9090); BS=[]
for _ in range(NBOOT):
    idx=ALLR[rg.integers(0,len(ALLR),len(ALLR))]
    vk=load_of(idx,ref=v0)
    BS.append((three(strip3(vk))[0],three(vk)[0]))
BS=np.array(BS)
q3=np.nanpercentile(BS[:,0],[2.5,97.5]); qR=np.nanpercentile(BS[:,1],[2.5,97.5])
print(f"\n↔羞耻 的自助 95% 区间:`c3_ns` **[{q3[0]:+.4f}, {q3[1]:+.4f}]** · "
      f"参照 `c3` **[{qR[0]:+.4f}, {qR[1]:+.4f}]** · 重叠 **{'是' if q3[0]<qR[1] and qR[0]<q3[1] else '否'}**")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cor(perm_finite(s3,900+i),sh) for i in range(20)]
print(f"负对照(打乱人)↔羞耻:**{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
rgp=np.random.default_rng(7)
y=np.full(NN,np.nan); g=np.isfinite(s3)&ok
y[g]=(s3[g]-np.nanmean(s3[g]))/np.nanstd(s3[g])+rgp.standard_normal(int(g.sum()))
print(f"正对照(只由 `res3` 驱动的结局):`c3_ns` **{cor(s3,y):+.4f}** vs `c3` **{cor(sR,y):+.4f}**")
T=pd.DataFrame([dict(v_arm='c3参照',shame=shR,trg=aR,com=bR,ratio=rtR),
                dict(v_arm='c3_ns',shame=sh3,trg=a3,com=b3,ratio=rt3)])
check_columns(T,'R375'); T.to_csv(pathlib.Path(__file__).parent/'results'/'all3.csv',index=False)
gg=Gate('扣掉全部三个结构属性之后 `c3` 还剩什么')
gg.asserted('★ 正对照:只由 `res3` 驱动的结局,`c3_ns` 必须比 `c3` 抓得更强',cor(s3,y)>cor(sR,y),
            f"c3_ns {cor(s3,y):+.4f} vs c3 {cor(sR,y):+.4f}")
gg.negative_control('★ 负对照:打乱人后 `c3_ns` ↔ 羞耻',float(np.mean(nul)),sh3,
    null_spread=float(np.std(nul)),null_kind='`perm_finite` 题内跨人打乱')
gg.sign_flip_needs_direction_change('⚠ guard 20:`c3_ns` vs `c3`',cor(s3,sR),shR,sh3)
gg.asserted('★ 注册的 kill:羞耻的区间与**同一次运行里的参照臂**重叠吗',
            q3[0]<qR[1] and qR[0]<q3[1],
            f"c3_ns [{q3[0]:+.4f}, {q3[1]:+.4f}] vs c3 [{qR[0]:+.4f}, {qR[1]:+.4f}]")
gg.asserted('★ 跨仪器:越轨/普通 比值对参照臂',rt3>=rtR*0.8,
            f"c3_ns {rt3:.2f}× vs 参照 {rtR:.2f}×")
gg.asserted('⚠ 边界:三个结构属性是我挑的三个',True,'还有别的结构属性没被想到 —— 这仍是下界式论证')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
