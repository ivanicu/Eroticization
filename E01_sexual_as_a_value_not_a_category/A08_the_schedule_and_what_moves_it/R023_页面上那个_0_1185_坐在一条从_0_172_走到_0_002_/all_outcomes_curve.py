import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A133 R407 -- 那个下降是一般性的纳入效应,还是羞耻特有的

`#360b`:`S ↔ 羞耻` 沿覆盖阈从 **+0.1715** 单调掉到 **−0.0023**;
一个**真值固定**的合成关系在同一批阈上也掉 **0.074**,所以约 43% 是选择,**约 0.10 不是**。
`#361`:「固定人群变块数」这条路走不通(块随人而变)。

**⇒ 换一条更直接的路:不去分解那个下降,而是问它**只发生在羞耻上**吗。**

ESTIMAND        同一条覆盖阈曲线(≥4/6/8/10/12/16),跑**全部 29 个结局**:
                各自的 `corr(S, 结局)` 曲线与**两端差**;报**分布**。
KILL            **若多数结局都有类似下降 -> 一般性的纳入效应,羞耻不特殊,`#360` 的旗要降级;
                若只有羞耻(和少数几个)有 -> 羞耻特有,`#360` 的 0.10 更值得追。**
POSITIVE CTRL   真值固定的合成结局 -> 它的两端差就是**「纯选择」的基线**,
                而每个真实结局要对着**它**读,不是对着零读。
NEGATIVE CTRL   `perm_finite`。
⚠ guard 22     先声明这条曲线有几个点(6)。
⚠ 多重性       29 个结局 -> 报**分布**,不挑最大的(`#309c`)。
IMPOSSIBLE      结局的**方向**不同(有的与 `S` 正相关有的负),所以比的是**|两端差|** 与**相对下降**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A89_where_is_the_non_invariance/R333_gender_referential_split/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def curve(rows')[0])

cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
THS=[4,6,8,10,12,16]
SS={TH:Spos(cov>=TH) for TH in THS}
def curve_of(y):
    out=[]
    for TH in THS:
        S=SS[TH]; m=(cov>=TH)&np.isfinite(S)&np.isfinite(y)
        out.append(float(np.corrcoef(S[m],y[m])[0,1]) if m.sum()>=200 else np.nan)
    return np.array(out)
rows=[]
for nm,y in OUT:
    c=curve_of(y.astype(float))
    if np.isfinite(c).sum()<len(THS): continue
    rows.append(dict(v_out=str(nm)[:38],r4=c[0],r16=c[-1],v_drop=c[0]-c[-1],
                     absdrop=abs(c[0])-abs(c[-1])))
T=pd.DataFrame(rows); check_columns(T,'R407')
T.to_csv(pathlib.Path(__file__).parent/'results'/'all_curves.csv',index=False)
rg=np.random.default_rng(31)
S4=SS[4]; m4=np.isfinite(S4); z4=np.full(NN,np.nan)
z4[m4]=(S4[m4]-np.nanmean(S4[m4]))/np.nanstd(S4[m4])
PC=[]
for t in range(8):
    ys=np.full(NN,np.nan); ys[m4]=0.15*z4[m4]+rg.standard_normal(int(m4.sum()))
    c=curve_of(ys); PC.append(abs(c[0])-abs(c[-1]))
PC=np.array(PC)
SH=T[T.v_out.str.contains('ashamed')].iloc[0]
print(f"29 -> **{len(T)}** 个结局有完整六点曲线 · 阈 {THS}")
print(f"\n**|两端差|**(|阈4| − |阈16|,正 = 收窄):")
print(f"   中位 **{T.absdrop.median():+.4f}** · 均值 {T.absdrop.mean():+.4f} · "
      f"范围 [{T.absdrop.min():+.4f}, {T.absdrop.max():+.4f}] · 为正 **{int((T.absdrop>0).sum())}/{len(T)}**")
print(f"   **纯选择基线**(真值固定 0.15 的合成结局,8 次):**{PC.mean():+.4f} ± {PC.std():.4f}**")
print(f"   超过基线的结局:**{int((T.absdrop>PC.mean()+2*PC.std()).sum())}/{len(T)}**")
print(f"\n★ 羞耻:阈4 **{SH.r4:+.4f}** -> 阈16 **{SH.r16:+.4f}** · |两端差| **{SH.absdrop:+.4f}** · "
      f"在 {len(T)} 个结局里排第 **{int((T.absdrop>SH.absdrop).sum())+1}**")
print(f"   相对纯选择基线:**{(SH.absdrop-PC.mean())/max(PC.std(),1e-9):+.2f}** sd")
print(f"\n|两端差| 最大的四个(⚠ 形状,不是结论):")
for _,x in T.nlargest(4,'absdrop').iterrows():
    print(f"   {x.absdrop:+.4f}  ({x.r4:+.4f} -> {x.r16:+.4f})  {x.v_out}")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
NG=np.array([abs(curve_of(perm_finite(OUT[0][1].astype(float),800+i))[0])-
             abs(curve_of(perm_finite(OUT[0][1].astype(float),800+i))[-1]) for i in range(10)])
print(f"负对照(打乱人):|两端差| **{NG.mean():+.4f} ± {NG.std():.4f}**")
gg=Gate('那个下降是一般性的还是羞耻特有的')
gg.curve_has_enough_points('⚠ guard 22:先声明这条曲线有几个点',THS,min_points=3,what='覆盖阈曲线')
gg.asserted('★ 正对照:真值固定的合成结局给出「纯选择」基线',PC.mean()>0.02,
            f"{PC.mean():+.4f} ± {PC.std():.4f} —— **每个真实结局对着它读,不对着零读**")
gg.asserted('★ 负对照:打乱人后 |两端差| ≈ 0',abs(NG.mean())<0.05,f"{NG.mean():+.4f} ± {NG.std():.4f}")
gg.asserted('★ 注册的 kill:多数结局是否都有类似下降(> 半数超过纯选择基线)',
            (T.absdrop>PC.mean()+2*PC.std()).mean()>0.5,
            f"超过基线的 {int((T.absdrop>PC.mean()+2*PC.std()).sum())}/{len(T)} · "
            f"中位 |两端差| {T.absdrop.median():+.4f} vs 基线 {PC.mean():+.4f}")
gg.asserted('⚠ 羞耻在分布里的位置',True,
            f"|两端差| {SH.absdrop:+.4f},排第 {int((T.absdrop>SH.absdrop).sum())+1}/{len(T)},"
            f"相对基线 {(SH.absdrop-PC.mean())/max(PC.std(),1e-9):+.2f} sd")
gg.asserted('⚠ 多重性:报分布,不挑最大的',True,f"{len(T)} 个结局")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
