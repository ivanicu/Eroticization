import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A129 R397 -- 那条次可加是集中在一格,还是遍布分布

`#350a` 把 `S × EARLY` 钉死了(跨人群复现,|t| 4.47),**但「次可加」只说了符号,没说机制。**
`#337b` 已排除人群时间表(格层扣掉后仍在)。

ESTIMAND        `S` 与 `EARLY` 各切**三分位** -> **九宫格**:每格的羞耻均值 · n ·
                **加性预测**(行效应 + 列效应,由边际给)· **观测 − 加性**。
KILL            **若差集中在「两样都最高」那一格 -> 它是一小撮人的效应,解释的尺度要改;
                若九格上分布均匀 -> 它是一个真正遍布分布的乘性调节。**
POSITIVE CTRL   合成一个**只在一格**的效应 -> 分解必须挑出那一格。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ 口径         用宽口径(`CALIBER.md` ⑩ 的「可用」栏:同一对象、加分辨率),
                **并报窄口径参照臂**(`#329b`),**且带上「两个口径是不同人群」那句**(`#346b`)。
⚠ MDE          九格每格 n≈900 -> 每格 se≈0.033;**逐格差必须对着它自己的 se 读。**
IMPOSSIBLE      三分位切分本身是一个旋钮;本轮只做一种切法,规格曲线留给下一轮。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
ONS=np.column_stack([d[c].map(BIN).values.astype(float) for c in onsc])
ncat=np.isfinite(ONS).sum(1); MO=np.where(ncat>=5,np.nanmean(ONS,1),np.nan)
cov=np.zeros(NN)
for M,ppl in MB: cov[ppl]+=1
def Spos(mask):
    cv=np.zeros(NN); ps=np.zeros(NN)
    for M,ppl in MB:
        rr=-np.log(np.clip(M.mean(0),1e-4,1.)); nn=M.sum(1)
        v=np.where(nn>0,(M@rr)/np.maximum(nn,1),np.nan); g=np.isfinite(v)
        cv[ppl[g]]+=1; ps[ppl[g]]+=v[g]
    return np.where(mask&(cv>=1),ps/np.maximum(cv,1),np.nan)
def cells(TH,y=None):
    mk=cov>=TH; S=Spos(mk); E=-MO
    m=mk&np.isfinite(S)&np.isfinite(E)&np.isfinite(sh)
    yy=(sh if y is None else y)[m]
    a,b=S[m],E[m]
    qa=np.quantile(a,[1/3,2/3]); qb=np.quantile(b,[1/3,2/3])
    ia=np.digitize(a,qa); ib=np.digitize(b,qb)
    G=np.zeros((3,3)); Nn=np.zeros((3,3)); SE=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            k=(ia==i)&(ib==j); G[i,j]=yy[k].mean(); Nn[i,j]=k.sum()
            SE[i,j]=yy[k].std()/np.sqrt(max(k.sum(),1))
    gm=yy.mean(); row=G.mean(1)-gm; col=G.mean(0)-gm
    ADD=gm+row[:,None]+col[None,:]
    return G,ADD,G-ADD,Nn,SE,int(m.sum())
for TH,tag in ((4,'宽口径 cov>=4'),(8,'窄口径 cov>=8(参照臂)')):
    G,ADD,D,Nn,SE,n=cells(TH)
    print(f"\n【{tag}】n={n:,} · 每格 n {int(Nn.min()):,}–{int(Nn.max()):,} · se 中位 {np.median(SE):.4f}")
    print(f"   羞耻均值(行=S 三分位,列=EARLY 三分位):")
    for i in range(3): print("      "+' '.join(f"{G[i,j]:+.3f}" for j in range(3)))
    print(f"   **观测 − 加性**:")
    for i in range(3): print("      "+' '.join(f"{D[i,j]:+.3f}" for j in range(3)))
    print(f"   最大 |偏差| **{np.abs(D).max():.4f}**(在格 {np.unravel_index(np.abs(D).argmax(),D.shape)})· "
          f"高高格 **{D[2,2]:+.4f}** · 偏差 sd **{D.std():.4f}** vs se 中位 **{np.median(SE):.4f}**")
    if TH==4: W=(G,ADD,D,Nn,SE,n)
    else: Nq=(G,ADD,D,Nn,SE,n)
G,ADD,D,Nn,SE,n=W
rg=np.random.default_rng(70)
mk=cov>=4; S=Spos(mk); E=-MO
m=mk&np.isfinite(S)&np.isfinite(E)&np.isfinite(sh)
a,b=S[m],E[m]; qa=np.quantile(a,[1/3,2/3]); qb=np.quantile(b,[1/3,2/3])
ia=np.digitize(a,qa); ib=np.digitize(b,qb)
ysyn=np.full(NN,np.nan); base=rg.standard_normal(int(m.sum()))
base[(ia==2)&(ib==2)]+=0.5
ysyn[m]=base
Gs,ADDs,Ds,_,_,_=cells(4,y=ysyn)
print(f"\n正对照(只在「两样都最高」那一格 +0.5):最大 |偏差| 在格 "
      f"**{np.unravel_index(np.abs(Ds).argmax(),Ds.shape)}** · 高高格 **{Ds[2,2]:+.4f}** · "
      f"其余格 |偏差| 中位 {np.median(np.abs(np.delete(Ds,8))):.4f}")
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=np.array([cells(4,y=perm_finite(sh,800+i))[2] for i in range(60)])
print(f"负对照(打乱人 60 次):偏差 sd **{nul.std(axis=0).mean():.4f}** · "
      f"高高格 **{nul[:,2,2].mean():+.4f} ± {nul[:,2,2].std():.4f}**")
zhi=(D[2,2]-nul[:,2,2].mean())/max(nul[:,2,2].std(),1e-12)
share=abs(D[2,2])/max(np.abs(D).sum(),1e-12)
print(f"\n★ 高高格偏差 **{D[2,2]:+.4f}** = **{zhi:+.2f}** 个零展布 · "
      f"占九格总 |偏差| 的 **{100*share:.1f}%**(均摊 = 11.1%)")
T=pd.DataFrame([dict(v_cell=f"S{i}E{j}",v_obs=float(G[i,j]),v_dev=float(D[i,j]),v_n=int(Nn[i,j]))
                for i in range(3) for j in range(3)])
check_columns(T,'R397'); T.to_csv(pathlib.Path(__file__).parent/'results'/'nine.csv',index=False)
gg=Gate('那条次可加集中在一格还是遍布分布')
gg.asserted('★ 正对照:只在一格的效应,分解必须挑出那一格',
            np.unravel_index(np.abs(Ds).argmax(),Ds.shape)==(2,2),
            f"最大在 {np.unravel_index(np.abs(Ds).argmax(),Ds.shape)},高高格 {Ds[2,2]:+.4f}")
gg.negative_control('★ 负对照:打乱人后的高高格偏差',float(nul[:,2,2].mean()),float(D[2,2]),
    null_spread=float(nul[:,2,2].std()),null_kind='`perm_finite` 题内跨人打乱')
gg.asserted('★ 注册的 kill:偏差是否集中在「两样都最高」那一格(占比 > 30%)',
            share>0.30,
            f"高高格占九格总 |偏差| 的 **{100*share:.1f}%**(均摊 11.1%)· "
            f"该格 {D[2,2]:+.4f} = {zhi:+.2f} 个零展布")
gg.asserted('⚠ 参照臂:窄口径的图样一不一样',
            np.sign(Nq[2][2,2])==np.sign(D[2,2]),
            f"宽 高高格 {D[2,2]:+.4f} · 窄 {Nq[2][2,2]:+.4f}")
gg.asserted('⚠⚠ 两个口径是不同的人群',True,'`#346b`:低覆盖的人系统上不一样')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
