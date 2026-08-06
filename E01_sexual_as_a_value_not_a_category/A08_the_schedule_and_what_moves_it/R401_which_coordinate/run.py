import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A131 R401 -- 那份广度共线,脏在哪一个坐标上

`#356a` 的边界写着「共同性不说是哪一边带着它」。**但那不是不可问的,只是需要一个不同的设计。**

⚠ **一个按构造的预测,写在跑之前**:
`S`(位置分)= 块内**所选项的平均**稀有度 —— **它除以了勾选数**,所以**本该与广度无关**;
`c1/c2/c3` 是块层剖面的**特征向量投影**,**没有除**;
`D` 是两组块的差(两边都除过);清晰度是两半剖面的相关(无量纲)。
**⇒ 预测:脏在 `c1/c2/c3`,`S` 干净。**

ESTIMAND        ① 逐坐标与两个广度代理的相关;
                ② 逐坐标 × 29 结局的**单坐标共同性**(该坐标与广度对该结局的共享成分),
                **报逐坐标的分布**。
KILL            **若共线集中在 `c1/c2/c3` 而 `S` 干净 -> 方向是「特征向量吸收了广度」,
                那是一个构造问题,可修;若 `S` 也脏 -> 除以勾选数没起作用,要回头查 `S` 的定义。**
POSITIVE CTRL   只由广度驱动的合成结局 -> 每个坐标的共同性都必须为正且随其与广度的相关排序。
NEGATIVE CTRL   与广度无关的纯噪声结局 -> 全部 ≈ 0。
⚠ 多重性       6 坐标 × 29 结局 -> **报逐坐标的分布,不报单格**(`#309c`)。
IMPOSSIBLE      共同性仍是记账项;本轮只定位**在哪个坐标上**,不定因果方向。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A89_where_is_the_non_invariance/R333_gender_referential_split/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def curve(rows')[0])

inv=pd.read_csv('data/derived/inventory.csv')
BINo={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
      '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
onsc=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BINo).notna().sum()>300]
ncat=np.column_stack([np.isfinite(d[c].map(BINo).values.astype(float)) for c in onsc]).sum(1).astype(float)
covb=np.zeros(NN)
for M,ppl in MB: covb[ppl]+=1
ALLR=np.flatnonzero(ok); CO0=coords(ALLR); CO=[CO0[0],CO0[1],CO0[2],CO0[3],-CO0[4],CO0[5]]
LAB=['S 位置','D 块间对比','c1','c2','c3⁻','清晰度']
base=ok.copy()
for q_ in CO: base&=np.isfinite(q_)
base&=np.isfinite(ncat)&np.isfinite(covb)
z=lambda v,m:(v[m]-v[m].mean())/max(v[m].std(),1e-12)
print(f"⚠ 跑前的构造预测:`S` 除以了勾选数 -> **本该干净**;`c1/c2/c3` 是特征向量投影 -> **本该脏**。\n")
print(f"① 逐坐标与广度代理的相关(n={int(base.sum()):,}):")
for l,q_ in zip(LAB,CO):
    print(f"   {l:<10} ↔起始类别数 **{np.corrcoef(z(q_,base),z(ncat,base))[0,1]:+.4f}** · "
          f"↔块覆盖数 **{np.corrcoef(z(q_,base),z(covb,base))[0,1]:+.4f}**")
def commonality(y,q_,br):
    m=base&np.isfinite(y)
    if m.sum()<300: return np.nan
    n=int(m.sum()); yy=z(y,m)
    def r2(cols):
        X=np.column_stack([np.ones(n)]+cols); b,*_=np.linalg.lstsq(X,yy,rcond=None)
        r=yy-X@b; return 1-float(r@r)/float(((yy-yy.mean())**2).sum())
    a=r2([z(q_,m)]); b_=r2([z(br,m)]); ab=r2([z(q_,m),z(br,m)])
    return 100*(a+b_-ab)
rows=[]
for nm,y in OUT:
    for l,q_ in zip(LAB,CO):
        v=commonality(y.astype(float),q_,ncat)
        if np.isfinite(v): rows.append(dict(v_out=str(nm)[:36],v_coord=l,v_com=v))
T=pd.DataFrame(rows); check_columns(T,'R401')
T.to_csv(pathlib.Path(__file__).parent/'results'/'percoord.csv',index=False)
print(f"\n② 逐坐标的**单坐标共同性**分布(与起始类别数,29 个结局,单位 pp):")
G=T.groupby('v_coord').v_com
for l in LAB:
    s=G.get_group(l)
    print(f"   {l:<10} 中位 **{s.median():+.4f}** · 均值 {s.mean():+.4f} · "
          f"为正 **{int((s>0).sum())}/{len(s)}** · 最大 {s.max():+.4f}")
med={l:float(G.get_group(l).median()) for l in LAB}
dirty=max(LAB[2:5],key=lambda l:med[l]); print(f"\n★ 最脏的是 **{dirty}**(中位 {med[dirty]:+.4f}pp)· "
      f"`S` 是 **{med['S 位置']:+.4f}pp** · 比 **{med[dirty]/max(abs(med['S 位置']),1e-9):.1f}×**")
rg=np.random.default_rng(19)
m=base; n=int(m.sum())
ypos=np.full(NN,np.nan); ypos[m]=0.4*z(ncat,m)+rg.standard_normal(n)
yneg=np.full(NN,np.nan); yneg[m]=rg.standard_normal(n)
print(f"\n正对照(只由广度驱动):" + ' · '.join(
    f"{l} **{commonality(ypos,q_,ncat):+.3f}**" for l,q_ in zip(LAB,CO)))
print(f"负对照(与广度无关的纯噪声):" + ' · '.join(
    f"{l} **{commonality(yneg,q_,ncat):+.4f}**" for l,q_ in zip(LAB,CO)))
pc=[commonality(ypos,q_,ncat) for q_ in CO]; ng=[commonality(yneg,q_,ncat) for q_ in CO]
gg=Gate('那份广度共线脏在哪个坐标')
gg.asserted('★ 正对照:只由广度驱动 -> 每个坐标的共同性都为正',all(x>0 for x in pc),
            ' · '.join(f"{LAB[i]} {pc[i]:+.3f}" for i in range(6)))
gg.asserted('★ 负对照:与广度无关的纯噪声 -> 全部 ≈ 0',max(abs(x) for x in ng)<0.15,
            ' · '.join(f"{LAB[i]} {ng[i]:+.4f}" for i in range(6)))
gg.asserted('★ 注册的 kill(跑前的构造预测):共线集中在 `c1/c2/c3` 而 `S` 干净',
            med[dirty]>3*abs(med['S 位置']),
            f"最脏 {dirty} 中位 {med[dirty]:+.4f}pp vs `S` {med['S 位置']:+.4f}pp "
            f"({med[dirty]/max(abs(med['S 位置']),1e-9):.1f}×)")
gg.asserted('⚠ 多重性:报逐坐标的分布,不报单格',True,'6 坐标 × 29 结局')
gg.asserted('⚠ 边界:共同性仍是记账项',True,'本轮只定位**在哪个坐标上**,不定因果方向')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
