import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A118 R374 -- 「扣掉流行度的 `c3`」当成一个人层分数

`#328b`:只扣掉流行度后,**块载荷**的图样还在(残差与原载荷相关 +0.9223,两端同一批块)。
**那就把它做成一个人层分数,看它带不带得动那些结论。**

⚠⚠ **跑之前写下的**:`res` 是从 `v0` 减出来的,所以 `c3_np` 与 `c3` **必然**高度相关 ——
**相关本身不是证据**。要看的是**羞耻**那一项掉不掉,以及**跨仪器**那一项掉不掉。

ESTIMAND        用残差载荷 `res` 投影得 `c3_np`,报三项:
                ① `corr(c3_np, c3)`(**已知会高,不作为证据**)
                ② `corr(c3_np, 羞耻)` vs `corr(c3, 羞耻)`
                ③ `c3_np` 与 `#303` 的**起始年龄仪器**上「越轨/普通」两半的相关
                   (**题目集与块仪器不相交** —— 脚本内断言)。
KILL            **若 ② ③ 都与原 `c3` 相当 -> 页面上关于 `c3` 的每一句都可以去掉「可能只是冷门」;
                若 ② 掉下去 -> 羞耻那条路部分经由「冷门」,那要写进页面。**
POSITIVE CTRL   合成一个只由 `res` 驱动的结局 -> 必须被 `c3_np` 抓到且被 `c3` 抓得更弱。
NEGATIVE CTRL   `perm_finite` 打乱人。
⚠ guard 20     |cos| 高而派生相关翻号 = 符号没对齐。
IMPOSSIBLE      `res` 只扣了**流行度**一个属性(`#328a` 的另外两个没扣)——
                所以这是「不是冷门度」的证据,不是「与所有结构属性无关」的证据。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A117_what_does_c3_look_like/R372_block_loadings/run.py').read_text()
exec(_SRC.split('"""',2)[2].split('NBOOT=300')[0])

PREVb=np.array([float(M.mean(0).mean()) for M,_ in MB])
Xp=np.column_stack([np.ones(NB),(PREVb-PREVb.mean())/PREVb.std()])
res=v0-Xp@np.linalg.lstsq(Xp,v0,rcond=None)[0]
res=res/np.linalg.norm(res)
c3_np=score_of(res); c3_ref=score_of(v0)
def cor(u,v,m=None):
    k=np.isfinite(u)&np.isfinite(v)&(ok if m is None else m)
    return (float(np.corrcoef(u[k],v[k])[0,1]),int(k.sum())) if k.sum()>200 else (np.nan,0)
if cor(c3_np,c3_ref)[0]<0: res=-res; c3_np=-c3_np       # ⚠ guard 20:符号对齐到参照
r_cc,_=cor(c3_np,c3_ref)
r_sh_np,n_=cor(c3_np,sh); r_sh_c3,_=cor(c3_ref,sh)
print(f"① `corr(c3_np, c3)` = **{r_cc:+.4f}** ⚠ **已知会高,不作为证据**")
print(f"② ↔羞耻:`c3_np` **{r_sh_np:+.4f}** vs `c3` **{r_sh_c3:+.4f}** · "
      f"保留 **{100*r_sh_np/max(r_sh_c3,1e-9):.1f}%**(n={n_:,})")
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if d[c].map(BIN).notna().sum()>300]
BCOL=[str(q.col) for _,q in keep2.iterrows()][:0]
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
print(f"③ 跨仪器(起始年龄,{NC} 个类别,与块仪器**零列重叠**已断言):")
for nm,v in (('c3_np',c3_np),('c3',c3_ref)):
    a,_=cor(v,trg,okO&ok); b,_=cor(v,com,okO&ok)
    print(f"   {nm:<7} ↔越轨半 **{a:+.4f}** · ↔普通半 **{b:+.4f}** · 比 **{abs(a)/max(abs(b),1e-9):.1f}×**")
    if nm=='c3_np': a_np,b_np=a,b
    else: a_c3,b_c3=a,b
def perm_finite(v,seed):
    z2=v.copy(); j=np.flatnonzero(np.isfinite(z2))
    z2[j]=z2[np.random.default_rng(seed).permutation(j)]; return z2
nul=[cor(perm_finite(c3_np,800+i),sh)[0] for i in range(20)]
print(f"负对照(打乱人)↔羞耻:**{np.mean(nul):+.4f} ± {np.std(nul):.4f}**")
rg=np.random.default_rng(44)
y=np.full(NN,np.nan); g=np.isfinite(c3_np)&ok
y[g]=(c3_np[g]-np.nanmean(c3_np[g]))/np.nanstd(c3_np[g])+rg.standard_normal(int(g.sum()))
p_np,_=cor(c3_np,y); p_c3,_=cor(c3_ref,y)
print(f"\n正对照(只由 `res` 驱动的合成结局):`c3_np` **{p_np:+.4f}** vs `c3` **{p_c3:+.4f}**")
T=pd.DataFrame([dict(v_q='c3_np↔c3',v_val=r_cc),dict(v_q='c3_np↔羞耻',v_val=r_sh_np),
                dict(v_q='c3↔羞耻',v_val=r_sh_c3),dict(v_q='c3_np↔越轨半',v_val=a_np),
                dict(v_q='c3_np↔普通半',v_val=b_np)])
check_columns(T,'R374'); T.to_csv(pathlib.Path(__file__).parent/'results'/'c3_np.csv',index=False)
mde=2.8/np.sqrt(max(n_,1))
gg=Gate('扣掉流行度的 `c3` 带不带得动那些结论')
gg.asserted('★ 正对照:只由 `res` 驱动的结局,`c3_np` 必须比 `c3` 抓得更强',p_np>p_c3,
            f"c3_np {p_np:+.4f} vs c3 {p_c3:+.4f}")
gg.negative_control('★ 负对照:打乱人后 `c3_np` ↔ 羞耻',float(np.mean(nul)),r_sh_np,
    null_spread=float(np.std(nul)),null_kind='`perm_finite` 题内跨人打乱')
gg.sign_flip_needs_direction_change('⚠ guard 20:`c3_np` vs `c3` 的方向一致性与派生量符号',
                                    r_cc,r_sh_c3,r_sh_np)
gg.asserted('★ 注册的 kill ②:羞耻那一项掉不掉(保留 > 80%)',
            r_sh_np/max(r_sh_c3,1e-9)>0.8,
            f"{r_sh_np:+.4f} / {r_sh_c3:+.4f} = {100*r_sh_np/max(r_sh_c3,1e-9):.1f}%")
gg.asserted('★ 注册的 kill ③:跨仪器那一项掉不掉',
            abs(a_np)/max(abs(b_np),1e-9)>3.0,
            f"c3_np 越轨/普通 = {abs(a_np)/max(abs(b_np),1e-9):.1f}× vs c3 的 {abs(a_c3)/max(abs(b_c3),1e-9):.1f}×")
gg.asserted('⚠ 边界:`res` 只扣了流行度一个属性',True,
            '这是「不是冷门度」的证据,不是「与所有结构属性无关」的证据(`#328a` 另外两个没扣)')
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
