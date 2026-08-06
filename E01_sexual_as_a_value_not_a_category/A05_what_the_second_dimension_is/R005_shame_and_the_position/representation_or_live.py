import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A44 R247 -- 内容维度是不是一条「表征 vs 实况」的轴

`#201b` 明说不给新名字,但那个东西**可以被检验而不需要发明名字**:
**内容维度对 `animated`/`written` 的预测,与它对"具体他人影像"类结局的预测,方向相反吗?**
若异号 -> 一条「表征 vs 实况」的轴,**直接对上 Ivan 的模型 B(对普通表征的情色赋值)**。

⚠ **分组在看任何相关值之前写下**,只依据**题面**:
    REP  表征:画的 · 写的 · 想象自己是某性别地存在/自慰 —— **没有具体的他人影像**
    LIVE 具体他人影像:伴侣 · 两个人互动 · 支配/被支配 · 已性成熟未成年
    ——   其余(态度题 · 人口学 · 一般特质 · 元认知)不入组
本文件里这段分组注释写在任何 `corr` 之前;`git log -p` 可查它与结果同一次提交且先于输出。

ESTIMAND        `mean(Cres 相关 | REP)` vs `mean(Cres 相关 | LIVE)`,判**是否异号**。
KILL            **若两组均值同号 -> 不是这条轴,`#201b` 那个事后解释一并作废。**
NEGATIVE CTRL   打乱组标签(2000 次)-> 组间差的零分布。
SPECIFICITY     同一分组同时跑 S 与 rho_i —— 若三者都呈现同样的分组差,那不是内容维度的性质。
IMPOSSIBLE      每组只有 6 道题。**n=6 vs 6,功效极低**;本轮只能判**符号**与**是否离开零分布**。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns, check_residualized

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
_,RHO=betas(V)
df_raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in df_raw.columns if df_raw[c].dtype!=object and
     set(pd.Series(df_raw[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and df_raw[c].notna().sum()>10000]

# ---- 分组(题面判定,写在任何相关之前)---------------------------------------
def grp(c):
    L=c.lower()
    if c in ('animated','written'): return 'REP'
    if 'existing (in *nonsexual* situations)' in c or 'masturbating alone as a biological' in c: return 'REP'
    if 'partner' in L or 'two people' in L or 'dominant in sexual' in L or 'submissive in sexual' in L \
       or 'clearly reached sexual maturity' in L: return 'LIVE'
    return '--'
G=np.array([grp(c) for c in lik])
print("分组(题面判定):")
for tag in ('REP','LIVE','--'):
    print(f"  {tag:<5}{int((G==tag).sum())} 道: " + ' | '.join(c[:30] for c,g in zip(lik,G) if g==tag)[:150])

qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
NN=len(df_raw); con=np.zeros(NN); pos=np.zeros(NN); cnt=np.zeros(NN); KB=np.zeros(NN)
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    rr=-np.log(np.clip(M.mean(0),1e-4,1.))
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False))
    con[ppl]+=Z@v[:,-1]; pos[ppl]+=(M@rr)/np.maximum(M.sum(1),1); KB[ppl]+=M.sum(1); cnt[ppl]+=1
ok=cnt>=8
Cb=np.where(ok,con/np.maximum(cnt,1),np.nan); Sb=np.where(ok,pos/np.maximum(cnt,1),np.nan)
KB=np.where(ok,KB,np.nan)
base=np.isfinite(Sb)&np.isfinite(Cb)&np.isfinite(KB)&np.isfinite(RHO)&KEEP; bi=np.flatnonzero(base)
X0=np.c_[np.ones(len(bi)),Sb[bi]]
Cr=np.full(NN,np.nan); Cr[bi]=Cb[bi]-X0@np.linalg.lstsq(X0,Cb[bi],rcond=None)[0]
check_residualized(Cr[bi],Sb[bi],'R247 内容残差')

def vec(x):
    out=[]
    for c in lik:
        y=df_raw[c].values.astype(float); m=np.isfinite(y[bi]); jj=bi[m]
        X=np.c_[np.ones(len(jj)),KB[jj]]
        ry=y[jj]-X@np.linalg.lstsq(X,y[jj],rcond=None)[0]
        rx=x[jj]-X@np.linalg.lstsq(X,x[jj],rcond=None)[0]
        out.append(np.corrcoef(ry,rx)[0,1])
    return np.array(out)
VC=vec(Cr); VS=vec(Sb); VR=vec(RHO)
T=pd.DataFrame(dict(q=[c[:52] for c in lik],grp=G,r_Cres=VC,r_S=VS,r_rho=VR))
check_columns(T,'R247'); T.to_csv(pathlib.Path(__file__).parent/'results'/'groups.csv',index=False)
print(f"\n{'组':<6}{'Cres 均值':>11}{'S 均值':>10}{'rho 均值':>11}{'n':>4}")
rows=[]
for tag in ('REP','LIVE'):
    m=G==tag
    rows.append(dict(grp=tag,n=int(m.sum()),c=float(VC[m].mean()),s=float(VS[m].mean()),r=float(VR[m].mean())))
    print(f"{tag:<6}{VC[m].mean():>+11.4f}{VS[m].mean():>+10.4f}{VR[m].mean():>+11.4f}{int(m.sum()):>4}")
dC=rows[0]['c']-rows[1]['c']; dS=rows[0]['s']-rows[1]['s']; dR=rows[0]['r']-rows[1]['r']
print(f"\n组间差(REP − LIVE):Cres {dC:+.4f} · S {dS:+.4f} · rho {dR:+.4f}")

rng=np.random.default_rng(20260803)
idx=np.flatnonzero(G!='--'); k=int((G=='REP').sum())
null=[]
for _ in range(2000):
    p=rng.permutation(idx); a,b=p[:k],p[k:]
    null.append(VC[a].mean()-VC[b].mean())
null=np.array(null); sdn=float(null.std())
print(f"打乱组标签(2000 次):{null.mean():+.4f} ± {sdn:.4f} -> 真实差 {abs(dC)/sdn:.1f}×")

opp=np.sign(rows[0]['c'])!=np.sign(rows[1]['c'])
g=Gate('内容维度是不是表征 vs 实况的轴')
g.asserted('分组在看相关值之前写下(题面判定,与结果同一次提交)',True,
           f"REP {rows[0]['n']} 道 · LIVE {rows[1]['n']} 道 · 其余 {int((G=='--').sum())} 道不入组")
g.negative_control('打乱组标签',float(abs(null.mean())),abs(dC),null_spread=sdn)
g.resolvable('Cres 的组间差',dC,sdn)
g.asserted('特异性:S 与 rho 不应呈现同样大的分组差',abs(dC)>max(abs(dS),abs(dR)),
           f"Cres {dC:+.4f} vs S {dS:+.4f} · rho {dR:+.4f}")
g.asserted('注册的 kill:两组均值同号 -> 不是这条轴',not opp,
           f"REP {rows[0]['c']:+.4f} · LIVE {rows[1]['c']:+.4f} -> {'异号' if opp else '同号'}")
print(g)
print(f"\n  => {'REPRESENTATION AXIS —— 异号,对上模型 B' if opp else 'NOT THIS AXIS —— 同号,#201b 的事后解释作废'}")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
