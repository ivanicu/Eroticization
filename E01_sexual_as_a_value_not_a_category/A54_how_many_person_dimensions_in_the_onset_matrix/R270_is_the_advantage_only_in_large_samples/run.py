import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A54 R270 -- `u_i` 的优势,是不是只存在于大样本

`#224c`:锋利度分档 —— `u_i` 在全样本上比 `rho_i` 灵敏,在 2,806 人的受限臂上反而更钝。
`#223a` 那张「`u_i` 15/31 vs `rho_i` 10/31」是在**全样本**上算的。
**若优势只存在于大样本,那多打中的 5 道可能全部来自"估得更准",而不是"测到了别的"。**

ESTIMAND        把 31 个结局的面板在**两个随机半样本**上各跑一次,
                判 `u_i` 相对 `rho_i` 的越阈值数之差,在半样本上还剩多少。
KILL            **若半样本上优势消失(平均差 ≤1)-> `#223b` 要收窄成
                「在全样本上更锋利,而这个项目的多数结论本来就跑在全样本上」。**
⚠ 阈值必须逐半重算  半样本 n 更小 -> |r| 的抽样噪声更大 -> 沿用全样本阈值会把噪声读成信号
                (`#207b` 的教训:大样本校准的阈值搬到小样本上会多读)。
NEGATIVE CTRL   每个结局在该半样本内打乱。
NOISE FLOOR     8 组随机对半。
IMPOSSIBLE      半样本上两把刀都变钝 —— 判的是**优势之差**,不是各自的绝对数。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

_SRC=(ROOT/'E01_sexual_as_a_value_not_a_category/A14_is_rare_affinity_a_start_or_a_destination'
          /'R173_does_the_map_radiate_outward'/'run.py').read_text()
exec(_SRC.split('"""',2)[2].split('def plant_u')[0])
V0=V.copy(); rar0=rar.copy(); N,Mc=V0.shape
d=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
lik=[c for c in d.columns if d[c].dtype!=object and
     set(pd.Series(d[c]).dropna().unique())<={-3.,-2.,-1.,0.,1.,2.,3.} and d[c].notna().sum()>10000]
AGE={'14-17':15.5,'18-20':19.0,'21-24':22.5,'25-28':26.5,'29-32':30.5}
ATT={'Significantly less attractive':-3,'Moderately less attractive':-2,'Slightly less attractive':-1,
     'About average attractiveness':0,'Slightly more attractive':1,'Moderately more attractive':2,
     'Significantly more attractive':3}
EXTRA={'age':d['age'].map(AGE),'openness':pd.to_numeric(d['opennessvariable'],errors='coerce'),
 'conscientiousness':pd.to_numeric(d['consciensiousnessvariable'],errors='coerce'),
 'extroversion':pd.to_numeric(d['extroversionvariable'],errors='coerce'),
 'neuroticism':pd.to_numeric(d['neuroticismvariable'],errors='coerce'),
 'agreeableness':pd.to_numeric(d['agreeablenessvariable'],errors='coerce'),
 'powerlessness':pd.to_numeric(d['powerlessnessvariable'],errors='coerce'),
 '关系风格':d['Personally, your preferred relationship style is: (4jib23m)'].map({'Monogamous':0,'Not monogamous':1}),
 '0–14 岁被打屁股':d['From the ages of 0-14, how often were you spanked as a form of discipline? (p957nyk)'].map({'Never':0,'Sometimes':1,'Often':2}),
 '自评吸引力':d['Compared to other people of your same gender and age range, you are (yh6d44s)'].map(ATT),
 '成长期性开放度':d['How "sexually liberated" was your upbringing? (fs700v2)'].map({'Repressed':-1,'Neutral':0,'Liberated':1})}
OUT=[(c,d[c].values.astype(float)) for c in lik]+[(k,v.values.astype(float)) for k,v in EXTRA.items()]

def demean_np(Aa,iters=200,tol=1e-10):
    D=np.where(np.isfinite(Aa),Aa,np.nan)
    for _ in range(iters):
        a=np.nanmean(D,axis=0,keepdims=True); D=D-np.where(np.isfinite(a),a,0)
        b=np.nanmean(D,axis=1,keepdims=True); D=D-np.where(np.isfinite(b),b,0)
        if np.nanmax(np.abs(a))<tol and np.nanmax(np.abs(b))<tol: break
    return D
def scores(rows,seed=0,iters=300):
    keep=np.zeros(N,bool); keep[rows]=True; keep&=(np.isfinite(V0).sum(1)>=8)
    D=demean_np(V0); W=np.isfinite(D)&keep[:,None]; Z=np.where(W,D,0.0)
    k=W.sum(1); rb=np.where(k>0,(W*rar0[None,:]).sum(1)/np.maximum(k,1),0.0)
    Xc=W*(rar0[None,:]-rb[:,None]); yb=np.where(k>0,Z.sum(1)/np.maximum(k,1),0.0); Yc=W*(Z-yb[:,None])
    num=(Yc*Xc).sum(1); den=np.sqrt((Xc*Xc).sum(1))*np.sqrt((Yc*Yc).sum(1))
    rho=np.full(N,np.nan); ok=(k>=8)&(den>1e-12)&keep; rho[ok]=num[ok]/den[ok]
    rng=np.random.default_rng(seed); x=rng.standard_normal(Mc)
    for _ in range(iters):
        Xf=W*x[None,:]; dn=(Xf*Xf).sum(1); u=np.where(dn>1e-12,(Z*Xf).sum(1)/np.maximum(dn,1e-12),0.0)
        Uc=W*u[:,None]; dn2=(Uc*Uc).sum(0); x=np.where(dn2>1e-12,(Z*Uc).sum(0)/np.maximum(dn2,1e-12),0.0)
        n=np.linalg.norm(x)
        if n>0: x=x/n
    Xf=W*x[None,:]; dn=(Xf*Xf).sum(1)
    U=np.full(N,np.nan); good=keep&(dn>1e-12); U[good]=(Z*Xf).sum(1)[good]/dn[good]
    if np.corrcoef(x,rar0)[0,1]<0: U=-U
    return U,rho
def panel(U,R,rng):
    bi=np.flatnonzero(np.isfinite(U)&np.isfinite(R))
    def cr(y,x,ii):
        m=np.isfinite(y[ii])&np.isfinite(x[ii]); jj=ii[m]
        return float(np.corrcoef(y[jj],x[jj])[0,1]) if len(jj)>200 else np.nan
    res={'u':[], 'r':[]}; nl={'u':[], 'r':[]}
    for nm,y in OUT:
        for key,x in (('u',U),('r',R)):
            v=cr(y,x,bi); res[key].append(v)
            ps=[]
            for _ in range(20):
                yp=y.copy(); yp[bi]=rng.permutation(y[bi]); w=cr(yp,x,bi)
                if np.isfinite(w): ps.append(abs(w))
            if len(ps)>=10: nl[key].append(ps)
    out={}
    for key in ('u','r'):
        L=min(len(z) for z in nl[key])
        thr=float(np.nanquantile(np.nanmax(np.array([z[:L] for z in nl[key]]),axis=0),0.95))
        out[key]=(int(np.nansum(np.abs(res[key])>thr)),thr)
    return out

rng=np.random.default_rng(20260803)
allrows=np.flatnonzero(np.isfinite(V0).sum(1)>=8)
Uf,Rf=scores(allrows,seed=1); full=panel(Uf,Rf,rng)
print(f"全样本(n={len(allrows):,}):u_i {full['u'][0]}/31(阈值 {full['u'][1]:.4f}) · "
      f"rho_i {full['r'][0]}/31(阈值 {full['r'][1]:.4f}) · 优势 **{full['u'][0]-full['r'][0]:+d}**")
rows=[]
for s in range(8):
    p=rng.permutation(allrows); h=len(p)//2
    for side,idx in (('A',p[:h]),('B',p[h:])):
        U,R=scores(idx,seed=100+s); r=panel(U,R,rng)
        rows.append(dict(seed=s,side=side,n=len(idx),u_pass=r['u'][0],rho_pass=r['r'][0],
                         adv=r['u'][0]-r['r'][0],thr_u=r['u'][1],thr_r=r['r'][1]))
T=pd.DataFrame(rows); check_columns(T,'R270'); T.to_csv(pathlib.Path(__file__).parent/'results'/'halves.csv',index=False)
print(f"\n半样本(每半约 {int(T.n.mean()):,} 人,16 个半):")
print(f"  u_i 越阈值 {T.u_pass.mean():.1f} ± {T.u_pass.std():.1f} · rho_i {T.rho_pass.mean():.1f} ± {T.rho_pass.std():.1f}")
print(f"  **优势 {T.adv.mean():+.2f} ± {T.adv.std():.2f}**(全样本是 {full['u'][0]-full['r'][0]:+d})")
print(f"  阈值:u_i {T.thr_u.mean():.4f} · rho_i {T.thr_r.mean():.4f}(全样本 {full['u'][1]:.4f} / {full['r'][1]:.4f})")

g=Gate('u_i 的优势是不是只在大样本')
g.asserted('可判前提:全样本复现 `#223a` 的 15 vs 10',abs(full['u'][0]-15)<=2 and abs(full['r'][0]-10)<=2,
           f"{full['u'][0]} vs {full['r'][0]}")
g.asserted('⚠ 阈值已逐半重算(`#207b` 的教训)',T.thr_u.mean()>full['u'][1],
           f"半样本阈值 {T.thr_u.mean():.4f} > 全样本 {full['u'][1]:.4f} —— 小样本噪声更大,阈值确实更高")
g.resolvable('半样本上的优势',float(T.adv.mean()),float(T.adv.std()/np.sqrt(len(T))))
g.offset_control('半样本优势 vs 全样本优势',float(T.adv.mean()),float(full['u'][0]-full['r'][0]),
                 float(T.adv.std()/np.sqrt(len(T))),
                 null_kind='同一面板在全样本上的优势 —— 不是零假设,是"若优势与样本量无关,半样本该落在哪"')
g.asserted('注册的 kill:半样本上优势消失(≤1)-> `#223b` 要收窄',T.adv.mean()<=1,
           f"半样本优势 {T.adv.mean():+.2f} vs 全样本 {full['u'][0]-full['r'][0]:+d}")
print(g)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
