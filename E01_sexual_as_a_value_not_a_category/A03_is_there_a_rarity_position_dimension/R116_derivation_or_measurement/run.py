import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A10 R01 -- IS THE ITEM MAIN EFFECT A MEASUREMENT, OR IS IT 1+1=2?

A09 established that the ITEM main effect is 3.5-11.7x the person x item interaction. Before that
becomes a claim about a CONTENT CATEGORY, it has to survive the arithmetic trap:

  "Before reporting any number, ask -- could this have come out otherwise? If the answer is no
   because the algebra forces it, the finding is a DERIVATION." (realstat)

The item main effect is the held-out R2 of column means. For a binary matrix that is very nearly a
FUNCTION OF THE PREVALENCE DISPERSION ALONE: with p_j the option prevalences, the between-column
variance is Var(p_j) and the total is E[p(1-p)] + Var(p_j). If I is fully predicted by that ratio,
then A09's headline is a restatement of "option base rates differ", not a discovery, and it says
NOTHING about a content detector -- because base rate and content are the same number here.

Three attacks, in the order of the ladder (cheapest first).

ATTACK 1 (arithmetic).  Regress the measured I on the closed-form prediction
                        I_hat = Var(p_j) / (E[p_j(1-p_j)] + Var(p_j)) across blocks.
                        R2 ~ 1 and slope ~ 1 -> DERIVATION. Residual -> a measurement.
ATTACK 2 (gauge).       Does the ORDERING of prevalence carry information beyond its DISPERSION?
                        Replace each block's prevalences with a rank-matched draw from a smooth
                        parametric family (same dispersion, different shape). If I is unchanged,
                        I is blind to which options are popular -- it measures spread, not content.
ATTACK 3 (artifact).    Presentation-order primacy. The release ALPHABETISES multi-select answers
                        (order-consistency across respondents = 1.0000 on 119 pairs), so display
                        order is destroyed. Under the assumption display == alphabetical -- which
                        cannot be verified here -- test prevalence vs alphabetical rank. A null is
                        admissible only after the positive control fires.

ESTIMAND        (1) R2 and slope of measured I on the closed-form I_hat, across blocks
                (2) |I(real prevalences) - I(dispersion-matched surrogate)|
                (3) Spearman(prevalence, alphabetical rank) per block, vs a planted-primacy control
IDENTIFICATION  (1) and (2) identified. (3) identified ONLY under an unverifiable assumption,
                and is reported as such rather than as a clean null.
SCOPE           the 23 blocks A09/R114 identified.
KILL            (1) R2 >= 0.95 AND slope in [0.9,1.1] -> the item effect is a DERIVATION and A09's
                    headline must be reworded from "content category" to "base-rate dispersion".
                    R2 < 0.95 -> there is residual content and the wording survives.
                (2) |delta| <= 2x seed spread -> I is blind to WHICH options are popular.
POSITIVE CTRL   (3) plant a primacy gradient of known size into the prevalences and require
                    detection; the test must NOT already fire at zero planted gradient.
NEGATIVE CTRL   (2) the surrogate IS the negative control for shape.
NOISE FLOOR     3 masks; spread on every quantity.
MULTIPLICITY    23 blocks x 3 attacks x 3 seeds, reported whole.
IMPOSSIBLE      presentation-order primacy as a CLEAN test -- it would require the survey
                instrument's option order, which this release does not carry. Reported N/A with
                that requirement, not as "planned" and not as a passed check.
"""
import pandas as pd, numpy as np, warnings, hashlib
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')

qm=pd.read_csv('data/derived/multiselect_questions.csv')
lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
MINN=20; RAW={}
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]
    vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=MINN].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    RAW[q.qi]=dict(M=M,opt=opt)
P3=pd.read_csv('E01_sexual_as_a_value_not_a_category/A09_does_the_epoch_title_survive/'
               'R114_fixed_margin_null/results/grid.csv')
d1=P3[P3.K==1].groupby(['q','f']).I.mean().unstack('f')
IDENT=sorted(d1.index[(d1[0.]-d1[5.]).abs()<=0.01])
print(f"identified blocks {len(IDENT)}",flush=True)

MASK=0.15; SEEDS=[11,29,47]

def item_r2(M,seed):
    """held-out R2 of column means alone -- the item main effect, first-in (no Shapley needed:
    attack 1 is about whether THIS quantity is forced by the marginals)."""
    rng=np.random.default_rng(seed); obs=rng.random(M.shape)>=MASK; he=~obs
    T=np.where(obs,M,np.nan); gm=np.nanmean(T)
    cm=np.nanmean(T,axis=0); cm=np.where(np.isnan(cm),gm,cm)
    base=np.mean((M[he]-gm)**2)
    pred=np.broadcast_to(np.clip(cm,0,1)[None,:],M.shape)
    return 1.-np.mean((M[he]-pred[he])**2)/base

def closed_form(p):
    v=np.var(p); return v/(np.mean(p*(1-p))+v)

# ---------------- ATTACK 1: arithmetic ----------------
rows=[]
for q in IDENT:
    M=RAW[q]['M']; p=M.mean(0)
    for sd in SEEDS:
        rows.append(dict(q=q,seed=sd,I=item_r2(M,sd),I_hat=closed_form(p),
                         m=M.shape[1],n=M.shape[0],disp=np.var(p)))
A1=pd.DataFrame(rows); g=A1.groupby('q').agg(I=('I','mean'),Isd=('I','std'),I_hat=('I_hat','first'),
                                             m=('m','first'),n=('n','first'))
b,a=np.polyfit(g.I_hat,g.I,1); pred=a+b*g.I_hat
R2=1-((g.I-pred)**2).sum()/((g.I-g.I.mean())**2).sum()
print("\n=== ATTACK 1 (arithmetic): measured item effect vs its closed form ===")
print(g.assign(pred=pred.round(4),resid=(g.I-pred).round(4)).round(4).sort_values('I').to_string())
print(f"\n  slope {b:.4f}   intercept {a:+.4f}   R2 {R2:.4f}   median |resid| {abs(g.I-pred).median():.4f}")
DERIV=(R2>=0.95) and (0.9<=b<=1.1)
print(f"  -> {'DERIVATION: the item effect is forced by prevalence dispersion' if DERIV else 'MEASUREMENT: residual structure beyond dispersion'}")

# ---------------- ATTACK 2: gauge -- dispersion-matched surrogate ----------------
print("\n=== ATTACK 2 (gauge): same dispersion, different prevalence SHAPE ===")
rows=[]
for q in IDENT:
    M=RAW[q]['M']; p=M.mean(0); mp,vp=p.mean(),p.var()
    for sd in SEEDS:
        rng=np.random.default_rng(4000+sd)
        # a smooth surrogate with the SAME mean and variance but a different shape (beta-like)
        z=np.sort(rng.normal(size=len(p))); z=(z-z.mean())/(z.std()+1e-12)
        ps=np.clip(mp+z*np.sqrt(vp),0.005,0.995)
        ps=ps*np.sqrt(vp/max(ps.var(),1e-12)); ps=np.clip(ps-ps.mean()+mp,0.005,0.995)
        Ms=(rng.random(M.shape)<ps[None,:]).astype(float)
        rows.append(dict(q=q,seed=sd,I_real=item_r2(M,sd),I_sur=item_r2(Ms,sd),
                         disp_real=vp,disp_sur=float(ps.var())))
A2=pd.DataFrame(rows); g2=A2.groupby('q').agg(I_real=('I_real','mean'),I_sur=('I_sur','mean'),
        sd_real=('I_real','std'),sd_sur=('I_sur','std'),dr=('disp_real','mean'),ds=('disp_sur','mean'))
g2['delta']=g2.I_real-g2.I_sur; g2['spread']=np.sqrt(g2.sd_real**2+g2.sd_sur**2)
g2['blind']=g2.delta.abs()<=2*g2.spread
print(g2.round(4).to_string())
print(f"\n  median |delta| {g2.delta.abs().median():.4f}   median 2x spread {2*g2.spread.median():.4f}")
print(f"  blocks where I is BLIND to prevalence shape: {int(g2.blind.sum())}/{len(g2)}")
print(f"  dispersion match check: median |Var_real - Var_sur| = {abs(g2.dr-g2.ds).median():.6f}")

# ---------------- ATTACK 3: primacy, with its positive control ----------------
print("\n=== ATTACK 3 (artifact): prevalence vs alphabetical rank ===")
print("  NOTE: the release alphabetises answers, so DISPLAY order is destroyed. This tests the")
print("  assumption display == alphabetical, which cannot be verified here.")
rows=[]
for q in IDENT:
    M=RAW[q]['M']; p=M.mean(0); r=np.arange(len(p))          # opt is sorted -> alphabetical rank
    rho=spearmanr(p,r).statistic
    for gpl in [0.0,0.05,0.15]:                              # planted primacy gradient, graded
        pg=np.clip(p+gpl*(1-r/max(len(r)-1,1)-0.5),0.005,0.995)
        rows.append(dict(q=q,g=gpl,rho=spearmanr(pg,r).statistic if gpl>0 else rho))
A3=pd.DataFrame(rows); s3=A3.groupby('g').rho.agg(['mean','std','count'])
print(s3.round(4).to_string())
obs=float(s3.loc[0.0,'mean']); sd0=float(s3.loc[0.0,'std'])
fires=abs(float(s3.loc[0.15,'mean']))>abs(obs)+2*sd0
print(f"\n  positive control fires at g=0.15: {'PASS' if fires else 'FAIL -- test is blind, null inadmissible'}")
print(f"  observed rho at g=0: {obs:+.4f} (sd {sd0:.4f}) over {int(s3.loc[0.0,'count'])} blocks")
if fires:
    print(f"  -> {'NO alphabetical-prevalence gradient detected' if abs(obs)<2*sd0/np.sqrt(len(IDENT)) else 'GRADIENT PRESENT -- alarming'}"
          f", UNDER the unverifiable assumption display==alphabetical.")
print("  N/A, with what it would require: the survey instrument's option order.")

print("\n"+"="*78)
print("  CONDITIONAL KILL")
print(f"   attack 1  R2 {R2:.4f}, slope {b:.4f}  -> {'DERIVATION' if DERIV else 'MEASUREMENT'}")
print(f"   attack 2  I blind to prevalence shape in {int(g2.blind.sum())}/{len(g2)} blocks")
if DERIV and g2.blind.sum()>=len(g2)*0.5:
    print("\n   -> A09's headline must be REWORDED. The item main effect is a DERIVATION from the")
    print("      prevalence dispersion, and it is blind to WHICH options are popular. What A09")
    print("      established is that OPTION BASE RATES VARY MORE THAN PEOPLE DO -- true, and not")
    print("      evidence of a content detector, because base rate and content are one number here.")
elif DERIV:
    print("\n   -> PARTIAL: forced by dispersion, but sensitive to prevalence shape. Reword the")
    print("      magnitude claim; the content reading is not yet excluded.")
else:
    print("\n   -> A09's wording SURVIVES: the item effect carries structure beyond dispersion.")
print(f"\nartifact sha1 {hashlib.sha1(pd.concat([A1,A2,A3]).to_csv(index=False).encode()).hexdigest()[:12]}")
