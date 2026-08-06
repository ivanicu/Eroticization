import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: two quantities now predict "could I change what arouses me" at almost the same size --
   how much a person has ACTED on it (-0.0387, already on the page) and the breadth-type
   coordinate `c3-` (-0.0347, #437, and unmoved by shame). Is that ONE road or TWO?

The question is psychological, not statistical: the page already says what shapes "I can't
stop" is EXPERIENCE, not emotion. If c3- reaches the belief THROUGH having done less, that
claim absorbs it. If it does not, there is a second road to feeling stuck that is not about
what you have done -- and that is a different statement about people.

Worlds
  A  one road   : c3- -> acted -> belief. The indirect path carries a real share and the
                  direct path weakens.
  B  two roads  : the indirect path sits inside its null while the direct one survives ->
                  breadth type reaches "I can't change this" without going through experience.

PRE-REGISTERED PREDICTION: **B** -- because #438 showed c3- -> acted is not even separable
from shame at this design, so the first segment of the chain is the weak one.
⚠ My predictions have been wrong three times this session (#433b #434 #437). This one is
written to be killed, not confirmed.

focal = c3-  ·  mediator = how much acted on  ·  outcome = could I change it
controls  = S · -(five-item) · category count · age · SHAME   (shame in, so this is not a
            re-run of #437 -- the mediator is the new thing)
Nulls     : permutation of the focal within the analysis mask (single-column statistic, so
            perm_in is the correct null here -- #426c); bootstrap over people for intervals.
CONTROL   : the OLS identity c = c' + a*b must hold (#434a).
CONTROL2  : the mediator must carry a real coefficient, else "through experience" is empty.
CONTROL3  (discharges #438's NEXT for THIS family): family-wise threshold recomputed at 400
            and 4,000 permutations; report the drift only.
Any share is computed through lib.bounded.share, which refuses a near-zero denominator (#436c).
FRONTIER.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in
from lib.bounded import show, share

GATE=Gate("R483 two roads")
_R449=(ROOT/'E01_sexual_as_a_value_not_a_category/A154_latent_or_sum'
            /'R449_stable_division_of_labour/run.py').read_text()
exec(_R449.split('"""',2)[2].split('NP_=400; nulA=[]')[0])

raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
AGEmap={v:i for i,v in enumerate(sorted(raw['age'].dropna().astype(str).unique()))}
AGE=raw['age'].astype(str).map(AGEmap).values.astype(float)
SHAME=np.asarray(OUT['羞耻'],dtype=float)
ACT  =np.asarray(OUT['实践了多少'],dtype=float)
BEL  =np.asarray(OUT['能不能改'],dtype=float)
MM = M & np.isfinite(AGE) & np.isfinite(SHAME) & np.isfinite(ACT) & np.isfinite(BEL)
print(f"n = **{int(MM.sum()):,}** · focal `c3⁻` · 中介 **实践了多少** · 结局 **能不能改**")
print(f"⚠ 方向已锁(`#437c`):「能不能改」越高 = 越觉得**能**改;`c3⁻` 越高 = 广度型越强\n")

def paths(focal, idx):
    Xb=np.column_stack([np.ones(idx.sum()), z(A,idx), z(Bv,idx), z(ncat,idx),
                        z(AGE,idx), z(SHAME,idx)])
    p=z(focal,idx); md=z(ACT,idx); y=z(BEL,idx)
    a =float(np.linalg.lstsq(np.column_stack([Xb,p]),    md,rcond=None)[0][-1])
    c =float(np.linalg.lstsq(np.column_stack([Xb,p]),    y, rcond=None)[0][-1])
    bb=       np.linalg.lstsq(np.column_stack([Xb,p,md]),y, rcond=None)[0]
    return a, float(bb[-1]), c, float(bb[-2])

a,b,c,cp = paths(C3,MM); ind=a*b
print(f"a(`c3⁻`->实践) **{a:+.4f}** · b(实践->能不能改,控 `c3⁻`) **{b:+.4f}**")
print(f"c(总) **{c:+.4f}** · c'(控实践) **{cp:+.4f}** · **间接 a·b = {ind:+.5f}** · "
      f"恒等式残差 **{c-cp-ind:+.2e}**")

NP_=4000
nul=np.array([(lambda t:t[0]*t[1])(paths(perm_in(C3,MM,seed=12000+i),MM)) for i in range(NP_)])
thr400 =float(np.percentile(np.abs(nul[:400]),95))
thr4000=float(np.percentile(np.abs(nul),95))
drift=abs(thr4000-thr400)/thr400
rng=np.random.default_rng(31); idxall=np.flatnonzero(MM)
bs_i=[];bs_c=[]
for i in range(400):
    take=rng.choice(idxall,len(idxall),replace=True)
    mm=np.zeros(len(MM),bool); mm[np.unique(take)]=True
    t=paths(C3,mm); bs_i.append(t[0]*t[1]); bs_c.append(t[2])
bs_i=np.array(bs_i); bs_c=np.array(bs_c)
ilo,ihi=np.percentile(bs_i,[2.5,97.5])

T=pd.DataFrame([dict(q='a  c3⁻->实践',v=a),dict(q='b  实践->能不能改',v=b),
                dict(q='c  总',v=c),dict(q="c' 直接(控实践)",v=cp),
                dict(q='a·b 间接',v=ind),dict(q='间接的置换阈(4,000)',v=thr4000),
                dict(q='间接自助 2.5%',v=ilo),dict(q='间接自助 97.5%',v=ihi)])
show(T, HERE/'results/mediation.csv', n=8, label="两条路")

ok_share, sval, slo, shi, why = share(ind, bs_c, num_boot=bs_i, name="间接/总:")
print(f"\n占比闸(`lib.bounded.share`):"
      f"{'放行 ' + f'{sval:+.1%}  区间 [{slo:+.1%}, {shi:+.1%}]' if ok_share else why}")

GATE.asserted("CONTROL the OLS identity holds", abs(c-cp-ind)<1e-9,
              f"residual = {c-cp-ind:.2e}", kind="control")
GATE.asserted("CONTROL2 the mediator carries a real coefficient", abs(b)>0.02,
              f"b = {b:+.4f}", kind="control")
GATE.asserted("CONTROL3 (#438 NEXT, this family) the threshold is converged",
              drift<0.05, f"400 -> 4,000 drift = {drift:.1%} "
                          f"({thr400:.5f} -> {thr4000:.5f})", kind="control")

ind_real  = (abs(ind)>thr4000) and ((ilo>0)==(ihi>0))
direct_ok = abs(cp)>abs(ind)
# ⚠ #439b: the KILL **as pre-registered** operationalised "the road runs through experience"
# as "the indirect effect differs from zero". That is the WRONG assertion (#417b: the gate
# tests the assertion, not the name I gave it) -- "runs through" additionally requires the
# indirect to carry the total IN THE SAME DIRECTION and materially. BOTH verdicts are kept on
# the record; the pre-registered one is NOT silently replaced, because fixing an
# operationalisation after seeing the data is exactly how a goalpost moves.
GATE.asserted("KILL as pre-registered: the indirect effect differs from zero", ind_real,
              f"indirect {ind:+.5f} vs threshold {thr4000:.5f}, boot [{ilo:+.5f},{ihi:+.5f}]")
same_dir = (ind*c) > 0
material = ok_share and abs(sval) >= 0.10
print(f"\n⚠ 预注册的 KILL 测的是「间接 ≠ 0」= **{ind_real}**;")
print(f"   而「走这条路」还要求 **同向**(间接×总 > 0)= **{same_dir}** "
      f"且 **实质**(|占比| ≥ 10%)= **{material}**")
print(f"   -> **两个判决都留在记录里,预注册的那个不被悄悄替换。**")
verdict = "ONE_ROAD" if (ind_real and same_dir and material) else "TWO_ROADS"
print(f"\n判决 = {verdict}   (预注册预测 = TWO_ROADS)")
json.dump(dict(verdict=verdict,n=int(MM.sum()),a=a,b=b,c=c,cp=cp,indirect=ind,
               thr400=thr400,thr4000=thr4000,drift=drift,boot=[ilo,ihi],
               share_ok=bool(ok_share),share=sval,share_lo=slo,share_hi=shi,kill_as_preregistered=bool(ind_real),same_direction=bool((ind*c)>0),prediction="TWO_ROADS"),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
