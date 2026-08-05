import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #463d noted the fork's two sides are named asymmetrically -- `c1+` has a CONTENT-level name
   (it anchors on a countable quantity, fetish categories, at 0.131) while `c3-` has only a
   STRUCTURAL one. #463's NEXT: does `c3-` have a countable external correlate at all?

Worlds
  A  yes -> both sides get a content-level name and the asymmetry I posited does not exist.
  B  no  -> `c3-` is a quantity with NO countable external correspondent, which would explain
     why five naming attempts have died since #203c, and that sentence could go on the page.

⚠ THE PREMISE OF MY OWN NEXT IS SUSPECT AND IS CHECKED FIRST: #442b's own table reports `c3`
   correlating **+0.46** with the endorsed-sex-act count -- far stronger than `c1`'s 0.131. If
   that reproduces, world B is dead before the round starts and the asymmetry I wrote into the
   NEXT was never there. Reading the record first has caught a wrong signpost of mine three
   times this session (#452a, #459a, #462b); this checks whether it is four.

Anchors: every countable quantity on this page whose direction is fixed by construction --
   endorsed sex acts, fetish categories, banded partner count, block coverage, total picks,
   onset-category count.
Rule: `lib.bounded.anchor_rule`, whose positive controls are MANDATORY (#462a) -- no hand-rolled
   bar here, that is what the tool was built for.
Mask: the block mask, the same as #463, or the two sides are not comparable.
CLOSURE unless world B fires.
"""
import numpy as np, pandas as pd, warnings, json as _json
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show, anchor_rule, anchor_rule_controls

GATE=Gate("R508 does c3 have a countable name")
GATE.asserted("CONTROL the rule's self-checks pass", all(anchor_rule_controls()),
              "3/3", kind="control")
_R371=(ROOT/'E01_sexual_as_a_value_not_a_category/A116_is_a_dimension_one_block'
            /'R371_leave_one_block/run.py').read_text()
exec(_R371.split('"""',2)[2].split('def cor(')[0])
keep=list(range(NB)); Ra,Rb=prof_(A,keep),prof_(B,keep)
C=np.zeros((NB,NB))
for i in keep:
    for j in keep:
        g=np.isfinite(Ra[i])&np.isfinite(Rb[j])&m
        if g.sum()>200: C[i,j]=np.corrcoef(Ra[i][g],Rb[j][g])[0,1]
C=(C+C.T)/2; w,V=np.linalg.eigh(C); o=np.argsort(-w); V=V[:,o]
Rr=(np.where(np.isfinite(Ra),Ra,0)+np.where(np.isfinite(Rb),Rb,0))/2
Rr=np.where(np.isfinite(Ra)|np.isfinite(Rb),Rr,np.nan)
Fm=np.isfinite(Rr); Zm=np.where(Fm,Rr,0.0)
def sc(k):
    nu=(V[:,k][:,None]*Zm).sum(0); de=(Fm*np.abs(V[:,k])[:,None]).sum(0)
    return np.where(de>1e-9,nu/np.maximum(de,1e-9),np.nan)
CAND={'c1':sc(0),'c2':sc(1),'c3':sc(2)}
PICKS=np.zeros(NN); _seen=np.zeros(NN,bool)
for _Mb,_ppl in MB:
    PICKS[_ppl]+=_Mb.sum(1); _seen[_ppl]=True
PICKS=np.where(_seen,PICKS,np.nan)
# #449a: privatise everything built so far BEFORE the next splice rebinds names.
_CAND={k:v.copy() for k,v in CAND.items()}; _PICKS=PICKS.copy(); _m=m.copy()
_R416=(ROOT/'E01_sexual_as_a_value_not_a_category/A137_did_the_rare_things_come_first'
            /'R416_is_it_just_earlier/run.py').read_text()
exec(_R416.split('"""',2)[2].split('def fitb')[0])
_COVB=np.array(COVB,dtype=float).copy(); _NCAT=np.array(ncat,dtype=float).copy()
CAND=_CAND; PICKS=_PICKS; m=_m
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
ANCH={'性行为计数':num('Totalsexacts'),'恋物类别数':num('totalfetishcategory'),
      '性伴数分档':raw['sexcount'].astype(str).map(BANDS).values.astype(float),
      '块覆盖数':_COVB,'总勾选项数':PICKS,'起始类别数':_NCAT}
PC={'pornhabit':num('pornhabit'),'bondageaverage':num('bondageaverage'),
    'extroversion':num('extroversionvariable')}
res,pc = anchor_rule(CAND, ANCH, PC, m)
print(f"正对照(块掩码)best |r|:{ {k:round(v[1],3) for k,v in pc.items()} } -> 全部过杠")
rows=[]
for k,(ok,best,rs) in res.items():
    top=max(rs.items(), key=lambda kv: abs(kv[1]) if kv[1]==kv[1] else -1)
    rows.append(dict(component=k, anchorable=ok, best=best, best_anchor=top[0],
                     **{f'r_{a}':v for a,v in rs.items()}))
T=pd.DataFrame(rows)
show(T[['component','anchorable','best','best_anchor']], HERE/'results/c3_anchors.csv',
     n=4, label="三个成分 × 六个可数锚")
T.to_csv(HERE/'results/full.csv',index=False)
for _,r in T.iterrows():
    print(f"   {r.component}: 最强 **{r.best:.4f}**({r.best_anchor})")

c3=T[T.component=='c3'].iloc[0]; c1=T[T.component=='c1'].iloc[0]
print(f"\n⚠ **我 NEXT 里的前提**:`c3⁻` 只有结构级的名字,`c1⁺` 有内容级的 ——")
print(f"   实测 `c3` 最强 **{c3.best:.4f}**,`c1` 最强 **{c1.best:.4f}** -> "
      f"**{'前提是错的:c3 的可数对应物比 c1 强' if c3.best>c1.best else '前提成立'}**")
GATE.asserted("CONTROL the record's own number reproduces (#442b: c3 vs sex acts ≈ 0.46)",
              abs(abs(c3['r_性行为计数'])-0.46)<0.03,
              f"c3 vs sex acts = {c3['r_性行为计数']:+.4f} vs #442b +0.4600", kind="control")
noname = not bool(c3.anchorable)
GATE.asserted("KILL c3 has NO countable external correlate (world B)", noname,
              f"c3 best |r| = {c3.best:.4f} on {c3.best_anchor}")
verdict = "NO_COUNTABLE_NAME" if noname else "PREMISE_WAS_WRONG"
print(f"\n判决 = **{verdict}**")
_json.dump(dict(verdict=verdict,rows=T.to_dict('records'),
                pc={k:v[1] for k,v in pc.items()}),
           open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
