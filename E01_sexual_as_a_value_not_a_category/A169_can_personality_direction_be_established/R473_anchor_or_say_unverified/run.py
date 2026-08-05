import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: three signed personality readings are on the public page; #428c established the numeric
   direction of those variables was never checked. Can it be established at all?

#419b's method (a same-family item whose WORDING locks the direction) is STRUCTURALLY
unavailable: the release carries only pre-computed scale scores, no Big Five items. So the
only remaining anchors are literature priors -- a WEAKER instrument, and it must be labelled
as one on the page.

Anchors, each a count (direction fixed by construction), each paired to the ONE trait it has
a strong, replicated relation to. An anchor establishes ONLY its own trait (#428's NEXT).
  neuroticism   <- TotalMentalIllness   (neuroticism ~ internalising psychopathology, r~.4-.5)
  extroversion  <- sexcount             (extraversion ~ number of partners, r~.15-.25)

Pre-registered BEFORE running:
  POWER GATE : an anchor with |r| < 0.10 CANNOT be used -> write "not anchorable", do not guess.
  KILL       : if no anchor clears the power gate, the three page readings must be marked
               UNVERIFIED on the page -- and that marking is itself required output.
  CONTROL    : the anchors must behave on a variable whose direction IS known --
               `pornhabit` (#428c: +0.181/+0.196 on the count anchors, direction confirmed).
  SIGN RULE  : the literature sign is stated BEFORE the number is read.
FRONTIER: world B forces three page claims to be downgraded.
"""
import pandas as pd, numpy as np
from lib.gates import Gate
from lib.nulls import perm_in

G=Gate("R473 personality direction")
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
# ⚠ #429a: two anchors were BROKEN and neither raised. `TotalMentalIllness` is a string with
# the single value 'Any' (a name saying "Total" that is not a number); `sexcount` is BANDED
# ('0','1-2','3-7','8-20','21+') and to_numeric(errors='coerce') silently kept only the 4,524
# rows that happen to parse -- not a failure, a BIASED SUBSAMPLE.
BANDS={'0':0,'1-2':1,'3-7':2,'8-20':3,'21+':4}
def num(c):
    v=df[c]
    # ⚠ #429b: `v.dtype==object` is FALSE in pandas 3 -- string columns carry dtype `str`.
    # The idiom silently fell through to to_numeric, keeping only the rows whose value is a
    # bare number ('0'), i.e. 4,524 rows ALL EQUAL ZERO -> zero variance -> a NaN correlation
    # that reads exactly like "no power". Test the CONTENT, never the dtype.
    if v.dropna().astype(str).isin(BANDS).mean()>0.9:
        return v.astype(str).map(BANDS)          # ordinal, direction fixed by construction
    return pd.to_numeric(v,errors='coerce')

# literature sign written down BEFORE any number is read
# neuroticism's only candidate anchor (TotalMentalIllness) has ZERO variance in this release
# -> it is STRUCTURALLY not anchorable here, and that is a finding, not a failure.
PAIRS=[("extroversionvariable","sexcount","+",         "外向↑ ⇒ 性伴数↑(r≈.15–.25)")]
NO_ANCHOR=["neuroticismvariable","opennessvariable","consciensiousnessvariable",
           "agreeablenessvariable","powerlessnessvariable"]
CTRL =[("pornhabit","Totalsexacts","+",                "#428c 已确认方向:高 = 用得多")]

rows=[]
for trait,anc,expect,why in PAIRS+CTRL:
    v=num(trait); a=num(anc); m=v.notna()&a.notna()
    r=float(np.corrcoef(v[m],a[m])[0,1])
    nl=[float(np.corrcoef(perm_in(v.values,m.values,seed=6100+i)[m.values],a[m])[0,1])
        for i in range(200)]
    sd=float(np.std(nl))
    rows.append(dict(trait=trait,anchor=anc,expect=expect,n=int(m.sum()),r=r,
                     null_sd=sd,z=r/sd if sd>0 else np.nan,
                     powered=abs(r)>=0.10,
                     implied=("as-coded" if (r>0)==(expect=="+") else "INVERTED"),
                     why=why))
T=pd.DataFrame(rows); T.to_csv(HERE/'results/anchor_power.csv',index=False)
print(f"{'变量':<26}{'锚':<22}{'文献符号':>8}{'r':>9}{'|z|':>7}{'功率':>6}   方向")
for _,x in T.iterrows():
    print(f"{x.trait:<26}{x.anchor:<22}{x.expect:>8}{x.r:>+9.4f}{abs(x.z):>7.1f}"
          f"{('通过' if x.powered else '**不足**'):>8}   {x.implied if x.powered else '—'}")

ctl=T[T.trait=="pornhabit"].iloc[0]
G.asserted("CONTROL the anchor method works on a variable whose direction IS known",
           bool(ctl.powered and ctl.implied=="as-coded"),
           f"pornhabit r = {ctl.r:+.4f}, implied {ctl.implied}", kind="control")

main=T[T.trait!="pornhabit"]
n_powered=int(main.powered.sum())
print(f"\n通过功率闸(|r| ≥ 0.10)的锚 = **{n_powered} / {len(main)}**")
for _,x in main.iterrows():
    if not x.powered:
        print(f"   ⚠ **{x.trait} 不可标定** —— 最强可用锚只有 |r| = {abs(x.r):.4f} < 0.10")

G.asserted("KILL at least one personality direction can be established",
           n_powered>0, f"powered anchors = {n_powered}/{len(main)}")

print(f"\n⚠ 结构上无锚可用(本释放版没有方向由构造固定的对应量)= {len(NO_ANCHOR)} 个:")
for t in NO_ANCHOR: print(f"   · {t}")
verdict = "ANCHORED" if n_powered==len(main) else ("PARTIAL" if n_powered else "NOT_ANCHORABLE")
print(f"\n判决 = {verdict}")
# what does it mean for the page?
if verdict!="ANCHORED":
    un=[x.trait for _,x in main.iterrows() if not x.powered]
    print(f"⇒ 页面上涉及 {un} 的带符号读法必须标 **UNVERIFIED**")
json.dump(dict(verdict=verdict, n_powered=n_powered, n_main=len(main),
               control_ok=bool(ctl.powered and ctl.implied=="as-coded"),
               rows=T.drop(columns=['why']).to_dict('records')),
          open(HERE/'results/verdict.json','w'), indent=1, default=str)
print(G.verdict())
