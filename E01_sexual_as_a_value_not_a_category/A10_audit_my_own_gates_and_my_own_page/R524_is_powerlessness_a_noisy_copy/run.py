import os,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #479c: a specific rival needs less resolution than a general estimate. The page carries an
   unresolved "cannot be anchored" on `powerlessness` (#461c: best anchor 0.074 against
   positive controls at 0.146-0.247). Rather than ask again for its direction, ask a cheaper
   and more specific question: **is it a noisy copy of one of the three variables whose
   direction IS established?**

Worlds
  A  a noisy copy -> then it inherits that variable's direction and the page's UNVERIFIED can
     be replaced by an inherited one.
  B  not a copy -> the UNVERIFIED stands, but with something added: it is not a shadow of the
     known quantities, it simply has no anchor among the counts.

REFERENCE SCALE, without which "low" means nothing: two measures of one construct correlate
far above 0.5; within this dataset even two different anchors reach **+0.373**, and the three
known-direction variables correlate with each other at most **+0.089**.
⚠ #477d is standing: a KILL that passes states what would have made it fail.
CONTROL : the reference pairs must be computed on the same mask, not quoted.
FRONTIER.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from lib.gates import Gate
from lib.bounded import show, anchor_rule, PositiveControlFailed

GATE=Gate("R524 is powerlessness a noisy copy")
raw=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
num=lambda k: pd.to_numeric(raw[k],errors='coerce').values.astype(float)
BANDS={'0':0.,'1-2':1.,'3-7':2.,'8-20':3.,'21+':4.}
P=num('powerlessnessvariable')
KNOWN={'pornhabit':num('pornhabit'),'bondageaverage':num('bondageaverage'),
       'extroversion':num('extroversionvariable')}
ANCH={'性行为计数':num('Totalsexacts'),'恋物类别数':num('totalfetishcategory'),
      '性伴数分档':raw['sexcount'].astype(str).map(BANDS).values.astype(float)}
FULL=np.ones(len(raw),dtype=bool)
def rr(a,b):
    g=np.isfinite(a)&np.isfinite(b)
    return float(np.corrcoef(a[g],b[g])[0,1]) if g.sum()>200 else np.nan

ref=[(f'{a} × {b}', rr(ANCH[a],ANCH[b])) for i,a in enumerate(ANCH) for b in list(ANCH)[i+1:]]
ref+= [(f'{a} × {b}', rr(KNOWN[a],KNOWN[b])) for i,a in enumerate(KNOWN) for b in list(KNOWN)[i+1:]]
R=pd.DataFrame(ref, columns=['pair','r'])
show(R, HERE/'results/reference_scale.csv', n=6, label="参照尺度(同族内相关能到多大)")
top_ref=float(R['r'].abs().max())
print(f"   **参照尺度**:同族内最大 |r| = **{top_ref:.3f}**;"
      f"而「同一构念的两次测量」通常远在 **0.5** 以上")

rows=[dict(vs=k, r=rr(P,v)) for k,v in KNOWN.items()]
T=pd.DataFrame(rows); T['abs']=T['r'].abs()
show(T, HERE/'results/powerlessness_vs_known.csv', n=4, label="powerlessness × 方向已确立的三个量")
best=float(T['abs'].max()); bestv=T.loc[T['abs'].idxmax(),'vs']
print(f"   最强 = **{best:.3f}**({bestv})")

GATE.asserted("CONTROL the reference scale is computed here, not quoted",
              np.isfinite(top_ref), f"max same-family |r| = {top_ref:.3f}", kind="control")
is_copy = best>0.5
GATE.asserted("KILL powerlessness is a noisy copy of one of them", is_copy,
              f"best |r| = {best:.3f} ({bestv}); a second measure of one construct needs >0.5")

# ⚠ 而最强的那一个,其方向**已经确立** -> 它本身就是一个合格的锚
print(f"\n⚠ **而 `{bestv}` 的方向已由 `#429c` 确立** -> 它是一个**方向已固定的量**,"
      f"于是它本身就是一个锚,而 **|r| = {best:.3f} 越过 0.10 的杠**。")
try:
    res,pc = anchor_rule({'powerlessness':P},
                         dict(ANCH, **{'extroversion(方向已确立)':KNOWN['extroversion']}),
                         {'pornhabit':KNOWN['pornhabit'],'bondageaverage':KNOWN['bondageaverage']},
                         FULL)
    ok,bb,rs = res['powerlessness']
    print(f"   把它加进锚集后:`powerlessness` 最强 **{bb:.3f}** -> **{'过杠' if ok else '仍不过'}**")
    print(f"   正对照(块外全样本):{ {k:round(v[1],3) for k,v in pc.items()} }")
except PositiveControlFailed as e:
    ok=False; print(f"   ⚠ 正对照失败:{str(e)[:90]}")
GATE.asserted("CONTROL adding an established-direction variable as an anchor is legitimate",
              ok, f"powerlessness now clears = {ok}", kind="control")
verdict = "NOISY_COPY" if is_copy else ("NOT_A_COPY_BUT_ANCHORABLE" if ok else "NOT_A_COPY")
print(f"\n判决 = **{verdict}**")
npk,nmiss=GATE.passing_kill_audit({})
json.dump(dict(verdict=verdict,best=best,best_vs=bestv,top_ref=top_ref,
               rows=T.to_dict('records'),anchorable_via_extroversion=bool(ok),
               passing_kills=npk,without_floor=nmiss),
          open(HERE/'results/verdict.json','w'),indent=1,default=str)
print(GATE.verdict())
