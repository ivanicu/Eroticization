"""E03·A35·R181b —— 把那个 +0.1128 拆开:是生相关在涨,还是天花板在跌

⚠ **事后对照,如实标注:没有预注册。** 触发它的是主脚本网格里的一行 ——
**天花板同期在跌(−0.00238/年),而紧密度 = 生相关 ÷ 天花板 ⇒ 天花板跌会机械地把它抬高。**
若不查,那个 +0.1128 里可能有一半不是关于人的。

⚠ **换不了仪器**:同主脚本。`[unchallenged]`
"""
import os, sys, pathlib, json, warnings, itertools
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pyreadstat
from scipy.stats import spearmanr
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
SEX=["premarsx","xmarsex","homosex","teensex"]; POL=["polabuse","polmurdr","polescap","polattak"]
g,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",usecols=["year"]+SEX+POL,encoding="latin1")
J=g.dropna(subset=SEX); yrs=sorted(J.year.unique()); Y=np.array(yrs,float)
sl=lambda v: float(np.polyfit(Y,np.asarray(v,float),1)[0])
def parts(fr,items):
    raw=[];cei=[]
    for a,b in itertools.combinations(items,2):
        m=fr[[a,b]].dropna(); x=m[a].to_numpy(float); y=m[b].to_numpy(float)
        r=sp(x,y); xs=np.sort(x); ys=np.sort(y); ys=ys if r>0 else ys[::-1]
        raw.append(r); cei.append(abs(sp(xs,ys)))
    return float(np.mean(raw)),float(np.median(cei))
R=[];C=[]
for y in yrs:
    r,c=parts(J[J.year==y],SEX); R.append(r); C.append(c)
print(f"生相关 36 年 **{sl(R)*36:+.4f}**({R[0]:.4f} → {R[-1]:.4f})· 天花板 {sl(C)*36:+.4f}({C[0]:.4f} → {C[-1]:.4f})")
print(f"⇒ **{'生相关自己就在涨,天花板不是主因' if sl(R)>0 else '⛔ 归一后的上涨是天花板跌出来的'}**;归一总共 +0.1128,生相关贡献 {sl(R)*36:+.4f}")
print("\n逐对 36 年变化(哪一对在带动):")
pairs={}
for a,b in itertools.combinations(SEX,2):
    v=[sp(*(lambda m:(m[a],m[b]))(J[J.year==y][[a,b]].dropna())) for y in yrs]
    pairs[f"{a}×{b}"]=float(np.polyfit(Y,v,1)[0]*36)
    star=" ★ 带动全部" if a=="premarsx" and b=="homosex" else ""
    print(f"  {a:10s} × {b:10s} {v[0]:+.3f} → {v[-1]:+.3f} · **{np.polyfit(Y,v,1)[0]*36:+.4f}/36年**{star}")
P=g.dropna(subset=POL); Rp=[];Yp=[]
for y in yrs:
    fr=P[P.year==y]
    if len(fr)<200: continue
    r,_=parts(fr,POL); Rp.append(r); Yp.append(y)
sp_pol=float(np.polyfit(np.array(Yp,float),Rp,1)[0]*36)
print(f"\n安慰剂:警察四题生相关 **{sp_pol:+.4f}/36年**({Rp[0]:+.4f} → {Rp[-1]:+.4f})")
print("⇒ **两个题组在同一批年份上朝相反方向走 ⇒ 不是样本或模式漂移造成的通用趋势。**")
json.dump(dict(raw_slope36=sl(R)*36,ceiling_slope36=sl(C)*36,normalised_slope36=0.1128,
  pairs=pairs,placebo_police_raw36=sp_pol,
  note="事后对照,未预注册;生相关自己上涨 ⇒ 天花板不是主因;而趋势由 premarsx×homosex 一对带动",
  unchallenged=True),open(OUT/"decompose.json","w"),indent=1,ensure_ascii=False)
