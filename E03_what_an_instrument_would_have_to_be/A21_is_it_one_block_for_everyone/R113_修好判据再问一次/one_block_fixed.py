"""E03·A21·R113 —— 修好判据再问一次:性道德的「一块」,紧不紧因人而异吗

**类型:FRONTIER**。`#670` 是 UNVERIFIED,两处判据是我自己设坏的。**`#111c`:第一次修正后重试。**

## 两处修法(`#670` 写死,不是看到结果后改的)
**① 不比「最高层是不是同一层」** —— GSS 年龄最高档 `65+`,NSFG 样本上限 45,**那个判据不可满足**。
   改比 **`ρ(层序, 该层最弱一环)`**:**不要求两具仪器有相同的分层标签。**
**② 安慰剂必须比同一个量** —— 非性道德题**不成块**(最弱一环全负),那个量不是同一件事。
   ⇒ 先找**同格式四档、且验证过确实成块**的一组。

## ② 的验证(已跑,写在这里)
`fefam`/`fechld`/`fepresch` **性别角色三题**,四档,与 NSFG 同为三题。
⚠ `fechld` 题干是「母亲工作**不会**伤害孩子」**含否定,方向与另两道相反** ⇒ 取反定向
(读题干是硬规则①,`#661`/`#662` 记过同型)。**定向后总体最弱一环 = +0.4167 ✅ 成块。**
`fehome` 是**二档**,格式不符,剔除。安乐死/自杀四题成块(+0.6469)**但是二值,格式不符**,不用。

⚠ **而这个安慰剂同时就是对手假说**:
**若同一条梯度也出现在性别角色题上,那就不是关于性的,是关于「这个人回答任何一组态度题有多一致」** ——
一种**作答风格特质**。**这是本轮的 meta-separator。**

⚠ **BASIN**:连续七次下注反对自己喜欢的结果、五次输。**W1 仍是我想要的** ⇒ **仍下注 W2/W3。**

## G1 ESTIMAND
每层内**最弱一环**(天花板归一后最小的那一对);**主量 = `ρ(层序, 最弱一环)`**。
## G2 CONTROLS
**正对照**:两具总体最弱一环复现 **0.4259 / 0.3452**(容差 0.03)。
**安慰剂/对手**:同一套分层作用在**性别角色三题**上的同一个 ρ。
  **这个零该不该是零?** 该 —— 若「一块的紧度」是性特有的,性别角色上不该有同样的梯度 ⇒ `negative_control`。
## G3/G4:GSS 四个分层 × {四档, 十分位} · NSFG 两个分层 × {三/四档, 十分位},全报。
## KILL(条件式)
if 正对照复现 and 安慰剂 ρ 的 |值| < 0.5×主量:
  **两具仪器同一分层的 ρ 同号 且 各自 bootstrap 区间不含零** -> **W1**
  任一具含零 -> **判不了** · 梯度最大层 == 构成差最大层 -> **W3**
else UNVERIFIED
## IMPOSSIBLE(不写 planned)
**NSFG 没有性别角色题** ⇒ **安慰剂只有 GSS 一侧** · NSFG 全为女性、单波、年龄上限 45 ·
**NSFG「大学+」n=389**,若它承担教育判决则该判决降级 · `[unchallenged]`
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, pyreadstat
from itertools import combinations
from scipy.stats import spearmanr
from lib.gates import Gate
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)
SEEDS=[20260806,7,991]; FLOOR=200
def sp(a,b): return float(spearmanr(np.asarray(a,float),np.asarray(b,float)).statistic)
def rmax(a,b,s=1):
    x=np.sort(np.asarray(a,float)); y=np.sort(np.asarray(b,float))
    if s<0: y=y[::-1]
    return sp(x,y)
def weakest(fr,items,byyear=None):
    vals=[]
    for a,b in combinations(items,2):
        per=[]
        for _,g in ([(None,fr)] if byyear is None else fr.groupby(byyear)):
            m=g[[a,b]].dropna()
            if len(m)<FLOOR or m[a].nunique()<2 or m[b].nunique()<2: continue
            r=sp(m[a],m[b])
            if not np.isfinite(r) or r==0: continue
            c=rmax(m[a],m[b],1 if r>0 else -1)
            if np.isfinite(c) and abs(c)>1e-9: per.append(r/abs(c))
        if per: vals.append(float(np.median(per)))
    return min(vals) if vals else np.nan

G=["premarsx","xmarsex","homosex","teensex"]; FE=["fefam","fechld","fepresch"]
df,_=pyreadstat.read_dta("data/external/gss/GSS_stata/gss7224_r3a.dta",
    usecols=["year","age","educ","sex","attend"]+G+FE, encoding="latin1")
d=df[df.year>=1990].copy(); d["fechld"]=-d["fechld"]
print(f"=== 正对照 ===\n  GSS 性四题总体最弱一环 = **{weakest(d,G,'year'):+.4f}**(`#670`: 0.4259)")
print(f"  GSS 性别角色三题(定向后)总体最弱一环 = **{weakest(d,FE,'year'):+.4f}**(安慰剂须成块)")

def grad(fr,items,keyser,byyear=None,q=None):
    """ρ(层序, 该层最弱一环)。q=None 用给定分档;q=k 用 k 分位。"""
    s=fr[keyser] if isinstance(keyser,str) else keyser
    lab = pd.qcut(s,q,labels=False,duplicates="drop") if q else s
    xs,ys=[],[]
    for i,lv in enumerate(sorted(pd.Series(lab).dropna().unique(),key=lambda x:x)):
        sub=fr[lab==lv]
        w=weakest(sub,items,byyear)
        if np.isfinite(w): xs.append(i); ys.append(w)
    if len(xs)<3: return np.nan,0,[]
    return sp(xs,ys), len(xs), list(zip(xs,ys))

d["年龄段"]=pd.cut(d.age,[17,34,49,64,99],labels=False)
d["教育档"]=pd.cut(d.educ,[-1,11,12,15,20],labels=False)
d["礼拜档"]=pd.cut(d.attend,[-1,1,3,5,8],labels=False)
print("\n=== G3/G4:GSS(主 = 性四题 · 对手 = 性别角色三题)===")
res={}
for s,q in [("年龄段",None),("教育档",None),("礼拜档",None),("age",10),("educ",5)]:
    key = s if q is None else s
    a,na,pa = grad(d,G,key,"year",q); b,nb,pb = grad(d,FE,key,"year",q)
    res[f"GSS·{s}{'' if q is None else f'·{q}分位'}"]=dict(main=float(a) if np.isfinite(a) else None,n=na,
                                                          placebo=float(b) if np.isfinite(b) else None,
                                                          pts=[[int(x),float(y)] for x,y in pa])
    print(f"  {s:6s}{'四档' if q is None else f'{q}分位':6s} 层数 {na}  主 ρ = **{a:+.4f}** · 对手(性别角色) ρ = {b:+.4f}")

# NSFG
NSp=pathlib.Path("data/external/nsfg")
pat=re.compile(r'_column\((\d+)\)\s+\w+\s+(\w+)\s+%(\d+)f\s+"([^"]*)"')
LAY={}
for line in open(NSp/"setup"/"2011_2013_FemRespSetup.dct",errors="replace"):
    mm=pat.search(line)
    if mm: LAY[mm.group(2).lower()]=(int(mm.group(1))-1,int(mm.group(3)))
N=["samesex","sxok18","sxok16"]; EXn=["ager","hieduc"]
cols={n:LAY[n] for n in N+EXn if n in LAY}; buf={n:[] for n in cols}
for line in open(NSp/"2011_2013_FemRespData.dat",errors="replace"):
    for n,(s,w) in cols.items():
        v=line[s:s+w].strip(); buf[n].append(float(v) if v not in ("",".") else np.nan)
X=pd.DataFrame(buf)
for c in N: X[c]=np.where(X[c].between(1,5),X[c],np.nan)
print(f"\n  NSFG 性三题总体最弱一环 = **{weakest(X,N):+.4f}**(`#670`: 0.3452)")
X["年龄段"]=pd.cut(X.ager,[14,24,34,45],labels=False); X["教育档"]=pd.cut(X.hieduc,[4,9,10,12,15],labels=False)
print("\n=== G3/G4:NSFG ===")
for s,q in [("年龄段",None),("教育档",None),("ager",5),("hieduc",4)]:
    a,na,pa = grad(X,N,s,None,q)
    res[f"NSFG·{s}{'' if q is None else f'·{q}分位'}"]=dict(main=float(a) if np.isfinite(a) else None,n=na,
                                                            placebo=None,pts=[[int(x),float(y)] for x,y in pa])
    print(f"  {s:8s}{'档' if q is None else f'{q}分位':6s} 层数 {na}  ρ = **{a:+.4f}**")

print("\n=== 跨仪器一致性(修好后的判据:不要求同一层标签)===")
for pair,(gk,nk) in {"年龄":("GSS·年龄段","NSFG·年龄段"),"教育":("GSS·教育档","NSFG·教育档")}.items():
    a=res[gk]["main"]; b=res[nk]["main"]; pl=res[gk]["placebo"]
    same = (a is not None and b is not None) and ((a>0)==(b>0))
    print(f"  {pair}: GSS ρ={a:+.4f}(对手 {pl:+.4f}) · NSFG ρ={b:+.4f} -> **{'同号' if same else '不同号'}**")
json.dump(dict(res=res,unchallenged=True),open(OUT/"one_block_fixed.json","w"),indent=1,ensure_ascii=False)
print(f"\nwrote {OUT/'one_block_fixed.json'}")
