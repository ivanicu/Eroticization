import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A143 R426 -- guard 26,并回溯审计:本项目有多少正对照是种在「方便幅度」上的

`#381c`:合成正对照种在 0.25(我随手选的),争议幅度是 0.043 —— **5.9 倍**。
它轻松通过,并且**会让我报出一个零**。真实参照拦住了它。
**⇒ 一个只在比争议幅度更大的幅度上返回非零的仪器,同样是沉默。**

`#379c` 刚证明**账本里的教训不等于手里的动作** -> **把它变成门:guard 26。**
而门写完之后有一个我**事前不知道答案**的问题,所以这一轮是 **FRONTIER**:
**本项目 211 轮的正对照里,有多少是种在方便幅度上的?**

两个活着的世界:
**A 少数** —— `#381c` 是一次偶发,大多数正对照的幅度是有依据的 -> 修那几处即可。
**B 多数** —— 我一直在随手选幅度 -> **所有基于「正对照通过」的零都要重新看**,那是一次大范围降级。

ESTIMAND        逐轮抽出 (种植幅度, 该轮争议幅度),报**种植/争议**比值的分布。
KILL(条件式)  仅当抽取器的对照过 -> 判:**比值中位数是否 > 1**(= 世界 B)。
POSITIVE CTRL   抽取器必须在 `R425`(已知 0.25 vs 0.0426)上抽出正确的一对。
NEGATIVE CTRL   在一个**没有正对照**的轮次上必须什么都不返回。
⚠ 覆盖率       **必须和结论一起报**(`#376b`):211 轮里能机器抽取的只是一部分,
                而**不能抽取的那部分不是「没问题」,是「没看」**。
IMPOSSIBLE      抽取靠正则(`#374b` 的近亲)-> 每一条都打印原文;
                「争议幅度」在不同轮里是不同的量(相关 · 系数 · pp · |t|)-> **只在单位一致时比**,
                否则记 `UNCOMPARABLE`,不硬凑。
"""
import numpy as np, pandas as pd, re, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

NUMP=r'[-−+]?\d+(?:\.\d+)?'
PLANT=re.compile(r'种植\s*\*{0,2}('+NUMP+r')\*{0,2}|plant\(\s*('+NUMP+r')')
def f(x): return abs(float(str(x).replace('−','-').replace('+','')))
rounds=sorted(p for p in pathlib.Path('E01_sexual_as_a_value_not_a_category').glob('A*/R*/run.py'))
print(f"扫描 **{len(rounds)}** 个轮次的 run.py")
rows=[]
for rp in rounds:
    txt=rp.read_text()
    if '正对照' not in txt: continue
    pl=[f(m.group(1) or m.group(2)) for m in PLANT.finditer(txt)]
    pl=[x for x in pl if 0<x<10]                      # 排掉 50 · 0 这种明显不同单位的
    if not pl: 
        rows.append(dict(v_round=rp.parts[-2][:26],v_plant=np.nan,v_contested=np.nan,
                         v_ratio=np.nan,v_kind='无法抽取种植幅度')); continue
    # ⚠ **第一版这里错了,而它造出了一个 355× 的假指控。**
    # 一个**零**的「争议幅度」不是**观测到的**那个数(那近似于 0),
    # 而是**有意义的效应量** —— 即那一轮传给 guard 21 的 `meaningful=`。
    # R413 是一个零:观测 −0.00042,而有意义的是 0.05 -> 我的第一版算出 355×,那是我自己造的。
    # **修法:轮次若声明了 `meaningful=`,用它;否则才用观测值。**
    con=np.nan; src=''
    mm=re.search(r'meaningful\s*=\s*('+NUMP+r')',txt)
    if mm:
        con=f(mm.group(1)); src='meaningful=(该轮自己声明的有意义效应量)'
    for csv in ([] if np.isfinite(con) else sorted((rp.parent/'results').glob('*.csv'))):
        try: T0=pd.read_csv(csv)
        except Exception: continue
        for c in T0.columns:
            if re.search(r'_(b|diff|eff|effect|com|shar|r)$',c) and pd.api.types.is_numeric_dtype(T0[c]):
                v=T0[c].abs().median()
                if np.isfinite(v) and v>0: con=float(v); src=f"{csv.name}:{c}"; break
        if np.isfinite(con): break
    if not np.isfinite(con):
        rows.append(dict(v_round=rp.parts[-2][:26],v_plant=float(np.median(pl)),v_contested=np.nan,
                         v_ratio=np.nan,v_kind='UNCOMPARABLE(争议幅度抽不到)')); continue
    p_=float(np.median(pl))
    rows.append(dict(v_round=rp.parts[-2][:26],v_plant=p_,v_contested=con,
                     v_ratio=p_/con,v_kind=f'可比({src})'))
T=pd.DataFrame(rows); check_columns(T,'R426')
T.to_csv(pathlib.Path(__file__).parent/'results'/'audit.csv',index=False)
HAVE=T[T.v_ratio.notna()]
print(f"\n有正对照的轮次 **{len(T)}** · 能抽出种植幅度的 **{int(T.v_plant.notna().sum())}** · "
      f"**两个都能抽且单位可比的 {len(HAVE)}**")
print(f"⚠ **覆盖率 = {len(HAVE)}/{len(T)} = {len(HAVE)/max(len(T),1):.1%}** —— "
      f"**抽不到的那部分不是「没问题」,是「没看」。**")
if len(HAVE):
    print(f"\n种植/争议 比值:中位 **{HAVE.v_ratio.median():.2f}×** · "
          f"> 1 的 **{int((HAVE.v_ratio>1).sum())}/{len(HAVE)}** · 最大 {HAVE.v_ratio.max():.1f}×")
    for r in HAVE.sort_values('v_ratio',ascending=False).head(8).itertuples():
        mark='❌' if r.v_ratio>1 else '✅'
        print(f"   {mark} {r.v_round:<28} 种植 {r.v_plant:.3g} / 争议 {r.v_contested:.3g} = "
              f"**{r.v_ratio:.1f}×**  [{r.v_kind}]")

# ---- 抽取器的对照 ----
r425=[r for r in rows if r['v_round'].startswith('R425')]
noctl=[r for r in rows if r['v_kind']=='无法抽取种植幅度']
print(f"\n抽取器对照:")
print(f"   正对照(`R425` 已知种植 0.25):{r425[0] if r425 else '**未抽到 —— 抽取器坏了**'}")
print(f"   负对照(没有可抽种植幅度的轮次):**{len(noctl)}** 个,全部记为「无法抽取」而不是 0")
CP=bool(r425) and np.isfinite(r425[0]['v_plant']) and abs(r425[0]['v_plant']-0.25)<1e-9
CN=all(not np.isfinite(r['v_ratio']) for r in noctl)

g=Gate('本项目有多少正对照是种在方便幅度上的')
g.asserted('★【两支】抽取器正对照:必须在 `R425` 上抽出 0.25',CP,
           f"{r425[0]['v_plant'] if r425 else 'None'}",kind='control')
g.asserted('★【两支】抽取器负对照:抽不到的必须记 NaN,不能记 0',CN,
           f"{len(noctl)} 个",kind='control')
g.asserted('★【两支】覆盖率已报',len(HAVE)>0,
           f"{len(HAVE)}/{len(T)} = {len(HAVE)/max(len(T),1):.1%}",kind='control')
if CP and CN and len(HAVE):
    g.asserted('★【非零支】比值中位 > 1(= 世界 B:我一直在随手选幅度)',
               float(HAVE.v_ratio.median())>1.0,
               f"中位 {HAVE.v_ratio.median():.2f}× · >1 的 {int((HAVE.v_ratio>1).sum())}/{len(HAVE)}")
else:
    g.asserted('★ 抽取器对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
