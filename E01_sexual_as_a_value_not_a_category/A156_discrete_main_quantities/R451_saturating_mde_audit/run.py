import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A156 R451 -- 有多少轮的主量是**离散**的,而它们的 MDE 可能结构上不可测

⚠ **CLOSURE**(§0.2)。不分离世界,**保护已有结论**。明确标注,不冒充发现。

`#406b`:符号一致性这种**离散**主量几乎没有动态范围 ——
**任何统一方向的推动都会把接近零的系数一起翻过来,MDE 扫描立刻饱和**,
于是它报出的「MDE 很小」其实是「**MDE 不可测**」。
**⇒ 本项目所有以「符号计数 / 格计数 / 越阈个数」为主量的轮次,都可能有同一个问题。**

ESTIMAND        ① 有多少轮的**主量是离散计数**(FRONTIER —— 分母我事前不知道);
                ② 其中**报过 MDE** 的,它的扫描是否**每一级都 100%**(**饱和的告密者**)。
判据(**先标支**)
                【两支】**正对照**:`R450`(已知饱和)必须被抓到;
                        **负对照**:一个已知**连续**主量的轮次(如 `R439`)必须不被标记。
                【非零支】存在饱和的轮次 -> 报名字与占比。
⚠ 覆盖率       **必须和结论一起报**(`#376b`):正则抓不到的不是「没问题」,是「没看」。
⚠ 过度指控     `#382c`/`#394c` 的教训:**正则会过度指控** -> **打印原文分类**,
                并把「**可能**(主量像是离散)」与「**确认**(扫描确实饱和)」**分开报**。
IMPOSSIBLE      ① 「主量是什么」在代码里没有显式标记 -> 只能从 `guard 26` 的调用与扫描的打印格式推断;
                ② 一个离散主量**不必然**饱和(格数多时仍有范围)-> 「可能」≠「有问题」。
"""
import numpy as np, pandas as pd, re, warnings, hashlib, collections
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns

rounds=sorted(pathlib.Path('E01_sexual_as_a_value_not_a_category').glob('A*/R*/run.py'))
DISC=re.compile(r'\(\s*np\.abs\([^)]*\)\s*>\s*\w+\s*\)\.mean\(\)|'
                r'sum\(1 for [^)]*if [^)]*>|'
                r'\.sum\(\)\s*/\s*len\(|max\(p,\s*len\(|'
                r'符号一致性|越阈率|格.*占比|计数')
MDESW=re.compile(r'检出 \*\*(\d+)/(\d+)')
rows=[]
for rp in rounds:
    t=rp.read_text()
    has_disc=bool(DISC.search(t))
    has_mde=('MDE 扫描' in t) or ('MDE =' in t)
    rows.append(dict(v_round=rp.parts[-2][:30],v_disc=has_disc,v_mde=has_mde,
                     v_path=str(rp)))
T=pd.DataFrame(rows); check_columns(T,'R451')
K=len(T); nd=int(T.v_disc.sum()); nm=int((T.v_disc&T.v_mde).sum())
print(f"⚠ **CLOSURE** —— 扫描 **{K}** 个轮次")
print(f"   主量**像是**离散计数的:**{nd}**({nd/K:.1%})")
print(f"   其中**同时报过 MDE** 的:**{nm}**")
print(f"\n⚠ **覆盖率**:这是一个**正则的**判断,而「主量是什么」在代码里**没有显式标记** ——")
print(f"   **抓不到的 {K-nd} 轮不是「没问题」,是「没看」。**")

# ---- ② 饱和的告密者:扫描是否每一级都 100% ----
def saturation(path):
    """从**结果**读(#398d 的规矩):跑一遍太贵 -> 读轮次 README 里印出的扫描表。"""
    rd=pathlib.Path(path).parent/'README.md'
    if not rd.exists(): return None
    txt=rd.read_text()
    hits=[(int(a),int(b)) for a,b in MDESW.findall(txt)]
    if len(hits)<3: return None
    return all(a==b for a,b in hits), hits
sat=[]; unsat=[]; nodata=[]
for r in T[T.v_disc&T.v_mde].itertuples():
    s=saturation(r.v_path)
    if s is None: nodata.append(r.v_round)
    elif s[0]: sat.append((r.v_round,s[1]))
    else: unsat.append(r.v_round)
print(f"\n② **饱和的告密者**(扫描的每一级都 100%):")
print(f"   **确认饱和 {len(sat)}** · 未饱和 {len(unsat)} · README 里没有可读的扫描表 {len(nodata)}")
for nm_,h in sat: print(f"      ⚠ **{nm_}** -> {h}")

# ---- 对照 ----
POS='R450_sign_across_outcomes'
NEG='R439_five_items_and_shame'
pos_flagged=bool(T[T.v_round.str.startswith(POS[:20])].v_disc.any())
neg_flagged=bool(T[T.v_round.str.startswith(NEG[:20])].v_disc.any())
print(f"\n对照:")
print(f"   **正对照** `{POS[:26]}`(`#406b` 已知饱和):{'**被抓到**' if pos_flagged else '**没抓到 —— 检测器坏了**'}")
print(f"   **负对照** `{NEG[:26]}`(主量是**连续**的偏系数):"
      f"{'⚠ 被误标' if neg_flagged else '**未被标记**'}")
T.drop(columns=['v_path']).to_csv(pathlib.Path(__file__).parent/'results'/'discrete_audit.csv',index=False)

g=Gate('有多少轮的主量是离散的,而 MDE 可能不可测')
g.asserted('★【两支】正对照:`R450`(已知饱和)必须被抓到',pos_flagged,
           f"{POS[:26]}",kind='control')
g.asserted('★【两支】负对照:连续主量的轮次不得被标记',not neg_flagged,
           f"{NEG[:26]}",kind='control')
g.asserted('★【两支】覆盖率已报',nd>0,f"{nd}/{K} = {nd/K:.1%}",kind='control')
if pos_flagged and not neg_flagged:
    g.asserted('★【非零支】存在**确认饱和**的轮次',len(sat)>0,
               f"确认 {len(sat)} · 可能(离散+有 MDE){nm} · 全部离散 {nd}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
