import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A157 R452 -- `S` 反号的那 8 个结局有共同点吗:用**卡片自己的**分类,不用我的读法

`#406a`:`S` 在 29 个结局里 21 正 **8 负**,而最负的三个是 `highenergy` · **宜人性** · **尽责性**。
**读起来像「关于自己性的陈述 vs 一般人格」—— 但那是我读出来的。**

⚠ 而这份数据里有一个**不由我给**的分类:**数据集卡片自己的一节**
`COMPUTED COLUMNS (not directly from GT questions)` ——
它逐条点名了 `opennessvariable` · `consciensiousnessvariable` · `extroversionvariable` ·
`neuroticismvariable` · `agreeablenessvariable` · `powerlessnessvariable` ·
`totalfetishcategory` · `bondageaverage` · `Total*`。
**⇒ 「合成量 vs 单题」这个分组由卡片给出,不由我给出。**(`#374a`:卡片的删除清单是可信的那一部分。)

两个活着的世界:
**A 分类解释它** -> 按卡片分组后,**组内符号一致性显著高于随机等大小分组**;
**B 不解释** -> 那 8 个反号仍无解,而**「我读出的那个共同点」不能当结论**。

ESTIMAND        29 个结局按**卡片的 `COMPUTED` 清单**分两组;
                主量 = **两组各自的符号一致性的加权平均**。
判据(**先标支**,`#379c`)
                【两支】负对照用**越阈率**;
                        **guard 26 显式传 `main_quantity='discrete_count'` 与 `sweep_detection`**
                        (`#407c` 刚做的接口,**第一次真用**)。
                【非零支】越过 offset 零 -> 世界 A;
                【零支】落在零里 -> 世界 B,而**我那个读法不能上页面**。
⚠ 零的种类     `offset_control`:**分组后一致性的零绝不是零** ——
                **任意**分组都会把组内一致性抬高(把同号的分到一起纯属偶然也会发生)。
                零 = **随机等大小分组**(组大小照旧)的一致性分布。
IMPOSSIBLE      ① 「合成量」组很小(≈6)-> 一致性在小组上天然更高,**而随机零用同样的组大小,已经吸收了这一点**;
                ② 卡片描述的不是这个文件(`#373`)-> 但**删除/合成清单那一部分是可信的**(`#374`);
                ③ 「分类解释符号」不等于「符号有心理学含义」。
"""
import numpy as np, pandas as pd, warnings, hashlib
warnings.filterwarnings('ignore')
from lib.gates import Gate, check_columns
from lib.nulls import perm_in, controls as null_controls
nc=null_controls(); nc.pop('_detail'); assert all(nc.values())
T=pd.read_csv('E01_sexual_as_a_value_not_a_category/A155_does_S_point_one_way/'
              'R450_sign_across_outcomes/results/S_signs.csv')     # ⚠ 从 results/ 读(#398d)
CARD=pathlib.Path('data/card/column_notes.txt').read_text()
sec=CARD.split('COMPUTED COLUMNS')[1].split('DROPPED GT QUESTIONS')[0]
TOK=[w for w in ('opennessvariable','consciensiousnessvariable','extroversionvariable',
                 'neuroticismvariable','agreeablenessvariable','powerlessnessvariable',
                 'totalfetishcategory','bondageaverage') if w in sec]
print(f"⚠ 分组由**卡片自己的 `COMPUTED COLUMNS` 一节**给出,不由我给出。")
print(f"   卡片在那一节点名的、且本轮用得上的标识:**{len(TOK)}** 个 -> {TOK}")
assert len(TOK)>=6, "卡片那一节没读到预期的标识 —— 停手"
ALIAS={'openness':'opennessvariable','conscientiousness':'consciensiousnessvariable',
       'extroversion':'extroversionvariable','neuroticism':'neuroticismvariable',
       'agreeableness':'agreeablenessvariable','powerlessness':'powerlessnessvariable'}
def is_computed(nm):
    s=str(nm); low=s.lower()
    if s in ALIAS: return True
    if low.startswith('total'): return True
    return any(t.lower() in low for t in TOK)
T['comp']=T.v_out.map(is_computed)
K=len(T); nc_=int(T.comp.sum())
print(f"\n29 个结局里:**合成量 {nc_}** · **单题 {K-nc_}**")
print(f"   合成量组:{T[T.comp].v_out.tolist()}")
def cons(sub):
    if len(sub)==0: return np.nan
    p=int((sub>0).sum()); return max(p,len(sub)-p)/len(sub)
def weighted(mask):
    a=cons(T.v_b[mask]); b=cons(T.v_b[~mask])
    na,nb=int(mask.sum()),int((~mask).sum())
    return (a*na+b*nb)/(na+nb)
OBS=weighted(T.comp.values)
print(f"\n分组后的加权一致性 = **{OBS:.4f}**")
print(f"   合成量组 **{cons(T.v_b[T.comp]):.4f}**({nc_} 个,为正 {int((T.v_b[T.comp]>0).sum())})· "
      f"单题组 **{cons(T.v_b[~T.comp]):.4f}**({K-nc_} 个,为正 {int((T.v_b[~T.comp]>0).sum())})")
NP_=2000
rg=np.random.default_rng(11); nul=[]
for _ in range(NP_):
    m=np.zeros(K,bool); m[rg.choice(K,nc_,replace=False)]=True
    nul.append(weighted(m))
nul=np.array(nul); HI=float(np.percentile(nul,95))
print(f"\n⚠ offset 零(**随机等大小分组** {NP_} 次;"
      f"**任意分组都会把组内一致性抬高 -> 这个零绝不是零**):")
print(f"   **{nul.mean():.4f} ± {nul.std():.4f}** · 95 分位 **{HI:.4f}**")
print(f"   实测 **{OBS:.4f}** -> **{(OBS-nul.mean())/max(nul.std(),1e-12):+.2f} sd** · "
      f"{'**越阈 -> 世界 A(卡片的分类解释它)**' if OBS>HI else '**落在零里 -> 世界 B**'}")
negs=[]
rg2=np.random.default_rng(77)
for _ in range(400):
    m=np.zeros(K,bool); m[rg2.choice(K,nc_,replace=False)]=True
    negs.append(weighted(m))
negs=np.array(negs); rate=float((negs>HI).mean())
print(f"\n负对照(**越阈率**,随机分组 400 次):**{100*rate:.1f}%**(合格 1–12%)")
print(f"\nguard 26 = **MDE 扫描**(⚠ 主量是**离散计数** -> 本轮**显式声明**,`#407c` 的接口):")
det=[]
for gg in (0.2,0.4,0.6,0.8):
    hit=0
    for s_ in range(30):
        rg3=np.random.default_rng(60+int(gg*10)*103+s_)
        b2=T.v_b.values.copy()
        idx=np.flatnonzero(T.comp.values)
        flip=rg3.random(len(idx))<gg
        b2[idx[flip]]=-np.abs(b2[idx[flip]])           # 把合成量组按比例推成同号(负)
        Tt=T.copy(); Tt['v_b']=b2
        a=cons(Tt.v_b[Tt.comp]); b=cons(Tt.v_b[~Tt.comp])
        w=(a*nc_+b*(K-nc_))/K
        if w>HI: hit+=1
    det.append(hit/30); print(f"   合成量组按 **{gg:.0%}** 推成同号 -> 检出 **{hit}/30 = {hit/0.3:>5.1f}%**")
MDE=next((g for g,dv in zip((0.2,0.4,0.6,0.8),det) if dv>=0.8),1.0)
print(f"   **MDE = {MDE:.0%}** · 有意义的分组效应 = **50%**")
T[['v_out','v_b','comp']].to_csv(pathlib.Path(__file__).parent/'results'/'by_card_class.csv',index=False)
NONNULL=OBS>HI
g=Gate('S 反号的那 8 个结局有共同点吗')
g.asserted('★【两支】负对照:**越阈率** ≈5%',0.01<=rate<=0.12,f"{100*rate:.1f}%",kind='control')
g.positive_control_at_the_contested_magnitude(
    '★【两支】guard 26(**显式声明离散主量 + 扫描**,`#407c` 的接口)',MDE,0.50,True,
    what='MDE 扫描',branch='non_null' if NONNULL else 'null',
    main_quantity='discrete_count',sweep_detection=det)
g.asserted('★【两支】offset 零非退化(任意分组都会抬高一致性)',nul.std()>0,
           f"{nul.mean():.4f} ± {nul.std():.4f}",kind='control')
if 0.01<=rate<=0.12:
    g.asserted('★【非零支】卡片的分类解释符号 -> 世界 A',NONNULL,
               f"{OBS:.4f} vs 上侧 {HI:.4f}")
else:
    g.asserted('★ 对照未过 -> 不判',False,'UNVERIFIED')
print(g)
print(f"\nSHA {hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()[:12]}")
