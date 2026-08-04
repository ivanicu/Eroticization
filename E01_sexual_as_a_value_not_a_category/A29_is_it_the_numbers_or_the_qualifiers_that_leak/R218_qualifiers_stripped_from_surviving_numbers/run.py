import os,sys,pathlib,subprocess
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'tools'))

"""
E01 A29 R218 -- 流失的是数字,还是限定语

`#172` 的「结构上做不到」:一句被删掉的**限定语**不带数字,那道闸看不见。
而 §2 的 η 规则说 `#167`-`#171` 这一整串的病根就是限定语。**所以我预期这里会有很多。**

ESTIMAND        最近 30 个改动 README 的提交里,**某个数活了下来、而限定它的那句话没跟过来**的次数。
KILL            条件式:**先要两个对照都过**(逻辑对照 + 管道对照),否则一个从未返回非零的
                仪器给出的零是**沉默,不是无罪**(P5 ★)。过了之后:
                **存在这样的提交 -> 这条线上流失的从来不是数字。**
POWER           规则要能开火,前提是"数的新家不带限定语"这件事本身有可能发生。
                量它的**基础率**:公开面上带数量的行里,有多大比例不带限定语。
                基础率若接近 0,零就是空洞的。
IMPOSSIBLE      一句**不含数字**的限定语这条规则看不见。
"""
import pandas as pd, hashlib, re
import readme_ledger_audit as A
from lib.gates import Gate, check_coverage
OUT=pathlib.Path(__file__).parent/'results'

# ---- 对照 ① 逻辑:纯文本,不碰 git -------------------------------------------
# #173c:第一版正对照走的是 git 路,失败在**管道**上,而我差点把它读成"逻辑不成立"。
old_t='The trait split-half is +0.432 (person bootstrap, n=7,316) and beyond that scope it is UNVERIFIED.'
n_pos=len(A.qualifiers_stripped_texts(old_t,'The trait split-half is +0.432 and that settles it.'))
n_neg=len(A.qualifiers_stripped_texts(old_t,'Nothing here.'))
print(f"对照①逻辑:正 {n_pos}(摘限定语留数字)· 负 {n_neg}(数字一起删,归 numbers_that_left 管)")

# ---- 对照 ② 管道:git 那条路真的读到文本了吗 ---------------------------------
_,pl=A.qualifiers_stripped('HEAD~5')
print(f"对照②管道:{ {k:v[0] for k,v in pl.items()} } 字符读入,新版 {list(pl.values())[0][1]} 字符")

# ---- 功效:基础率 -----------------------------------------------------------
alltxt='\n'.join(pathlib.Path(f).read_text() for f in ['README.md','README_zh.md'])
numlines=[l for l in alltxt.split('\n') if l.strip() and A._MAGNUM.search(l)]
naked=[l for l in numlines if not A._QUAL.search(l)]
base=len(naked)/len(numlines)
print(f"\n功效:带数量的行 {len(numlines)},其中**不带任何限定语** {len(naked)} = **{100*base:.1f}%**")
print(f"      -> 若限定语被剥掉,数落在无限定语的行上的概率约 {100*base:.0f}%,规则**不是**被堵死的")

# ---- 回测 -------------------------------------------------------------------
commits=subprocess.run(['git','log','--format=%h','-40','--','README.md'],
                       capture_output=True,text=True,check=True).stdout.split()[:30]
rows=[]
for c in commits:
    par=subprocess.run(['git','rev-parse','--short',f'{c}^'],capture_output=True,text=True).stdout.strip()
    if not par: continue
    D,_=A.qualifiers_stripped(par,rev_to=c)
    rows.append(dict(commit=c,n_stripped=len(D),
                     tokens=', '.join(sorted(set(D.surviving))[:4]) if len(D) else ''))
T=pd.DataFrame(rows); T.to_csv(OUT/'backtest.csv',index=False)
pd.DataFrame(dict(line=[l[:160] for l in naked])).to_csv(OUT/'naked_number_lines.csv',index=False)
hit=T[T.n_stripped>0]
print(f"\n回测 {len(T)} 个提交:**{len(hit)} 个有「数活着、限定语没跟过来」**")
check_coverage(len(T),len(commits),'R218 回测样本',tol=0.05)

g=Gate('流失的是数字还是限定语')
g.asserted('对照①逻辑:摘限定语留数字必须抓到,数字一起删必须不报',n_pos>0 and n_neg==0,
           f"正 {n_pos} / 负 {n_neg}")
g.asserted('对照②管道:git 路确实读到了非空文本',all(v[0]>0 for v in pl.values()) and
           list(pl.values())[0][1]>0, f"{ {k:v[0] for k,v in pl.items()} }")
g.asserted('功效:基础率不接近 0(否则零是空洞的)',base>0.15,f"{100*base:.1f}%")
g.asserted('注册的 kill:存在「数活着、限定语没跟过来」的提交',len(hit)>0,f"{len(hit)}/{len(T)}")
print(g)
print(f"\n  => 两个对照都过、功效 {100*base:.0f}%,而回测 {len(hit)}/{len(T)}。"
      f"\n     **限定语没有被剥离。真正的数字是那个 {100*base:.1f}% —— 它们从来就没有过限定语。**")
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
