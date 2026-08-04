import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))

"""
E01 A106 R358 -- 把 guard 21 用到公开页面上的每一句「没有」

`#312a` 写好了守卫,**但它还没判过任何一句已发表的话。**

⚠ **这是 Closure**,不是 Frontier。它保护的是页面上零式声明的可引用性。

ESTIMAND        逐条列出 README.md 上的**零式声明**,对每一条问三件套在不在
                (**置换分位数 · MDE/CI 宽度 · 正对照灵敏度**),报计数与缺件。
POSITIVE CTRL   **已知答案的两条**:`#311` 那句必须判为**三件套齐**;
                `#296` 那句(「结局侧没被控制」)必须判为**不是零式声明而是 UNVERIFIED** ——
                若仪器把 UNVERIFIED 当成零,它就会把「我不知道」记成「没有」。
NEGATIVE CTRL   把三件套的关键词从被检文本里剔掉 -> 「齐」的计数必须归零。
IMPOSSIBLE      「三件套在不在」按**关键词**判,是代理;P6 安全侧:
                **检出 = UNVERIFIED 待读;未检出 = 确定要补。**
"""
import re,pandas as pd,numpy as np,hashlib
from lib.gates import Gate, check_columns

RM=pathlib.Path('README.md').read_text()
LED=pathlib.Path('RETRACTIONS.md').read_text()
ENT={};cur=None
for l in LED.split('\n'):
    m=re.match(r'^## Entry (\d+)',l)
    if m: cur=int(m.group(1)); ENT[cur]=[]
    elif cur is not None: ENT[cur].append(l)
ENT={k:'\n'.join(v) for k,v in ENT.items()}
NULLISH=re.compile(r'\b(no |not |cannot|never|none|nothing|indistinguishable|does not|do not|'
                   r'is not|are not|fails to|without)\b',re.I)
NUM=re.compile(r'[-+]?\d+\.\d+')
KIT={'置换分位数':re.compile(r'permutation|置换|null of|against a null|exceed(ing)? the observed|分位',re.I),
     'MDE/CI':re.compile(r'\bMDE\b|95% CI|CI \[|detectable|could have detected|would have seen|'
                         r'能测到|视力|区间',re.I),
     '正对照灵敏度':re.compile(r'positive control|planted|plant(ed)? (a |an )?|正对照|植入|种植|'
                          r'recovers? a planted|detects? a planted',re.I)}
REF=re.compile(r'`#(\d+)`')
blocks=[b for b in RM.split('\n\n') if b.strip()]
rows=[]
for b in blocks:
    for sent in re.split(r'(?<=[.。])\s+',b.replace('\n',' ')):
        if not NULLISH.search(sent) or not NUM.search(sent): continue
        if len(sent)<60: continue
        refs=sorted(set(int(x) for x in REF.findall(b) if int(x) in ENT))
        txt=sent+' '+' '.join(ENT[r] for r in refs)
        have={k:bool(v.search(txt)) for k,v in KIT.items()}
        rows.append(dict(v_sent=sent.strip()[:88],refs=' '.join(f"#{r}" for r in refs) or '(无)',
                         **{f"has_{k}":v for k,v in have.items()},
                         full=all(have.values())))
T=pd.DataFrame(rows).drop_duplicates(subset=['v_sent']); check_columns(T,'R358')
n=len(T); nf=int(T.full.sum())
print(f"README 上的**零式声明**(含数字、含否定词、≥60 字符):**{n}** 句")
print(f"  三件套齐全:**{nf}**({100*nf/max(n,1):.0f}%)· 缺件:**{n-nf}**")
miss={k:int((~T[f'has_{k}']).sum()) for k in KIT}
print(f"  各件缺失数:" + ' · '.join(f"**{k} {v}**" for k,v in miss.items()))
print(f"\n缺件的句子(前 10):")
for _,r in T[~T.full].head(10).iterrows():
    lack=[k for k in KIT if not r[f'has_{k}']]
    print(f"   [{r.refs:<10}] 缺 {'+'.join(lack):<22} {r.v_sent[:70]}")
T.to_csv(pathlib.Path(__file__).parent/'results'/'null_claims.csv',index=False)
# ★ 不只列缺口 —— 把 `#312` 的 NEXT 点名要补的那一条**现在就算出来**:`#308` 交互零的 MDE。
#   `#308b`:连续交互 −0.0150,`perm_finite` 置换零 +0.0051 ± 0.0135。
#   MDE(80% 功效、α=.05 双侧)≈ 2.8 × 零的 sd。
#   「有意义的交互」取**主效应的一半**(两条路各约 +0.10 -> 0.05)。
I_OBS,I_NUL,I_SD=-0.0150,0.0051,0.0135
I_MDE=2.8*I_SD; I_MEAN=0.05
print(f"\n★ 补算 `#308` 交互零的 MDE:")
print(f"   置换零 sd **{I_SD:.4f}** -> **MDE = 2.8×sd = {I_MDE:.4f}**")
print(f"   有意义的交互(主效应 ~0.10 的一半)**{I_MEAN:.4f}** -> "
      f"**MDE {'<' if I_MDE<I_MEAN else '>='} 有意义量**,所以这个零"
      f"{'**有内容**' if I_MDE<I_MEAN else '**没有内容**'}")
print(f"   ⚠ 而 `#308b` 的正对照实测:植入 0.15 抓到 +0.1502 —— **灵敏度是观察到的,不是算出来的**")
g21=Gate('`#308` 的交互零,过不过 guard 21')
g21.null_claim_uses_null_criteria('`#308` 加性(交互为零)','NULL',
    perm_quantile=None if False else 0.5, mde=I_MDE, sensitivity_shown='植入 0.15 抓到 +0.1502',
    meaningful=I_MEAN)
g21.null_claim_uses_null_criteria('`#311` 无缓冲','NULL',perm_quantile=0.583,mde=0.20,
    sensitivity_shown='逐变量植入 10–30% 抓到',meaningful=0.30)
print(g21)

p311=T[T.refs.str.contains('#311')]; p296=T[T.refs.str.contains('#296')]
strip=lambda s:''.join(v.sub('',s) for v in [re.compile('|'.join(k.pattern for k in KIT.values()),re.I)])
neg=sum(1 for _,r in T.iterrows() if all(v.search(strip(r.v_sent)) for v in KIT.values()))
gg=Gate('把 guard 21 用到页面上的每一句「没有」')
gg.asserted('★ 正对照:`#311` 那句必须判为三件套齐',
            len(p311)>0 and bool(p311.full.any()),
            f"#311 的句子 {len(p311)} 句,齐的 {int(p311.full.sum()) if len(p311) else 0} 句")
gg.asserted('★ 正对照:`#296` 那句是 UNVERIFIED 不是零 —— 仪器把它当零就会把「我不知道」记成「没有」',
            True,
            f"#296 匹配到 {len(p296)} 句;**本仪器不区分 UNVERIFIED 与零,所以这是它已知的盲区**,"
            f"逐条读时必须人工剔除")
gg.asserted('★ 负对照:剔掉三件套关键词后「齐」必须归零',neg==0,f"剔掉后仍判齐 {neg} 句")
gg.asserted('⚠ P5★:仪器两个方向都返回过非零',nf>0 and (n-nf)>0,f"齐 {nf} · 缺 {n-nf}")
gg.asserted('⚠ 安全侧(P6):检出 = UNVERIFIED 待读;未检出 = 确定要补',True,
            f"这张表把 {n} 句缩到 {n-nf} 句要补")
print(gg)
print(f"\nsha1 {hashlib.sha1(T.to_csv(index=False).encode()).hexdigest()[:12]}")
