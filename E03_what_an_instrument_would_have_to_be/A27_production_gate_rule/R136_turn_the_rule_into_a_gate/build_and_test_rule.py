"""E03·A27·R136 —— 把「每条结论必须带自己的零」变成闸上的一条规则

**类型:PRODUCTION(诚实标注)。不推进对象,不检验任何声明。**
`#692` 与 `#693` 连续两轮的失败都指向同一件事:**缺的不是分析,是这条规矩没有被机械执行。**
`#693` 量出:账本后段 93 条里只有 **14–22%** 报了自己的零 ⇒ **再清点第三次徒劳,数据不存在。**

## 仪器(硬规则②)
**唯一的仪器是账本文本自身。没有第二具仪器** —— 本轮不读任何数据集,「跨仪器」不适用。

## ③ 合入判据(`#693` 在跑之前写死,本脚本先自证再改闸)
**新规则必须在两个种入的伪条目上分别开火与放行,否则不合入。**
## ⑤ 最强混淆(`#693` 预注册)
「零的表述」至少四种写法(零 95% 分位 / 打乱…的零 / 置换零 / 自助区间)
⇒ **规则必须四种都认,否则会把合格条目误判** —— **先用 `#693` 那 20 条全量回测,误报一条都不许有。**
"""
import os, sys, pathlib, json, re, warnings
ROOT = pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")
OUT=pathlib.Path(__file__).parent/"results"; OUT.mkdir(exist_ok=True)

EFFECT = re.compile(r'\*\*[-+]\d*\.\d{3,4}\*\*')
NULLS  = re.compile(r'零\s*(?:的)?\s*95%|打乱[^。\n]{0,20}零|置换零|'
                    r'自助[^。\n]{0,10}区间|95%\s*(?:自助)?区间|零的种类|null_kind|'
                    r'零\s*95%\s*分位|经验\s*p|p\s*=\s*\*{0,2}\d')
OPTOUT = re.compile(r'本轮不报效应|没有效应|不检验任何声明|无零可报|结构性拿不到零')

def effect_without_null(text, cutoff=693):
    """返回 (blocking, grandfathered)。
    PROPERTY   一条带了效应的结论,却没报它自己的零
    PROXY      正文含 `**±0.xxxx**` 而不含任何一种零的表述,且无豁免语
    IMPLICATION 只有一个方向可靠:**含效应且不含任何零表述 -> 它确实没报零**(可靠)。
               反过来不成立:含零表述**不**证明那个零配的是这个效应(`#693`:猜配对不是测量)。
    SAFE SIDE  只报「没报零」;从不报「这一条的零配对正确」。
    """
    marks=[(int(m.group(1)),m.start()) for m in re.finditer(r'^## Entr(?:y|ies) (\d+)',text,re.M)]
    blocking,old=[],[]
    for i,(n,s0) in enumerate(marks):
        body=text[s0:(marks[i+1][1] if i+1<len(marks) else len(text))]
        if not EFFECT.search(body): continue
        if NULLS.search(body) or OPTOUT.search(body): continue
        (blocking if n>=cutoff else old).append(n)
    return blocking,old

led=open("RETRACTIONS.md").read()
print("=== ③ 种入测试:两条伪条目,必须一拦一放 ===")
BAD  = "\n## Entry 9001 · `E9·A9·R9` — 伪条目:报了效应没报零\n效应是 **+0.1234**,就这样。\n"
GOOD = "\n## Entry 9002 · `E9·A9·R9` — 伪条目:报了效应也报了零\n效应是 **+0.1234**,打乱 educ 的零 95% 分位 0.03。\n"
b1,_=effect_without_null(led+BAD);  hit_bad  = 9001 in b1
b2,_=effect_without_null(led+GOOD); hit_good = 9002 in b2
print(f"  伪条目 A(有效应无零)必须被拦:{'✅ 拦住' if hit_bad else '⛔ 没拦住'}")
print(f"  伪条目 B(有效应有零)必须放行:{'✅ 放行' if not hit_good else '⛔ 误拦'}")
print("\n=== ⑤ 全量回测:#693 认定报了零的那 20 条,一条都不许被误报 ===")
P4={"strict":r'零\s*(?:的)?\s*95%',"shuffle":r'打乱[^。\n]{0,20}零',"perm":r'置换零',
    "boot":r'自助[^。\n]{0,10}区间|95%\s*(?:自助)?区间'}
marks=[(int(m.group(1)),m.start()) for m in re.finditer(r'^## Entr(?:y|ies) (\d+)',led,re.M)]
have20=[]
for i,(n,s0) in enumerate(marks):
    if n<600: continue
    body=led[s0:(marks[i+1][1] if i+1<len(marks) else len(led))]
    if any(re.search(v,body) for v in P4.values()): have20.append(n)
blk,gr=effect_without_null(led)
false_pos=[n for n in have20 if n in blk or n in gr]
print(f"  #693 认定报了零的条目:**{len(have20)} 条**")
print(f"  其中被新规则标记的(= 误报):**{len(false_pos)} 条** {false_pos if false_pos else '—'}")
print(f"  ⇒ {'✅ 零误报' if not false_pos else '⛔ 有误报,规则不合入'}")
print(f"\n=== ④ 基线处置 ===")
print(f"  阻断(#693 起):**{len(blk)} 条** {blk}")
print(f"  grandfathered(#693 之前,只点名):**{len(gr)} 条**")
ok = hit_bad and (not hit_good) and (not false_pos)
print(f"\n**⇒ {'三项判据全过 —— 合入' if ok else '有判据未过 —— 不合入,如实记'}**")
json.dump(dict(plant_bad_blocked=bool(hit_bad),plant_good_passed=bool(not hit_good),
               n_reported_null=len(have20),false_positives=false_pos,
               blocking=blk,grandfathered_n=len(gr),merge=bool(ok)),
          open(OUT/"rule_test.json","w"),indent=1,ensure_ascii=False)
