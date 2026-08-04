#!/usr/bin/env python3
"""
tools/round_status.py -- 给定一个轮次,它支撑的声明**当前**的账本状态是什么?

由 #154c 触发。本会话三次犯同一个错:**去看轮次,而不是去读账本**。
  #143  重跑 `A02/R10` 去"重新定价"一条 `Entry 24` 在 118 条之前就已作废的声明
  #148e 两轮找一个 `#130b` 已经定位过的机制
  #154  重抽 `A11/R14` 的零 —— 而 `Entry 101` 写着 `R15` 已经取代了它

`readme_ledger_audit.py` 做的是 README → 账本。这个做**反方向**:轮次 → 账本,
并把**最后一条提到它的条目**顶出来。有了它,"我要重跑 R14"时一眼能看到
"Entry 101 说 R15 取代了它"。

⚠ P6 代理账:
  PROPERTY   这个轮次的结论,当前是否仍然是那条声明的现行版本
  PROXY      账本里最后一条**提到该轮次**的条目,及其是否含取代/修正类词汇
  IMPLICATION 只有一个方向可靠:**最后一条提到它的条目含取代语言 -> 必须去读那条**(可靠)。
             反过来"没有取代语言 -> 它是现行的"**不可靠** —— 取代可能没点名这一轮。
  SAFE SIDE  输出是**必读清单**,不是判决。
"""
import re,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
SUP=('取代','supersede','superseded','修法','已作废','撤回','降级','WITHDRAWN','retract',
     'RETRACTED','的修法','不再','改用','换成')
RID=re.compile(r'E01[·.]?(A\d{2})[·.]?(R\d{2})')
HEAD=re.compile(r'^## Entry (\d+),\s*added by\s*(.+)$')

def index(ledger='RETRACTIONS.md'):
    lines=(ROOT/ledger).read_text().splitlines()
    cur=None; added={}; cited={}
    for i,l in enumerate(lines,1):
        h=HEAD.match(l)
        if h:
            cur=int(h.group(1))
            for a,r in RID.findall(h.group(2)):
                added.setdefault(f"{a}/{r}",[]).append(cur)
            # 头里的 `+`R15`` 形式
            for r in re.findall(r'\+\s*`?(R\d{2})`?',h.group(2)):
                for a,_ in RID.findall(h.group(2)):
                    added.setdefault(f"{a}/{r}",[]).append(cur)
            continue
        if cur is None: continue
        for a,r in RID.findall(l):
            cited.setdefault(f"{a}/{r}",[]).append((cur,i,l.strip()))
        for r in re.findall(r'`(R\d{2})[_`]',l):
            pass
    return added,cited,lines

def bodies(lines):
    """每条 Entry 的正文行范围。"""
    out={}; cur=None; start=None
    for i,l in enumerate(lines):
        h=HEAD.match(l)
        if h:
            if cur is not None: out[cur]=(start,i)
            cur=int(h.group(1)); start=i
    if cur is not None: out[cur]=(start,len(lines))
    return out

def sibling_supersession(rd,ent,lines,B):
    """⚠ #155:救不了我的那一条 —— `Entry 101` 正文写的是**裸 `R15`**,不是 `E01·A11·R15`。
    所以要扫这一轮**自己那条条目的正文**,找里面提到的**兄弟轮次** + 取代类词汇。
    这正是当时会救我的信号:`Entry 101` 正文里同时有 `R15` 和「fix」「failed」。"""
    import re as _re
    if ent not in B: return None
    a,b=B[ent]; me=rd.split('/')[1]
    hits=[]
    for l in lines[a:b]:
        sibs=set(_re.findall(r'`?(R\d{2})`?',l))-{me}
        if sibs and any(w in l for w in SUP+('fix','failed','replacement','instead')):
            hits.append((sorted(sibs),l.strip()[:150]))
    return hits or None

def status(paths=None):
    added,cited,lines=index()
    rounds=sorted(set(added)|set(cited))
    out=[]
    for rd in rounds:
        ev=sorted(set([e for e in added.get(rd,[])]+[c for c,_,_ in cited.get(rd,[])]))
        if not ev: continue
        last=max(ev)
        ctx=[t for c,_,t in cited.get(rd,[]) if c==last]
        sup=any(any(w in t for w in SUP) for t in ctx)
        B=bodies(lines); sib=None
        for e in added.get(rd,[]):
            s=sibling_supersession(rd,e,lines,B)
            if s: sib=(e,s); break
        out.append(dict(round=rd,added=added.get(rd,[]),last=last,n_entries=len(ev),
                        superseded_hint=sup or bool(sib),
                        ctx=(ctx[0][:150] if ctx else (sib[1][0][1] if sib else '')),
                        sib=(f"Entry {sib[0]} 正文提到兄弟轮次 {sib[1][0][0]}" if sib else '')))
    return out

if __name__=='__main__':
    rows=status()
    want=sys.argv[1:] if len(sys.argv)>1 else None
    print(f"账本里出现过的轮次:{len(rows)} 个\n")
    print(f"  {'轮次':<10}{'产出条目':<14}{'最后提到':>8}{'条目数':>7}  取代提示")
    for r in sorted(rows,key=lambda x:-x['last']):
        if want and not any(w in r['round'] for w in want): continue
        mark='⚠ 必读' if r['superseded_hint'] else ''
        print(f"  {r['round']:<10}{str(r['added'])[:13]:<14}{r['last']:>8}{r['n_entries']:>7}  {mark}")
        if r['superseded_hint']:
            if r['sib']: print(f"       └─ {r['sib']}")
            print(f"       └─ {r['ctx']}")
    print("\n⚠ SAFE SIDE(#P6):只在**命中**方向可读。没有取代提示**不等于**这一轮是现行的。")
