"""tools/label_sweep.py — 把两件工具接起来,扫全库

造于 `E02·A221·R593`(`#547` 的 NEXT)。行动类型:**PRODUCTION + 一次可能触发大范围重查的检验**。
`inventory.py` 判「有没有意思」→ 本脚本对判为**语义在**的来源,把 `label_audit` 的判据扫到**全部变量**。

⚠ `P2` 成本表:GSS 有 6,943 个变量、620 MB。**先只读元数据**(变量标签),
   **再只把「标签含取值词」的列读进来** —— 否则一次全读会把内存吃掉。

**预注册(先于结果写死):** 若某来源的**不相容率 > 5%**,
   则**本项目此前用过该来源的每一个变量都必须重查**,清单在本轮列出。
   ⚠ 沿用 `#546c`/`#547b` 的纪律:**不发「一定对」**;分母只计**判据适用**的列。
"""
import json, pathlib, re, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from label_audit import verdict            # 判据只有一份,不重写

# 本项目此前真正用过的变量(供 >5% 时的重查清单)
USED = {
    "gss": ["homosex", "premarsx", "xmarsex", "teensex", "pornlaw", "attend", "zodiac", "age",
            "educ", "sex", "year", "sexsex5", "cappun", "grass", "letdie1", "suicide1", "fepol",
            "natfare", "abany", "abdefect", "abnomore", "abhlth", "abpoor", "abrape", "absingle"],
    "nsfg": ["samesex", "sxok18", "sxok16", "staytog", "chunless", "chsuppor", "gayadopt", "okcohab",
             "marrfail", "chcohab", "prvntdiv", "samesexany", "oppsexany", "hadsex", "evrmarry",
             "cohever", "lifprtnr", "ager"],
}


def sweep_stata(path, cap=None):
    import pandas as pd
    it = pd.read_stata(path, iterator=True)
    vl = it.variable_labels()
    cand = [c for c, l in vl.items() if l and verdict(l, [0.0, 1.0])[0] != "N/A"]
    if cap: cand = cand[:cap]
    print(f"    变量 {len(vl)} 个 · 标签含取值词的候选 {len(cand)} 个(**只读这些列**)")
    out = []
    CH = 200
    for i in range(0, len(cand), CH):
        blk = cand[i:i + CH]
        d = pd.read_stata(path, columns=blk, convert_categoricals=False)
        for c in blk:
            vs = d[c].dropna().unique().tolist()
            v, why = verdict(vl[c], vs)
            out.append(dict(var=c, label=str(vl[c])[:60], verdict=v, why=why,
                            n=int(d[c].notna().sum()),
                            inclusion=["全样本非缺失唯一值", "标签内嵌于 .dta"]))
    return out


def sweep_dct(dct, dat):
    pat = re.compile(r'_column\((\d+)\)\s+\S+\s+(\S+)\s+%(\d+)\w?f\s+"([^"]*)"')
    lay = {}
    for line in open(dct, errors="replace"):
        m = pat.search(line)
        if m: lay[m.group(2).lower()] = (int(m.group(1)) - 1, int(m.group(3)), m.group(4))
    cand = {n: v for n, v in lay.items() if verdict(v[2], [0.0, 1.0])[0] != "N/A"}
    print(f"    变量 {len(lay)} 个 · 候选 {len(cand)} 个(**只解析这些列**)")
    buf = {n: set() for n in cand}
    for line in open(dat, errors="replace"):
        for n, (s, w, _) in cand.items():
            x = line[s:s + w].strip()
            if x not in ("", "."):
                if len(buf[n]) < 40: buf[n].add(float(x))
    out = []
    for n, (s, w, lab) in cand.items():
        v, why = verdict(lab, sorted(buf[n]))
        out.append(dict(var=n, label=lab[:60], verdict=v, why=why, n=len(buf[n]),
                        inclusion=["取值集合上限 40 个", "标签来自 .dct"]))
    return out


def main():
    res = {}
    print("=== gss(.dta 内嵌标签)===")
    res["gss"] = sweep_stata(ROOT / "data/external/gss/GSS_stata/gss7224_r3a.dta")
    print("=== nsfg 2011-2013 女性 ===")
    res["nsfg"] = sweep_dct(ROOT / "data/external/nsfg/setup/2011_2013_FemRespSetup.dct",
                            ROOT / "data/external/nsfg/2011_2013_FemRespData.dat")
    print("\n=== 逐来源的不相容率(分母 = 判据适用的列)===")
    summary, trip = {}, []
    for src, rows in res.items():
        appl = [r for r in rows if r["verdict"] != "N/A"]
        bad = [r for r in appl if r["verdict"] == "INCOMPATIBLE"]
        rate = len(bad) / len(appl) if appl else None
        summary[src] = dict(n_vars=len(rows), n_applicable=len(appl), n_incompatible=len(bad),
                            rate=rate, examples=[r["var"] for r in bad[:8]])
        print(f"  {src:6s} 判据适用 {len(appl):5d} 列 · **不相容 {len(bad):4d}** · "
              f"率 = **{rate:.4f}**" if rate is not None else f"  {src}: 无适用列")
        if rate is not None and rate > 0.05:
            hit = [v for v in USED.get(src, []) if any(r["var"] == v and r["verdict"] == "INCOMPATIBLE"
                                                       for r in rows)]
            trip.append(dict(source=src, rate=rate, used_and_incompatible=hit,
                             all_used=USED.get(src, [])))
            print(f"    ⚠ **> 5% ⇒ 预注册触发:本项目用过的 {len(USED.get(src,[]))} 个变量全部需重查**")
            print(f"    其中已判不相容的:{hit or '无'}")
    out = ROOT / "E02_condemnation_is_not_rarity/A221_semantics_present_is_not_semantics_right/R593_sweep/results"
    out.mkdir(parents=True, exist_ok=True)
    json.dump(dict(summary=summary, triggered=trip, rows=res,
                   prereg="不相容率 > 5% -> 本项目用过该来源的每个变量都要重查",
                   note="不发「一定对」;分母只计判据适用的列", unchallenged=True),
              open(out / "label_sweep.json", "w"), indent=1)
    print(f"\nwrote {out/'label_sweep.json'}")


if __name__ == "__main__":
    main()
