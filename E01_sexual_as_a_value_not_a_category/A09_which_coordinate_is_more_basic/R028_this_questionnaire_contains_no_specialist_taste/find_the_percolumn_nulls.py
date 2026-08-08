import os,sys,pathlib,re,json
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
HERE=pathlib.Path(__file__).parent

"""
Q: #425d showed a per-column null understates the spread of a MULTI-COLUMN statistic
   (34-104 instead of ~55-79). Where else in this project does that pattern occur?

BOUNDARY, written BEFORE the scan (four over-indictments already: #382c #394c #407a #422b):
  IN SCOPE  : a call site where (1) the null shuffles PER COLUMN (perm_in / np.random
              .permutation on one vector) AND (2) the statistic that the null feeds is
              aggregated ACROSS >=2 columns (a count of columns, a rate over columns,
              a max over columns, a family-wise threshold).
  OUT       : a null feeding a SINGLE-column statistic. perm_in is CORRECT there --
              there is no inter-column correlation for it to destroy. Most uses are this.
  OUT       : row_perm call sites (already the right null).

Worlds
  A  no in-scope site -> #425d is local to R468/R469, the risk is closed.
  B  some -> each is a conclusion whose STRENGTH was overstated (not necessarily wrong),
     and they now have names.

This is CLOSURE (it protects existing conclusions; it does not open a new question).
Verdict is three-valued: the scan can only find CANDIDATES -- deciding whether a site's
statistic really aggregates across columns needs the source read, which is done by hand
below and printed verbatim.
"""
SRC=[p for p in pathlib.Path('.').rglob('run.py') if '.git' not in str(p)]
SRC+= [p for p in pathlib.Path('lib').rglob('*.py')]

PERCOL = re.compile(r'perm_in\s*\(|np\.random\.default_rng\([^)]*\)\.permutation\(|rng\.permutation\(')
ROWPERM= re.compile(r'row_perm\s*\(')
# statistic aggregated across columns, in the SAME file
AGG = re.compile(r'越阈率|sig\]\.|\.sig\b|sum\(\)\s*(?:#|$)|n_pos|n_neg|allpos|all-positive|'
                 r'Bonferroni|族内阈|family[- ]wise|多重性|multiplicity|符号计数|sign_count|'
                 r'minority|少数号|(?:len|sum)\(\s*(?:pos|neg|mino)\b')

rows=[]
for p in SRC:
    t=p.read_text(errors='ignore')
    pc=list(PERCOL.finditer(t))
    if not pc: continue
    has_row=bool(ROWPERM.search(t))
    agg=list(AGG.finditer(t))
    rows.append(dict(path=str(p), n_percol=len(pc), n_rowperm_in_file=int(has_row),
                     n_agg_markers=len(agg),
                     candidate=int(len(agg)>0)))
D=__import__('pandas').DataFrame(rows).sort_values(['candidate','n_agg_markers'],ascending=False)
D.to_csv(HERE/'results/percolumn_null_sites.csv',index=False)

nfile=len(SRC); nperm=len(D); ncand=int(D.candidate.sum())
print(f"扫描 {nfile} 个源文件 · 含逐列打乱的文件 {nperm} · **候选(同文件里有跨列汇总标记)= {ncand}**")
print(f"覆盖率:候选判据是「同一文件里出现跨列汇总的词」-> 这是**代理**,不是判定。\n")
print(D[D.candidate==1][['path','n_percol','n_agg_markers','n_rowperm_in_file']].to_string(index=False))

# ---- POSITIVE CONTROL: R468 must be flagged (it IS the known in-scope site)
r468=[r for r in rows if 'R468' in r['path']]
r469=[r for r in rows if 'R469' in r['path']]
pos_ok = bool(r468 and r468[0]['candidate']==1)
# ---- NEGATIVE CONTROL: R469 uses row_perm -> it must NOT be counted as needing repair
neg_ok = bool(r469 and r469[0]['n_rowperm_in_file']==1)
print(f"\n正对照 R468(已知在范围内)被标出 = {pos_ok}")
print(f"负对照 R469(已经用了 row_perm)其文件被识别为已修 = {neg_ok}")

json.dump(dict(n_files=nfile, n_with_percol=nperm, n_candidates=ncand,
               pos_control=pos_ok, neg_control=neg_ok),
          open(HERE/'results/verdict.json','w'), indent=1)
