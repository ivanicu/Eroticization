import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

import pandas as pd
qm = pd.read_csv('data/derived/multiselect_questions.csv')
lg = pd.read_parquet('data/derived/endorsements_long.parquet')
keep = qm[(~qm.single_pick)&(qm.n_options>=8)&(qm.n_respondents>=800)]
opts = (lg[lg.qi.isin(keep.qi)].groupby(['qi','option']).size().rename('n').reset_index())
opts = opts.merge(keep[['qi','col','n_respondents']], on='qi')
opts['rate'] = (opts.n/opts.n_respondents).round(3)
opts['item_id'] = opts.qi.astype(str)+'::'+opts.option
opts[['item_id','qi','col','option','n','n_respondents','rate']].to_csv('data/derived/options.csv', index=False)
print("options written:", len(opts))
for qi in [67, 29, 36, 62]:
    b = opts[opts.qi==qi].sort_values('rate', ascending=False)
    print(f"\n===== qi={qi}  {b.col.iloc[0][:70]}  (n={b.n_respondents.iloc[0]}) =====")
    for _,r in b.iterrows(): print(f"   {r.rate:5.3f}  {r.option[:105]}")
