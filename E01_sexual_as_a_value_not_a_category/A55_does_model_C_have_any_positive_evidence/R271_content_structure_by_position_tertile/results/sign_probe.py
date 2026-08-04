import os,sys,pathlib
ROOT=pathlib.Path('/home/ivan/research.psychology.eroticization-operator.operate.md.private.editable')
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
import numpy as np,pandas as pd,warnings; warnings.filterwarnings('ignore')
qm=pd.read_csv('data/derived/multiselect_questions.csv'); lg=pd.read_parquet('data/derived/endorsements_long.parquet')
keep=qm[(~qm.single_pick)&(qm.n_options>=10)&(qm.n_respondents>=1200)&(qm.mean_picks>1.5)]
signs=[]
for _,q in keep.iterrows():
    s=lg[lg.qi==q.qi]; vc=s.option.value_counts(); s=s[s.option.isin(set(vc[vc>=20].index))]
    ppl=np.array(sorted(s.person.unique())); opt=np.array(sorted(s.option.unique()))
    if len(ppl)<1200 or len(opt)<8: continue
    pi={p:i for i,p in enumerate(ppl)}; oi={o:i for i,o in enumerate(opt)}
    M=np.zeros((len(ppl),len(opt))); M[s.person.map(pi).values,s.option.map(oi).values]=1
    Z=M-M.mean(0,keepdims=True); w,v=np.linalg.eigh(np.cov(Z,rowvar=False)); pc=v[:,-1]
    signs.append((q.qi,float(pc.mean()),float(np.sign(pc).mean()),int(len(opt))))
T=pd.DataFrame(signs,columns=['qi','mean_loading','sign_balance','k'])
print(T.to_string(index=False))
print(f"\n块数 {len(T)};平均载荷为正的块 {int((T.mean_loading>0).sum())} / {len(T)}")
print(f"载荷全同号(|sign_balance|=1)的块 {int((T.sign_balance.abs()>0.999).sum())} / {len(T)}")
