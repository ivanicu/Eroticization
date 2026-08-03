import os,sys,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[3]
os.chdir(ROOT); sys.path.insert(0,str(ROOT))
from lib.rounds import round_path

"""
ITER 3 separator. Is a coordinate ACQUIRED AS A UNIT, or is it a descriptive summary of
separately-acquired interests?
  World G (unit acquisition): categories a person likes for the same underlying reason were
     acquired at the same time -> within-person onset similarity tracks preference similarity.
  World I (independent acquisition): onset similarity is unrelated to preference structure;
     it is only recall style and base rates.
RSA between two 31x31 matrices built on the SAME pairwise-complete people, both double-centered
(person recall style and category mean removed). Current age partialled out (a 19-year-old
cannot report onset 26+). Null permutes category labels.
"""
import pandas as pd, numpy as np, re, warnings
from numpy.linalg import lstsq
warnings.filterwarnings('ignore'); rng=np.random.default_rng(6180)
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
ons=[c for c in inv[inv['kind']=='AGE_ONSET']['col'] if df[c].map(BIN).notna().sum()>300]
rate=inv[inv['kind']=='RATING_0_5']['col'].tolist()

# map onset category -> its arousal rating column, by keyword
KEY={'bondage':'bondage','humiliation':'humiliation','nonconsent':'nonconsent','sadomasochism':'sadomasochism',
 'sensory':'sensory','transformation':'transform','specific roles':'roles','mental alteration':'mentalalteration',
 'pregnancy':'pregnancy','genderplay':'genderplay','exhibitionism':'exhibitionself','multiple partners':'multiplepartners',
 'incest':'incest','bestiality':'bestiality','abnormal body':'abnormalbody','bodily-secretions':'secretions',
 'mythical':'mythical','creepy':'creepy','brutality':'brutality','vore':'vore','clothing':'clothing',
 'body parts':'appearance','gentleness':'gentleness','power dynamics':'powerdynamic','dirtiness':'dirty',
 'eagerness':'eagerness','objects':'objects','toys':'toys'}
pairs=[]
for c in ons:
    lc=c.lower()
    for k,v in KEY.items():
        if k in lc and v in df.columns:
            pairs.append((c,v)); break
pairs=list(dict.fromkeys(pairs))
print(f"matched onset<->rating categories: {len(pairs)}")
names=[re.sub(r'.*interest in |.*sexual interest in |\?.*','',c)[:26] for c,_ in pairs]

O=pd.DataFrame({v:df[c].map(BIN) for c,v in pairs})
P=pd.DataFrame({v:pd.to_numeric(df[v],errors='coerce') for _,v in pairs})
AGEMAP={'14-17':0,'18-20':1,'21-24':2,'25-28':3,'29-32':4}; age=df['age'].map(AGEMAP)

def resid_mat(M, extra=None):
    """double-centre (person recall style + category mean) and partial out current age"""
    Z=M.copy()
    Z=Z.sub(Z.mean(axis=1),axis=0)              # person mean
    Z=Z.sub(Z.mean(axis=0),axis=1)              # category mean
    if extra is not None:
        for c in Z.columns:
            m=Z[c].notna()&extra.notna()
            if m.sum()>200:
                X=np.c_[np.ones(m.sum()),extra[m].values]
                b,*_=lstsq(X,Z.loc[m,c].values,rcond=None)
                Z.loc[m,c]=Z.loc[m,c].values-X@b
    return Z
Ores=resid_mat(O, age); Pres=resid_mat(P)

k=len(pairs)
CO=np.full((k,k),np.nan); CP=np.full((k,k),np.nan); NP=np.zeros((k,k))
for i in range(k):
    for j in range(i+1,k):
        m=Ores.iloc[:,i].notna()&Ores.iloc[:,j].notna()&Pres.iloc[:,i].notna()&Pres.iloc[:,j].notna()
        if m.sum()<150: continue
        CO[i,j]=CO[j,i]=np.corrcoef(Ores.iloc[:,i][m],Ores.iloc[:,j][m])[0,1]
        CP[i,j]=CP[j,i]=np.corrcoef(Pres.iloc[:,i][m],Pres.iloc[:,j][m])[0,1]
        NP[i,j]=NP[j,i]=m.sum()
iu=np.triu_indices(k,1); ok=~np.isnan(CO[iu])&~np.isnan(CP[iu])
co,cp=CO[iu][ok],CP[iu][ok]
r=np.corrcoef(co,cp)[0,1]
print(f"\ncategory pairs with n>=150 : {ok.sum()}   median pair n = {int(np.median(NP[iu][ok]))}")
print(f"mean within-person onset correlation across categories : {np.nanmean(co):+.3f}")
print(f"mean residual preference correlation                   : {np.nanmean(cp):+.3f}")
print(f"\nRSA  corr(ACQUIRED-TOGETHER, LIKED-TOGETHER) = {r:+.3f}")
# null: permute category labels of one matrix
nl=[]
for _ in range(2000):
    p=rng.permutation(k); M=CP[np.ix_(p,p)]
    v=M[iu][ok]; g=~np.isnan(v)
    nl.append(np.corrcoef(co[g],v[g])[0,1])
nl=np.array(nl)
print(f"label-permutation null: mean {nl.mean():+.3f}  sd {nl.std():.3f}  p(|null|>=|obs|) = {np.mean(np.abs(nl)>=abs(r)):.4f}")
np.save('data/derived/CO.npy',CO); np.save('data/derived/CP.npy',CP)
pd.Series(names).to_csv('data/derived/onset_cats.csv',index=False,header=['cat'])
print("\ntop 'acquired together' pairs:")
o=np.argsort(-CO[iu][ok]); a,b=iu[0][ok][o],iu[1][ok][o]
for x,y in list(zip(a,b))[:8]: print(f"   r_onset={CO[x,y]:+.3f}  r_pref={CP[x,y]:+.3f}   {names[x]} + {names[y]}")
