import pandas as pd, numpy as np, re
df=pd.read_csv('data/raw/BKSPublic.csv',low_memory=False)
inv=pd.read_csv('data/derived/inventory.csv')
ons=inv[inv['kind']=='AGE_ONSET']['col'].tolist()
BIN={'0-4yo':2,'5-6yo':5.5,'7-8yo':7.5,'9-10yo':9.5,'11-12yo':11.5,'13-14yo':13.5,
     '15-16yo':15.5,'17-18yo':17.5,'19-25yo':22,'26yo+':28}
O=pd.DataFrame({c:df[c].map(BIN) for c in ons})
print(f"onset columns: {len(ons)}")
print(f"respondents with >=2 onsets : {(O.notna().sum(1)>=2).sum():,}")
print(f"respondents with >=5 onsets : {(O.notna().sum(1)>=5).sum():,}")
print(f"median onsets per person    : {O.notna().sum(1).median():.0f}")
print("\nper-column: n, mean onset")
s=pd.DataFrame({'n':O.notna().sum(),'mean':O.mean().round(1)}).sort_values('mean')
s.index=[re.sub(r'\s*\(\w+\)$','',i)[:78] for i in s.index]
print(s.to_string())
O.to_csv('data/derived/onset.csv',index=False)
