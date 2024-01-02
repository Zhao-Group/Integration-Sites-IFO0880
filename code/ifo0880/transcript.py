#Modules
import numpy as np
import pandas as pd

#Data
df = pd.read_csv('data/ifo0880/output_10_5.csv')
rna_seq_raw = pd.read_excel('data/ifo0880/253_2021_11549_MOESM2_ESM.xlsx', sheet_name = 'All data')

#Annotation consistency
rna_seq_raw['Synonym'] = 'AAT19DRAFT_' + rna_seq_raw['Protein ID'].astype(str)

#Mean normalized expression value log2(counts per million)
rna_seq_raw['Ac'] = rna_seq_raw[['Rhoto.Ac.1', 'Rhoto.Ac.2', 'Rhoto.Ac.3']].mean(axis=1)
rna_seq_raw['Glc'] = rna_seq_raw[['Rhoto.Glc.1', 'Rhoto.Glc.2', 'Rhoto.Glc.3']].mean(axis=1)
rna_seq_raw['So'] = rna_seq_raw[['Rhoto.So.1', 'Rhoto.So.2', 'Rhoto.So.3']].mean(axis=1)
rna_seq_raw['Xyl'] = rna_seq_raw[['Rhoto.Xyl.1', 'Rhoto.Xyl.2', 'Rhoto.Xyl.3']].mean(axis=1)
rna_seq_raw['YP'] = rna_seq_raw[['Rhoto.YP.1', 'Rhoto.YP.2', 'Rhoto.YP.3']].mean(axis=1)

rna_seq = rna_seq_raw[['Synonym', 'FC.Glc_vs_YP', 'FC.Xyl_vs_Glc', 'FC.Ac_vs_Glc', 'FC.So_vs_Glc', 'Ac', 'Glc', 'So', 'Xyl', 'YP']]

##Neutrality condition: Stability (2 > FC > -2) and adjacent gene expression (growth on glucose) > 1st quartile 
stable_genes = list(rna_seq[
    (rna_seq['FC.Glc_vs_YP'].between(-2, 2)) &
    (rna_seq['FC.Xyl_vs_Glc'].between(-2, 2)) &
    (rna_seq['FC.Ac_vs_Glc'].between(-2, 2)) &
    (rna_seq['FC.So_vs_Glc'].between(-2, 2))
]['Synonym'])

left_tr = []
right_tr = []
left_stab = []
right_stab = []
for i in range(np.shape(df)[0]):
    #Adding transcriptomic information in Glucose media
    l_tr = rna_seq.loc[rna_seq['Synonym'] == df['Left Gene'][i], 'Glc'] #Glucose media for mean normalized value
    r_tr = rna_seq.loc[rna_seq['Synonym'] == df['Right Gene'][i], 'Glc'] #Glucose media for mean normalized value
    
    if not l_tr.empty:
        left_tr.append(l_tr.iloc[0])
    else:
        left_tr.append(np.nan)
        
    if not r_tr.empty:
        right_tr.append(r_tr.iloc[0])
    else:
        right_tr.append(np.nan)

    #Adding stability index (0: Unstable, 1: Stable)
    if df['Left Gene'][i] in stable_genes:
        left_stab.append(1)
    else:
        left_stab.append(0)

    if df['Right Gene'][i] in stable_genes:
        right_stab.append(1)
    else:
        right_stab.append(0)
        
df['Left Gene Expr'] = left_tr
df['Right Gene Expr'] = right_tr
df['Left Gene Stab'] = left_stab
df['Right Gene Stab'] = right_stab

cutoff_n = np.percentile(rna_seq.Glc, 25) 
print(cutoff_n)
df_neutral = df.loc[(df['Left Gene Expr'] >= cutoff_n) & df['Left Gene Expr'].notna() &
                       (df['Right Gene Expr'] > cutoff_n) & df['Right Gene Expr'].notna() & 
                       (df['Left Gene Stab'] == 1) & 
                       (df['Right Gene Stab'] == 1)].reset_index(drop=True)

pd.DataFrame(df_neutral).to_csv('data/ifo0880/selected_sites_neutral_10_5.csv', index = False)

##Expression condition: Adjacent gene expression (growth on glucose) > Median 
cutoff_e = np.percentile(rna_seq.Glc, 75)
print(cutoff_e)
df_express = df.loc[(df['Left Gene Expr'] >= cutoff_e) & df['Left Gene Expr'].notna() &
                       (df['Right Gene Expr'] > cutoff_e) & df['Right Gene Expr'].notna()].reset_index(drop=True)
pd.DataFrame(df_express).to_csv('data/ifo0880/selected_sites_expression_10_5.csv', index = False)
