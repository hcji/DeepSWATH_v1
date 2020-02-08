# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 08:10:25 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm


def count_features(ms2):
    features = pd.read_csv(ms2)
    features = features[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    return len(features)


def get_str_for_venn(ms2):
    features = pd.read_csv(ms2)
    features = features[features['intensity'] > 100]
    v = [str(round(features.loc[r,'precursor_mz'],2)) + ' ' +
         str(round(features.loc[r,'precursor_rt'],-1)) + ' ' +
         str(round(features.loc[r,'mz'],1)) for r in features.index] 
    return set(v)


def DB_DIA_compare(db, f_deepdia, f_msdial, mztol=0.05):
    db_res = pd.read_csv(db)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    
    features = db_res[['Name', 'PrecursorMz']]
    features = features.drop_duplicates()
    
    output = pd.DataFrame(columns=['Name', 'precursor_mz', 'precursor_intensity', 'DeepDIA_corr', 'MSDIAL_corr'])
    for i in tqdm(range(features.shape[0])):
        f = features.iloc[i,:]
        standard = db_res[db_res['PrecursorMz'] == f['PrecursorMz']]
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        precursor_intensity = np.mean(deepdia['precursor_intensity'])
        
        mzs = standard['ProductMz']
        std_int = np.array(standard['LibraryIntensity'])
        deepdia_int, msdial_int = [], []
        for mz in mzs:
            deepdia_int.append(np.sum(deepdia['intensity'][ np.abs(deepdia['mz'] - mz) < 0.1 ]))
            msdial_int.append(np.sum(msdial['intensity'][ np.abs(msdial['mz'] - mz) < 0.1 ]))
        deepdia_corr = np.nanmax([0, pearsonr(std_int, deepdia_int)[0]])
        msdial_corr = np.nanmax([0, pearsonr(std_int, msdial_int)[0]])
        output.loc[i] = [f['Name'], f['PrecursorMz'], precursor_intensity, deepdia_corr, msdial_corr]
    return output
        


def DDA_DIA_compare(f_dda, f_deepdia, f_msdial, mztol=0.05, rttol=60):
    dda_res = pd.read_csv(f_dda)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    # for fair comparesion, exclude low intensity peaks in msdial
    msdial_res = msdial_res[msdial_res['intensity'] > 100]
    
    features = dda_res[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    
    output = pd.DataFrame(columns=['precursor_mz', 'precursor_rt', 'precursor_intensity', 'DeepDIA_corr', 'MSDIAL_corr'])
    for i in tqdm(range(features.shape[0])):
        f = features.iloc[i,:]
        dda = dda_res[ np.abs(dda_res['precursor_mz'] - f['precursor_mz']) < mztol ]
        dda = dda[ np.abs(dda['precursor_rt'] - f['precursor_rt']) < rttol ]
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['precursor_mz']) < mztol ]
        deepdia = deepdia[ np.abs(deepdia['precursor_rt'] - f['precursor_rt']) < rttol ]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['precursor_mz']) < mztol ]
        msdial = msdial[ np.abs(msdial['precursor_rt'] - f['precursor_rt']) < rttol ]
        
        if dda.shape[0] <= 3:
            continue  
        mzs = list(set( list(np.round(dda['mz'], 1)) + list(np.round(msdial['mz'], 1)) + list(np.round(deepdia['mz'], 1)) ))      
        dda_int, deepdia_int, msdial_int = [], [], []
        for mz in mzs:
            dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.15 ]))
            deepdia_int.append(np.sum(deepdia['intensity'][ np.abs(deepdia['mz'] - mz) < 0.15 ]))
            msdial_int.append(np.sum(msdial['intensity'][ np.abs(msdial['mz'] - mz) < 0.15 ]))
        deepdia_corr = np.nanmax([0, pearsonr(dda_int, deepdia_int)[0]])
        msdial_corr = np.nanmax([0, pearsonr(dda_int, msdial_int)[0]])
        output.loc[i] = [f['precursor_mz'], f['precursor_rt'], f['precursor_intensity'], deepdia_corr, msdial_corr]
    return output

        
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    # MetDIA Comparision
    db = 'Comparision/MetDIA_data/results/30STD_Database.csv'
    f_deepdia = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_DeepDIA.csv'
    f_msdial = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_MSDIAL.csv'
    metdia = DB_DIA_compare(db, f_deepdia, f_msdial, mztol=0.05)
    
    # MetaboDIA Comparision
    p_dda = 'Comparision/MetaboDIA_data/results/PH697097_pos_IDA.csv'
    n_dda = 'Comparision/MetaboDIA_data/results/PH697097_neg_IDA.csv'
    p_deepdia = 'Comparision/MetaboDIA_data/results/PH697097_pos_DeepDIA.csv'
    n_deepdia = 'Comparision/MetaboDIA_data/results/PH697097_neg_DeepDIA.csv'
    p_msdial = 'Comparision/MetaboDIA_data/results/PH697097_pos_MSDIAL.csv'
    n_msdial = 'Comparision/MetaboDIA_data/results/PH697097_neg_MSDIAL.csv'
    
    p_metabodia = DDA_DIA_compare(p_dda, p_deepdia, p_msdial)
    n_metabodia = DDA_DIA_compare(n_dda, n_deepdia, n_msdial)

    # Correlation Violin Plot
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 6))
    axes[0,0].violinplot([list(metdia['DeepDIA_corr'])], [1], showmeans=False, showmedians=True)
    axes[0,0].violinplot([list(metdia['MSDIAL_corr'])], [2], showmeans=False, showmedians=True)
    axes[0,0].set_xticks(range(4))
    axes[0,0].set_xticklabels(['', 'DeepEI', 'MSDIAL', ''])
    axes[0,0].set_ylabel('Correlation')
      
    axes[0,1].violinplot([list(p_metabodia['DeepDIA_corr']), list(n_metabodia['DeepDIA_corr'])], [1,5], showmeans=False, showmedians=True)
    axes[0,1].violinplot([list(p_metabodia['MSDIAL_corr']), list(n_metabodia['MSDIAL_corr'])], [3,7], showmeans=False, showmedians=True)
    axes[0,1].set_xticks(range(8))
    axes[0,1].set_xticklabels(['', 'DeepDIA', '\nPositive', 'MSDIAL', '', 'DeepDIA', '\nNegative', 'MSDIAL'])
    axes[0,1].set_ylabel('Correlation')
    
    exp_mz = 626.356
    exp_rt = 739.134
    mztol = 0.01
    rttol = 15
    dda_res = pd.read_csv(p_dda)
    deepdia_res = pd.read_csv(p_deepdia)
    msdial_res = pd.read_csv(p_msdial) 
    
    dda = dda_res[ np.abs(dda_res['precursor_mz'] - exp_mz) < mztol ]
    dda = dda[ np.abs(dda['precursor_rt'] - exp_rt) < rttol ]
    deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - exp_mz) < mztol ]
    deepdia = deepdia[ np.abs(deepdia['precursor_rt'] - exp_rt) < rttol ]
    msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - exp_mz) < mztol ]
    msdial = msdial[ np.abs(msdial['precursor_rt'] - exp_rt) < rttol ]
    
    axes[1,0].vlines(dda['mz'], 0, dda['intensity'] / np.max(dda['intensity']), color='red', alpha=0.7, label='DDA')
    axes[1,0].vlines(deepdia['mz'], 0, -deepdia['intensity'] / np.max(deepdia['intensity']), color='blue', alpha=0.7, label='DeepDIA')
    axes[1,0].text(80, 0.8, 'MS/MS \n precursor:'+str(exp_mz))
    axes[1,0].axhline(0, color='black')
    axes[1,0].set_xlabel('m/z')
    axes[1,0].set_ylabel('Abundance')
    axes[1,0].legend()
    
    axes[1,1].vlines(dda['mz'], 0, dda['intensity']/ np.max(dda['intensity']), color='red', alpha=0.8, label='DDA')
    axes[1,1].vlines(msdial['mz'], 0, -msdial['intensity']/ np.max(msdial['intensity']), color='orange', alpha=0.8, label='MS-DIAL')
    axes[1,1].text(60, 0.8, 'MS/MS \nprecursor:'+str(exp_mz))
    axes[1,1].axhline(0, color='black')
    axes[1,1].set_xlabel('m/z')
    axes[1,1].set_ylabel('Abundance')
    axes[1,1].legend()
    
    plt.show()
    