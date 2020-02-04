# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 08:10:25 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm


def get_dda_consensus(dda_xcms, dda_msdial):
    xcms = pd.read_csv(dda_xcms)
    msdial = pd.read_csv(dda_msdial)
    xcms_fet = xcms[['precursor_mz', 'precursor_rt', 'precursor_intensity']].drop_duplicates()
    msdial_fet = msdial[['precursor_mz', 'precursor_rt', 'precursor_intensity']].drop_duplicates()
    consensus_fet = pd.DataFrame(columns=xcms_fet.columns)
    for i in xcms_fet.index:
        pmz = xcms_fet['precursor_mz'][i]
        prt = xcms_fet['precursor_rt'][i]
        pint = xcms_fet['precursor_intensity'][i]
        for j in msdial_fet.index:
            if (abs(msdial_fet['precursor_mz'][j] - pmz) < 0.01) and (abs(msdial_fet['precursor_rt'][j] - prt) < 10):
                consensus_fet.loc[len(consensus_fet)] = [pmz, prt, pint]
    consensus = pd.DataFrame(columns=['precursor_mz', 'precursor_rt', 'precursor_intensity', 'mz', 'intensity'])
    for i in consensus_fet.index:
        pmz = consensus_fet['precursor_mz'][i]
        prt = consensus_fet['precursor_rt'][i]
        pint = consensus_fet['precursor_intensity'][i]
        for j in xcms.index:
            if (xcms['precursor_mz'][j] == pmz) and (xcms['precursor_rt'][j] == prt):
                consensus.loc[len(consensus)] = [pmz, prt, pint, xcms['mz'][j], xcms['intensity'][j]]
    return consensus         


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


def DDA_DIA_compare(f_dda, f_deepdia, f_msdial, mztol=0.01, rttol=30):
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

        mzs = list(set( list(np.round(dda['mz'], 1)) + list(np.round(msdial['mz'], 1)) + list(np.round(deepdia['mz'], 1)) ))
        if len(mzs) < 3:
            continue
        
        dda_int, deepdia_int, msdial_int = [], [], []
        for mz in mzs:
            dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.1 ]))
            deepdia_int.append(np.sum(deepdia['intensity'][ np.abs(deepdia['mz'] - mz) < 0.1 ]))
            msdial_int.append(np.sum(msdial['intensity'][ np.abs(msdial['mz'] - mz) < 0.1 ]))
        deepdia_corr = np.nanmax([0, pearsonr(dda_int, deepdia_int)[0]])
        msdial_corr = np.nanmax([0, pearsonr(dda_int, msdial_int)[0]])
        output.loc[i] = [f['precursor_mz'], f['precursor_rt'], f['precursor_intensity'], deepdia_corr, msdial_corr]
    return output

        
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    
    p_dda_msdial = 'Comparision/MetaboDIA_Data/results/DDA_results/MS_DIAL_PH697097_pos_IDA.csv'
    p_dda_xcms = 'Comparision/MetaboDIA_Data/results/DDA_results/XCMS_PH697097_pos_IDA.csv'
    n_dda_msdial = 'Comparision/MetaboDIA_Data/results/DDA_results/MS_DIAL_PH697097_neg_IDA.csv'
    n_dda_xcms = 'Comparision/MetaboDIA_Data/results/DDA_results/XCMS_PH697097_neg_IDA.csv'
    p_dda_consensus = get_dda_consensus(p_dda_xcms, p_dda_msdial)
    n_dda_consensus = get_dda_consensus(n_dda_xcms, n_dda_msdial)
    p_dda_consensus.to_csv('Comparision/MetaboDIA_data/results/DDA_results/PH697097_pos_IDA.ms2.csv')
    n_dda_consensus.to_csv('Comparision/MetaboDIA_data/results/DDA_results/PH697097_neg_IDA.ms2.csv')
    
    p_dda = 'Comparision/MetaboDIA_data/results/DDA_results/PH697097_pos_IDA.ms2.csv'
    n_dda = 'Comparision/MetaboDIA_data/results/DDA_results/PH697097_neg_IDA.ms2.csv'
    p_deepdia = 'Comparision/MetaboDIA_data/results/DeepDIA_results/PH697097_pos_SWATH.ms2.csv'
    n_deepdia = 'Comparision/MetaboDIA_data/results/DeepDIA_results/PH697097_neg_SWATH.ms2.csv'
    p_msdial = 'Comparision/MetaboDIA_data/results/MSDIAL_results/PH697097_pos_SWATH.csv'
    n_msdial = 'Comparision/MetaboDIA_data/results/MSDIAL_results/PH697097_neg_SWATH.csv'
    
    p_dda_dia = DDA_DIA_compare(p_dda, p_deepdia, p_msdial)
    n_dda_dia = DDA_DIA_compare(n_dda, n_deepdia, n_msdial)
    
    num_p_dda, num_n_dda = count_features(p_dda), count_features(n_dda)
    num_p_deepdia, num_n_deepdia = count_features(p_deepdia), count_features(n_deepdia)
    num_p_msdial, num_n_msdial = count_features(p_msdial), count_features(n_msdial)
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        
    axes[0,0].bar([1,5], [num_p_dda, num_n_dda], width=1, label='DDA', color='red', alpha=0.5)
    axes[0,0].bar([2,6], [num_p_deepdia, num_n_deepdia], width=1, label='DeepDIA', alpha=0.5)
    axes[0,0].bar([3,7], [num_p_msdial, num_n_msdial], width=1, label='MS-DIAL', alpha=0.5)
    axes[0,0].set_xticklabels(['', '', 'Positive', '', '','','Negative'])
    axes[0,0].set_yticks(np.arange(0,14001,2000))
    axes[0,0].set_ylabel('Number')
    axes[0,0].legend()
    
    axes[0,1].violinplot([list(p_dda_dia['DeepDIA_corr']), list(n_dda_dia['DeepDIA_corr'])], [1,5], showmeans=False, showmedians=True)
    axes[0,1].violinplot([list(p_dda_dia['MSDIAL_corr']), list(n_dda_dia['MSDIAL_corr'])], [3,7], showmeans=False, showmedians=True)
    axes[0,1].set_xticks(range(8))
    axes[0,1].set_xticklabels(['', 'DeepEI', '\nPositive', 'MSDIAL', '', 'DeepEI', '\nNegative', 'MSDIAL'])
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
    
    