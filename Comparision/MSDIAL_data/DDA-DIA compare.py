# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:04:52 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr


def DDA_DIA_compare(f_dda, f_deepdia, f_msdial, mztol=0.05, rttol=20):
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
    
    n_dda = 'Comparision/MSDIAL_data/results/MSDIAL_results/Nega_Ida_QC_1_1.csv'
    n_deepdia = 'Comparision/MSDIAL_data/results/DeepDIA_results/Nega_Swath_QC_1_1.ms2.csv'
    n_msdial = 'Comparision/MSDIAL_data/results/MSDIAL_results/Nega_Swath_QC_1_1.csv'
    

    n_dda_dia = DDA_DIA_compare(n_dda, n_deepdia, n_msdial)
    
    