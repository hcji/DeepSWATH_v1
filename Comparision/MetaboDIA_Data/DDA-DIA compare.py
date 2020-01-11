# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 08:10:25 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm
'''
f_deepdia = 'D:/MetaboDIA_data/PH/PH697021_neg_SWATH.ms2.csv'
f_msdial = 'D:/MetaboDIA_data/MS-DIAL Result/PH697021_neg_SWATH.csv'
f_dda = 'D:/MetaboDIA_data/PH/PH697021_neg_IDA.ms2.csv'
'''
def DDA_DIA_compare(f_dda, f_deepdia, f_msdial, mztol=0.01, rttol=10):
    dda_res = pd.read_csv(f_dda)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    
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
        # for fair comparesion, exclude low intensity peaks in msdial
        msdial = msdial[msdial['intensity'] > 100]
        
        mzs = list(set( list(np.round(dda['mz'], 2)) ))
        if len(mzs) < 3:
            continue
        
        dda_int, deepdia_int, msdial_int = [], [], []
        for mz in mzs:
            dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.01 ]))
            deepdia_int.append(np.sum(deepdia['intensity'][ np.abs(deepdia['mz'] - mz) < 0.01 ]))
            msdial_int.append(np.sum(msdial['intensity'][ np.abs(msdial['mz'] - mz) < 0.01 ]))
        deepdia_corr = pearsonr(dda_int, deepdia_int)[0]
        msdial_corr = pearsonr(dda_int, msdial_int)[0]
        output.loc[i] = [f['precursor_mz'], f['precursor_rt'], f['precursor_intensity'], deepdia_corr, msdial_corr]