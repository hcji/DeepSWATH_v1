# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 16:26:17 2020

@author: hcji
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import pearsonr

comp = pd.read_csv('Comparision/MetDIA_Data/results/comparision.csv')
db = 'Comparision/MetDIA_data/results/30STD_Database.csv'
f_deepdia = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_DeepDIA.csv'
f_msdial = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_MSDIAL.csv'

def split_precursor(mslist):
    f1 = list(set(zip(mslist['precursor_mz'], mslist['precursor_rt'])))
    output = []
    for i in f1:
        ms = mslist[mslist['precursor_mz'] == i[0]]
        ms = ms[ms['precursor_rt'] == i[1]]
        output.append(ms)
    return output


def ms_output(db, f_deepdia, f_msdial, mztol=0.05):
    db_res = pd.read_csv(db)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    
    features = db_res[['Name', 'PrecursorMz']]
    features = features.drop_duplicates()
    
    for i in tqdm(range(features.shape[0])):
        f = features.iloc[i,:]
        standard = db_res[db_res['PrecursorMz'] == f['PrecursorMz']]
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        
        deepdia = split_precursor(deepdia)
        msdial = split_precursor(msdial)
        
        mzs = standard['ProductMz']
        std_int = np.array(standard['LibraryIntensity'])
        
        deepdia_corrs, msdial_corrs = [], []
        for j in range(len(deepdia)):
            deepdia_int = []
            for mz in mzs:
                deepdia_int.append(np.sum(deepdia[j]['intensity'][ np.abs(deepdia[j]['mz'] - mz) < 0.1 ]))
            deepdia_corrs.append(np.nanmax([0, pearsonr(std_int, deepdia_int)[0]]))
        
        for j in range(len(msdial)):
            msdial_int = []
            for mz in mzs:
                msdial_int.append(np.sum(msdial[j]['intensity'][ np.abs(msdial[j]['mz'] - mz) < 0.1 ]))
            msdial_corrs.append(np.nanmax([0, pearsonr(std_int, msdial_int)[0]]))

        w1 = np.argmax(deepdia_corrs)
        w2 = np.argmax(msdial_corrs)
        
        deepdia = deepdia[w1][['mz', 'intensity']]
        msdial = msdial[w2][['mz', 'intensity']]
               
        plt.vlines(deepdia['mz'], 0, deepdia['intensity'] / np.max(deepdia['intensity']), color='red', alpha=0.7, label='DeepSWATH')
        plt.vlines(msdial['mz'], 0, -msdial['intensity'] / np.max(msdial['intensity']), color='blue', alpha=0.7, label='MSDIAL')
        plt.axhline(0, color='black')
        plt.xlabel('m/z')
        plt.ylabel('Abundance')
        plt.legend()
        
        deepdia.to_csv('MS/'+f['Name']+'_deepdia.csv', index=False, header=False)
        msdial.to_csv('MS/'+f['Name']+'_msdial.csv', index=False, header=False)
        
        
if __name__ == '__main__':
        
    ms_output(db, f_deepdia, f_msdial, mztol=0.05)
    