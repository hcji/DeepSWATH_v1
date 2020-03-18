# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 16:26:17 2020

@author: hcji
"""

import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import pearsonr
from PyCFMID.PyCFMID import cfm_id_database


comp = pd.read_csv('Comparision/MetaboDIA_Data/results/comparision_pos.csv')
chebi = pd.read_csv('Data/chebi.csv')
f_dda = 'Comparision/MetaboDIA_data/results/PH697097_pos_IDA.csv'
f_deepdia = 'Comparision/MetaboDIA_data/results/PH697097_pos_DeepDIA.csv'
f_msdial = 'Comparision/MetaboDIA_data/results/PH697097_pos_MSDIAL.csv'


def split_precursor(mslist):
    f1 = list(set(zip(mslist['precursor_mz'], mslist['precursor_rt'])))
    output = []
    for i in f1:
        ms = mslist[mslist['precursor_mz'] == i[0]]
        ms = ms[ms['precursor_rt'] == i[1]]
        output.append(ms)
    return output


def DDA_DIA_identification(f_dda, f_deepdia, f_msdial, mztol=0.05, rttol=40):
    dda_res = pd.read_csv(f_dda)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    # for fair comparesion, exclude low intensity peaks in msdial
    msdial_res = msdial_res[msdial_res['intensity'] > 150]
    
    features = dda_res[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    
    precursor_mz, precursor_rt = [], []
    rank_1, rank_2 = [], []
    for i in tqdm(range(comp.shape[0])):
        f = comp.iloc[i,:]        
        standard = dda_res[dda_res['precursor_mz'] == f['precursor_mz']]
        
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['precursor_mz']) < mztol]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['precursor_mz']) < mztol]
        
        deepdia = split_precursor(deepdia)
        msdial = split_precursor(msdial)
        
        mzs = standard['precursor_mz']
        std_int = np.array(standard['intensity'])
        
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
        standard = standard[['mz', 'intensity']]
        
        if np.min(np.abs(chebi['Mass'] + 1.003 - f['precursor_mz'])) > mztol:
            continue
        else:
            candidates = chebi[np.abs(chebi['Mass'] + 1.003 - f['precursor_mz']) <= mztol]
        
        os.mkdir('Input')
        candidates[['ChEBI ID', 'SMILES']].to_csv('Input/candidate.txt', index=False, header=False, sep=' ')
        
        dda_id = cfm_id_database(standard, '', database=os.getcwd() + '/Input/candidate.txt', abs_mass_tol=mztol)['result']
        deepdia_id = cfm_id_database(deepdia, '', database=os.getcwd() + '/Input/candidate.txt', abs_mass_tol=mztol)['result']
        msdial_id = cfm_id_database(msdial, '', database=os.getcwd() + '/Input/candidate.txt', abs_mass_tol=mztol)['result']
        
        shutil.rmtree('Input')
        shutil.rmtree('Output')
        
        dda_comp = dda_id.iloc[0]['ID']
        deepdia_rank = len(np.where(deepdia_id['Score'] > list(deepdia_id['Score'])[list(deepdia_id['ID']).index(dda_comp)])[0]) + 1
        msdial_rank = len(np.where(msdial_id['Score'] > list(msdial_id['Score'])[list(msdial_id['ID']).index(dda_comp)])[0]) + 1
        
        rank_1.append(deepdia_rank)
        rank_2.append(msdial_rank)
        precursor_mz.append(f['precursor_mz'])
        precursor_rt.append(f['precursor_rt'])
        
        output = pd.DataFrame({'precursor_mz': precursor_mz, 'precursor_rt': precursor_rt, 'rank_deepdia': rank_1, 'rank_msdial': rank_2})
        output.to_csv('Comparision/MetaboDIA_Data/results/identification_pos.csv')
    return output


if __name__ == '__main__':
        
    output = DDA_DIA_identification(f_dda, f_deepdia, f_msdial)
    output.to_csv('Comparision/MetaboDIA_Data/results/identification_pos.csv')
    