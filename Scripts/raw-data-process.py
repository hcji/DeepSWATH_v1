# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 09:15:20 2019

@author: yn
"""

import numpy as np
import pandas as pd
from random import choice, randint
from DeepDIA.utils import parser_mzxml, extract_eic, fragment_eic, get_ms2

if __name__ == '__main__':
    
    import os
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    
    files = os.listdir('D:/MetaboDIA_data/CS')
    swath = [f for f in files if 'pos_SWATH.mzXML' in f]
    
    all_precursor_eics = []
    all_fragment_eics = []
    all_decoy_eics = []
    
    for f in swath:
        f = 'D:/MetaboDIA_data/CS/' + f
        x = f.replace('SWATH.mzXML', 'IDA.ms2.csv')
        features = pd.read_csv(x)
        peaks = parser_mzxml(f)
        
        peaks1 = [p for p in peaks if p.getMSLevel()==1]
        peaks2 = [p for p in peaks if p.getMSLevel()==2]
        precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
        del(peaks)
        
        exs = features[['precursor_mz', 'precursor_rt']]
        exs = exs.drop_duplicates()
        
        for i in tqdm(exs.index):
            exrt = exs['precursor_rt'][i]
            exmz = exs['precursor_mz'][i]
            
            frags = features[features['precursor_rt']==exrt]
            frags = frags[frags['precursor_mz']==exmz]
            
            ms2 = get_ms2(peaks2, precursors, exmz, exrt)
            cid = np.where(np.logical_and(np.abs(exmz - ms2[0]) > 0, ms2[1] > 200))[0]
            candidate_mz, candidate_abund = ms2[0][cid], ms2[1][cid]
            decoy_mzs = [m for m in candidate_mz if np.min(np.abs(m - frags['mz'])) > 5]
            exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)
            # plt.plot(exeic[0], exeic[1])
            
            for j in frags.index:
                fragmz = frags['mz'][j]
                abund = frags['intensity'][j]
                
                if abund < 50:
                    continue
                
                if len(decoy_mzs) > 0:
                    decoy_mz = choice(decoy_mzs)
                    decoy_mzs.remove(decoy_mz)
                else:
                    decoy_mz = randint(0, int(exmz*1000)) / 1000
                    
                frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=35)
                decoyeic = fragment_eic(peaks2, precursors, exmz, exrt, decoy_mz, rtlength=35)
                # plt.plot(frageic[0], frageic[1])
                # plt.plot(decoyeic[0], decoyeic[1])
            
                std_rt = np.linspace(exeic[0][0], exeic[0][-1], 50)
                std_ex = np.interp(std_rt, exeic[0], exeic[1])
                std_fg = np.interp(std_rt, frageic[0], frageic[1])
                std_dy = np.interp(std_rt, decoyeic[0], decoyeic[1])
                # plt.plot(std_rt, std_ex)
                # plt.plot(std_rt, std_fg)
            
                all_precursor_eics.append(std_ex)
                all_fragment_eics.append(std_fg)
                all_decoy_eics.append(std_dy)
        print(f + ' is finished')
    
    all_precursor_eics = np.asarray(all_precursor_eics)
    all_fragment_eics = np.asarray(all_fragment_eics)
    np.save('Data/all_precursor_eics.npy', all_precursor_eics)
    np.save('Data/all_fragment_eics.npy', all_fragment_eics)
    np.save('Data/all_decoy_eics.npy', all_decoy_eics)