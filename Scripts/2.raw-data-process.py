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
    swath = [f for f in files if 'SWATH.mzXML' in f]
    
    all_precursor_eics = []
    all_fragment_eics = []
    all_decoy_eics = []
    
    for f in swath:
        f = 'D:/MetaboDIA_data/CS/' + f
        x = f.replace('SWATH.mzXML', 'IDA.ms2.csv')
        y = f.replace('SWATH.mzXML', 'SWATH.feature.csv')
        features = pd.read_csv(x)
        swathfet = pd.read_csv(y)
        peaks = parser_mzxml(f)
        
        peaks1 = [p for p in peaks if p.getMSLevel()==1]
        peaks2 = [p for p in peaks if p.getMSLevel()==2]
        precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
        del(peaks)
        
        exdda = features[['precursor_mz', 'precursor_rt']]
        exdda = exdda.drop_duplicates()
        
        for i in tqdm(exdda.index):
            exrtdda = exdda['precursor_rt'][i]
            exmzdda = exdda['precursor_mz'][i]
            
            exdia = swathfet[np.abs(swathfet['mz'] - exmzdda) < 0.01]
            if len(exdia) < 1:
                continue
            if min(np.abs(exdia['rt'] - exrtdda)) > 5:
                continue
            else:
                exrt = exdia['rt'][np.argmin(np.abs(exdia['rt'] - exrtdda))]
                exmz = exdia['mz'][np.argmin(np.abs(exdia['rt'] - exrtdda))]
            
            frags = features[features['precursor_rt']==exrtdda]
            frags = frags[frags['precursor_mz']==exmzdda]
            
            ms2 = get_ms2(peaks2, precursors, exmz, exrt)
            cid = np.where(np.logical_and(np.abs(exmz - ms2[0]) > 0, ms2[1] > 100))[0]
            if len(cid) == 0:
                continue
            candidate_mz, candidate_abund = ms2[0][cid], ms2[1][cid]
            decoy_mzs = [m for m in candidate_mz if np.min(np.abs(m - frags['mz'])) > 5]
            exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)
            # plt.plot(exeic[0], exeic[1])
            
            for j in frags.index:
                fragmz = frags['mz'][j]
                if np.min(np.abs(fragmz - candidate_mz)) > 0.05:
                    continue
                
                abund = frags['intensity'][j]
                if abund < 200:
                    continue
                
                if len(decoy_mzs) > 0:
                    decoy_mz = choice(decoy_mzs)
                    decoy_mzs.remove(decoy_mz)
                else:
                    break
                    
                frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=35)
                decoyeic = fragment_eic(peaks2, precursors, exmz, exrt, decoy_mz, rtlength=35)
                # plt.plot(frageic[0], frageic[1])
                # plt.plot(decoyeic[0], decoyeic[1])
            
                std_rt = np.linspace(exeic[0][0], exeic[0][-1], 100)
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
    all_decoy_eics = np.asarray(all_decoy_eics)
    np.save('Data/all_precursor_eics.npy', all_precursor_eics)
    np.save('Data/all_fragment_eics.npy', all_fragment_eics)
    np.save('Data/all_decoy_eics.npy', all_decoy_eics)