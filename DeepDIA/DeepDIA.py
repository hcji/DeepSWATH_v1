# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 07:37:03 2020

@author: yn
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from tensorflow.keras.models import load_model
from sklearn.preprocessing import normalize
from DeepDIA.utils import parser_mzxml, extract_eic, fragment_eic, get_ms2

def DeepDIA_process(file, features, noise=200):
    # file = 'Example/CS52684_neg_SWATH.mzXML'
    # features = pd.read_csv('Example/CS52684_neg_SWATH.features.csv')
    mod = load_model('Model/DeepDIA_Model.h5')
    peaks = parser_mzxml(file)
    all_precursor_mz, all_precursor_rt = [], []
    all_exid, all_frag_mz, all_frag_abund = [], [], []
    all_precursor_eics, all_fragment_eics = [], []
    
    peaks1 = [p for p in peaks if p.getMSLevel()==1]
    peaks2 = [p for p in peaks if p.getMSLevel()==2]
    precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
    del(peaks)
    
    for i in tqdm(features.index):
        exid = features.iloc[i, 0]
        exrt = features['rt'][i]
        exmz = features['mz'][i]
        exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)
        # plt.plot(exeic[0], exeic[1])
        
        ms2 = get_ms2(peaks2, precursors, exmz, exrt)
        cid = np.where(np.logical_and(np.abs(exmz - ms2[0]) > 0, ms2[1] > noise))[0]
        candidate_mz, candidate_abund = ms2[0][cid], ms2[1][cid]
        if len(candidate_mz) < 1:
            continue
        for j in range(len(candidate_mz)):
            fragmz = candidate_mz[j]
            fragabund = candidate_abund[j]
            frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=35)
            # plt.plot(frageic[0], frageic[1])
            
            std_rt = np.linspace(exeic[0][0], exeic[0][-1], 50)
            std_ex = np.interp(std_rt, exeic[0], exeic[1])
            std_fg = np.interp(std_rt, frageic[0], frageic[1])
            
            all_exid.append(exid)
            all_precursor_mz.append(exmz)
            all_precursor_rt.append(exrt)
            all_frag_mz.append(fragmz)
            all_frag_abund.append(fragabund)
            all_precursor_eics.append(std_ex)
            all_fragment_eics.append(std_fg)
    X1 = np.asarray(all_precursor_eics)
    X2 = np.asarray(all_fragment_eics)
    
    X1 = normalize(X1, axis=1, norm='max')
    X2 = normalize(X2, axis=1, norm='max')
    X1 = np.expand_dims(X1,-1)
    X2 = np.expand_dims(X2,-1)
    
    prediction = mod.predict([X1, X2])
    pos = np.where(prediction[:,1] > 0.5)[0]
    # plt.plot(all_precursor_eics[pos[21]])
    # plt.plot(all_fragment_eics[pos[21]])
    output = pd.DataFrame({'exid': np.asarray(all_exid)[pos], 'precursor_mz': np.asarray(all_precursor_mz)[pos], 
                           'precursor_rt': np.asarray(all_precursor_rt)[pos], 'mz': np.asarray(all_frag_mz)[pos],
                           'intensity': np.asarray(all_frag_abund)[pos]})
    return output
        