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
from DeepDIA.utils import parser_mzxml, parser_mzml, extract_eic, fragment_eic, get_ms2

def DeepDIA_process(file, features, noise=100):
    # file = 'Comparision/MetDIA_Data/30STD_mix 330ppb-1.mzML'
    # features = pd.read_csv('Comparision/MetDIA_data/results/xcms_ms1_feature.csv')
    # features = features[['mz', 'rt', 'maxo']]
    # features.columns = ['mz', 'rt', 'intensity']
    mod = load_model('Model/DeepDIA_Model.h5')
    if file.split('.')[-1] == 'mzXML':
        peaks = parser_mzxml(file)
    elif file.split('.')[-1] == 'mzML':
        peaks = parser_mzml(file)
    else:
        raise IOError ('invalid input file')
    
    peaks1 = [p for p in peaks if p.getMSLevel()==1]
    peaks2 = [p for p in peaks if p.getMSLevel()==2]
    precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
    del(peaks)
    
    output = []
    # output = pd.DataFrame(columns = ['exid', 'precursor_mz', 'precursor_rt', 'precursor_intensity', 'mz', 'intensity'])
    # output = open (output_path, 'a+')
    for i in tqdm(range(len(features.index))):
        exid = features.index[i]
        exrt = features['rt'][exid]
        exmz = features['mz'][exid]
        exint = features['intensity'][exid]
        exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)
        # plt.plot(exeic[0], exeic[1])
        
        ms2 = get_ms2(peaks2, precursors, exmz, exrt)
        cid = np.where(np.logical_and(np.abs(exmz - ms2[0]) > 0, ms2[1] > noise))[0]
        candidate_mz, candidate_abund = ms2[0][cid], ms2[1][cid]
        if len(candidate_mz) < 1:
            continue
        
        temp_frag_mz, temp_frag_abund = [], []
        temp_precursor_eics, temp_fragment_eics = [], []
        
        for j in range(len(candidate_mz)):
            fragmz = candidate_mz[j]
            fragabund = candidate_abund[j]
            frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=35)
            # plt.plot(frageic[0], frageic[1])
            
            std_rt = np.linspace(exeic[0][0], exeic[0][-1], 50)
            std_ex = np.interp(std_rt, exeic[0], exeic[1])
            std_fg = np.interp(std_rt, frageic[0], frageic[1])
            
            temp_frag_mz.append(fragmz)
            temp_frag_abund.append(fragabund)
            temp_precursor_eics.append(std_ex)
            temp_fragment_eics.append(std_fg)
        
        X1 = np.asarray(temp_precursor_eics)
        X2 = np.asarray(temp_fragment_eics)
    
        X1 = normalize(X1, axis=1, norm='max')
        X2 = normalize(X2, axis=1, norm='max')
        X1 = np.expand_dims(X1,-1)
        X2 = np.expand_dims(X2,-1)
    
        prediction = mod.predict([X1, X2])
        pos = np.where(prediction[:,1] > 0.5)[0]
        # plt.plot(all_precursor_eics[pos[21]])
        # plt.plot(all_fragment_eics[pos[21]])
        temp_frag_mz = np.asarray(temp_frag_mz)[pos]
        temp_frag_abund = np.asarray(temp_frag_abund)[pos]
        temp_output = pd.DataFrame({'exid': exid, 'precursor_mz': exmz, 'precursor_rt': exrt, 'precursor_intensity': exint,
                                    'mz': temp_frag_mz, 'intensity': temp_frag_abund})
        output.append(temp_output)
        # output.write(str(temp_output))
    # output.close()
    return output
        