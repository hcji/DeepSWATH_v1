# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 07:37:03 2020

@author: yn
"""

import pymzml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from tensorflow.keras.models import load_model
from sklearn.preprocessing import normalize
from DeepSWATH.utils import extract_eic, fragment_eic, get_ms2
from DeepSWATH.utils import parser_mzxml, parser_mzml

def DeepSWATH_process(file, features, min_int=150):
    # file = 'Comparision/MetDIA_Data/data/30STD_mix 330ppb-1.mzML'
    # features = pd.read_csv('Comparision/MetDIA_data/results/xcms_ms1_feature.csv')
    # features = features[['mz', 'rt', 'maxo']]
    # features.columns = ['mz', 'rt', 'intensity']
    features = features.drop_duplicates()
    mod = load_model('Model/DeepSWATH_Model_with_simu.h5')
    
    if file.split('.')[-1] == 'mzXML':
        reader = parser_mzxml(file)
    elif file.split('.')[-1] == 'mzML':
        reader = parser_mzml(file)
    else:
        raise IOError ('invalid input file')
    
    precursors = [-1.0]
    for p in reader:
        try:
            pmz = p.getPrecursors()[0].getMZ()
        except:
            continue
        if pmz < precursors[-1]:
            break
        else:
            precursors.append(pmz)
    precursors = np.array(precursors)
    
    all_peaks = dict(zip(precursors, [[]]*len(precursors)))
    for p in reader:
        if p.getMSLevel() == 1:
            all_peaks[-1.0].append(p)
        else:
            j = precursors[np.argmin(np.abs(p.getPrecursors()[0].getMZ() - precursors))]
            all_peaks[j].append(p)
    del(reader)
    
    output = []
    # output = pd.DataFrame(columns = ['exid', 'precursor_mz', 'precursor_rt', 'precursor_intensity', 'mz', 'intensity'])
    # output = open (output_path, 'a+')
    for exid in tqdm(features.index):
        exrt = features['rt'][exid]
        exmz = features['mz'][exid]
        exint = features['intensity'][exid]
        exeic = extract_eic(all_peaks[-1.0], exmz, exrt, rtlength=30)
        # plt.plot(exeic[0], exeic[1])

        k = precursors[np.argmin(np.abs(exmz - precursors))]
        peaks = all_peaks[k]
        
        ms2_i = get_ms2(peaks, precursors, exmz, exrt)      
        cid = np.where(np.logical_and(np.abs(exmz - ms2_i[0]) > -2.008, ms2_i[1] > min_int))[0]       
        candidate_mz, candidate_abund = ms2_i[0][cid], ms2_i[1][cid]
        if len(candidate_mz) < 1:
            continue
        
        temp_frag_mz, temp_frag_abund = [], []
        temp_precursor_eics, temp_fragment_eics = [], []
        
        for j in range(len(candidate_mz)):
            fragmz = candidate_mz[j]
            fragabund = candidate_abund[j]
            # remove very near mz
            if len(temp_frag_mz) > 0:
                if fragmz - temp_frag_mz[-1] < 0.025:
                    continue
            
            # fragabund = candidate_abund[j]
            frageic = fragment_eic(peaks, precursors, exmz, exrt, fragmz, rtlength=35)
            if len(np.where(np.array(frageic[1]) > 0)[0]) <= 3:
                continue
            # plt.plot(frageic[0], frageic[1])
            
            std_rt = np.linspace(exeic[0][0], exeic[0][-1], 100)
            std_ex = np.interp(std_rt, exeic[0], exeic[1])
            std_fg = np.interp(std_rt, frageic[0], frageic[1])
            
            temp_frag_mz.append(fragmz)
            temp_frag_abund.append(fragabund)
            temp_precursor_eics.append(std_ex)
            temp_fragment_eics.append(std_fg)
        
        if len(temp_frag_mz) == 0:
            continue
        
        X1 = np.asarray(temp_precursor_eics)
        X2 = np.asarray(temp_fragment_eics)
    
        X1 = normalize(X1, axis=1, norm='max')
        X2 = normalize(X2, axis=1, norm='max')
        X1 = np.expand_dims(X1,-1)
        X2 = np.expand_dims(X2,-1)
    
        prediction = mod.predict([X1, X2])
        pos = np.where(prediction[:,1] > 0.5)[0]
        '''
        plt.figure(dpi = 300)
        plt.plot(temp_precursor_eics[0], color = 'blue')
        for a in range(len(temp_fragment_eics)):
            plt.plot(temp_fragment_eics[a], color = 'gray')
        
        for p in pos:
            plt.plot(temp_fragment_eics[p], color='red')
        plt.show()
        '''
        temp_frag_mz = np.asarray(temp_frag_mz)[pos]
        temp_frag_abund = np.asarray(temp_frag_abund)[pos]
        temp_output = pd.DataFrame({'exid': exid, 'precursor_mz': exmz, 'precursor_rt': exrt, 'precursor_intensity': exint,
                                    'mz': temp_frag_mz, 'intensity': temp_frag_abund})
        output.append(temp_output)
        # output.write(str(temp_output))
    output = pd.concat(output)
    # output.close()
    return output
