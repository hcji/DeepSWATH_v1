# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 08:04:29 2020

@author: hcji
"""

import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process

file = 'Example/CS52684_neg_SWATH.mzXML'

'''
# untargeted
features = pd.read_csv('Example/CS52684_neg_SWATH.features.csv', index_col=0)    
dia_result = DeepDIA_process(file, features)
dia_result.to_csv('Example/CS52684_neg_SWATH.ms2.csv')
'''

# targeted
dda_result = pd.read_csv('Example/CS52684_neg_IDA.ms2.csv')
features = dda_result[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
features = features.drop_duplicates()
features.columns = ['mz', 'rt', 'intensity']
dia_result = DeepDIA_process(file, features)

###############################################################################
# discussion 1
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.models import load_model
from DeepDIA.utils import parser_mzxml, extract_eic, fragment_eic, get_ms2
from DeepDIA.plot import plot_ms

peaks = parser_mzxml(file)
peaks1 = [p for p in peaks if p.getMSLevel()==1]
peaks2 = [p for p in peaks if p.getMSLevel()==2]
precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
del(peaks)

exmz = 464.312
exrt = 945.448
dda = dda_result[np.abs( dda_result['precursor_mz'] - exmz ) < 0.01]
dda = dda[np.abs( dda['precursor_rt'] - exrt ) < 15]
exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)

ms2 = get_ms2(peaks2, precursors, exmz, exrt)

decoy = [284.261, 464.287, 466.32, 265.251, 152.996, 388.153]
# precursor eic
plt.plot(exeic[0], exeic[1]) 

# DDA MS2
plot_ms(dda['mz'], dda['intensity']) 

# SWATH MS2
plot_ms(ms2[0], ms2[1])

# Frag MS2
plot_ms(ms2[0], ms2[1], dda['mz'])

# Decoy MS2
plot_ms(ms2[0], ms2[1], None, decoy)

# Frag eics
plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, label='precursor') 
for fragmz in dda['mz']:
    frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=30)
    plt.plot(frageic[0], frageic[1], color='red')
    plt.legend()

# decoy eics
plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, label='precursor') 
for decoymz in decoy:
    decoyeic = fragment_eic(peaks2, precursors, exmz, exrt, decoymz, rtlength=30)
    plt.plot(decoyeic[0], decoyeic[1], color='blue')

#################################################################################
# discussion 2

exmz = 501.276
exrt = 843.297
dda = dda_result[np.abs( dda_result['precursor_mz'] - exmz ) < 0.01]
dda = dda[np.abs( dda['precursor_rt'] - exrt ) < 15]
exeic = extract_eic(peaks1, exmz, exrt, rtlength=30)

ms2 = get_ms2(peaks2, precursors, exmz, exrt)

# precursor eic
plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.9) 

# SWATH MS2
plot_ms(ms2[0], ms2[1])

# all eics
plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.9, label='precursor') 
for fragmz in ms2[0]:
    frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=30)
    plt.plot(frageic[0], frageic[1], color='red', alpha = 0.5)
    plt.legend()

# select eics
dia = dia_result[np.abs( dia_result['precursor_mz'] - exmz ) < 0.01]
dia = dia[np.abs( dia['precursor_rt'] - exrt ) < 15]
plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.9, label='precursor') 
for fragmz in ms2[0]:
    frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=30)
    plt.plot(frageic[0], frageic[1], color='blue', alpha = 0.5)
for fragmz in dia['mz']:
    frageic = fragment_eic(peaks2, precursors, exmz, exrt, fragmz, rtlength=30)
    plt.plot(frageic[0], frageic[1], color='red', alpha = 0.9)
    plt.legend()

# select ions
plt.vlines(dda['mz'], 0, dda['intensity']/ np.max(dda['intensity']), color='red', alpha=0.8, label='DDA')
plt.vlines(dia['mz'], 0, -dia['intensity']/ np.max(dia['intensity']), color='blue', alpha=0.8, label='DeepDIA')
plt.axhline(0, color='black')
plt.xlabel('m/z')
plt.ylabel('Abundance')
plt.legend()