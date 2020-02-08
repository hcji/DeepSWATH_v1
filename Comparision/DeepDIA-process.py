# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:55:44 2020

@author: hcji
"""

import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process

# MetDIA data process
f = 'Comparision/MetDIA_data/data/30STD_mix 330ppb-1.mzML'
features = pd.read_csv('Comparision/MetDIA_data/results/xcms_ms1_feature.csv')
features = features[['mz', 'rt', 'maxo']]
features.columns = ['mz', 'rt', 'intensity']
dia_result = DeepDIA_process(f, features)
dia_result.to_csv('Comparision/MetDIA_data/results/30STD_mix_330ppb_1_DeepDIA.csv')

# MetaboDIA data process
f = 'Comparision/MetaboDIA_data/data/PH697097_pos_SWATH.mzML'
features = pd.read_csv('Comparision/MetaboDIA_data/results/xcms_ms1_feature_pos.csv')
features = features[['mz', 'rt', 'maxo']]
features.columns = ['mz', 'rt', 'intensity']
dia_result = DeepDIA_process(f, features)
dia_result.to_csv('Comparision/MetaboDIA_data/results/PH697097_pos_DeepDIA.csv')

f = 'Comparision/MetaboDIA_data/data/PH697097_neg_SWATH.mzML'
features = pd.read_csv('Comparision/MetaboDIA_data/results/xcms_ms1_feature_neg.csv')
features = features[['mz', 'rt', 'maxo']]
features.columns = ['mz', 'rt', 'intensity']
dia_result = DeepDIA_process(f, features)
dia_result.to_csv('Comparision/MetaboDIA_data/results/PH697097_neg_DeepDIA.csv')
