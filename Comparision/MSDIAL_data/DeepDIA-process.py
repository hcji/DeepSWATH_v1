# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 10:55:44 2020

@author: hcji
"""

import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process

f = 'E:/project/MSDIAL_data/20140809_MSDIAL_DemoFiles_Swath (wiff.wiffscan)/Nega_Swath_QC_1_1.mzML'
features = pd.read_csv('E:/project/MSDIAL_data/20140809_MSDIAL_DemoFiles_Swath (wiff.wiffscan)/Nega_Swath_QC_1_1.feature.csv')
features = features[['mz', 'rt', 'maxo']]
features.columns = ['mz', 'rt', 'intensity']
dia_result = DeepDIA_process(f, features)
dia_result.to_csv('Comparision/MSDIAL_data/results/DeepDIA_results/Nega_Swath_QC_1_1.ms2.csv')
