# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 08:04:29 2020

@author: hcji
"""

import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process

# untargeted
file = 'Example/CS52684_neg_SWATH.mzXML'
features = pd.read_csv('Example/CS52684_neg_SWATH.features.csv', index_col=0)    
dia_result = DeepDIA_process(file, features)
dia_result.to_csv('Example/CS52684_neg_SWATH.ms2.csv')

# targeted
dda_result = pd.read_csv('Example/CS52684_neg_IDA.ms2.csv')
features = dda_result[['precursor_mz', 'precursor_rt']]
features = features.drop_duplicates()
features.columns = ['mz', 'rt']
dia_result = DeepDIA_process(file, features)