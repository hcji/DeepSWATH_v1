# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 08:04:29 2020

@author: hcji
"""

import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process


file = 'Example/CS52684_neg_SWATH.mzXML'
features = pd.read_csv('Example/CS52684_neg_SWATH.features.csv')
    
dia_result = DeepDIA_process(file, features)
dia_result.to_csv('Example/CS52684_neg_SWATH.ms2.csv')
dda_result = pd.read_csv('Example/CS52684_neg_IDA.ms2.csv')