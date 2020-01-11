# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 20:57:03 2020

@author: hcji
"""

import os
import pandas as pd
from DeepDIA.DeepDIA import DeepDIA_process


Dir = 'D:/MetaboDIA_data/PH'
files = os.listdir(Dir)
files = [f for f in files if '_SWATH.mzXML' in f]

def process(f):
    f = Dir + '/' + f
    ft, fo = f, f
    dda_result = pd.read_csv(ft.replace('.mzXML', '.feature.csv'))
    features = dda_result[['mz', 'rt', 'maxo']]
    features.columns = ['precursor_mz', 'precursor_rt', 'precursor_intensity']
    features = features.drop_duplicates()
    features.columns = ['mz', 'rt', 'intensity']
    dia_result = DeepDIA_process(f, features)
    dia_result.to_csv(fo.replace('.mzXML', '.ms2.csv'))

if __name__ == '__main__':
    
    from joblib import Parallel, delayed
    Parallel(n_jobs=6, verbose=5)(delayed(process)(f) for f in files)
    