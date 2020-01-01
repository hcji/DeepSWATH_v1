# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 07:37:03 2020

@author: yn
"""

from DeepDIA.utils import parser_mzxml, extract_eic, fragment_eic, get_ms2
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
numpy2ri.activate()
robjects.r('''source('DeepDIA/xcms.R')''')
get_fingerprint = robjects.globalenv['get_ms1_features']


def DeepDIA_process(file):
    # file = 'D:/MetaboDIA_data/CS/CS53088_neg_SWATH.mzXML'
    features = get_ms1_features(file)