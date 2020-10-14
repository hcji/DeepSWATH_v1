# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 17:33:28 2020

@author: hcji
"""

import numpy as np

def write_msp(ms, name, file,  precursor_type='[M+H]+', ion_mode='Positive'):
    # file = 'MS/unknown.msp'
    msp = open(file, 'w+')
    msp.write('NAME: ' + name + '\n')
    msp.write('PRECURSORMZ: ' + str(np.mean(ms['precursor_mz'])) + '\n')
    msp.write('PRECURSORTYPE: ' + precursor_type + '\n')
    msp.write('IONMODE: ' + ion_mode + '\n')
    msp.write('Num Peaks: ' + str(len(ms)) + '\n')
    for i in ms.index:
        msp.write(str(ms['mz'][i]) + '\t' + str(ms['intensity'][i]) + '\n')
    msp.close()
        