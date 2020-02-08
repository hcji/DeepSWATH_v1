# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:25:49 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from libmetgem import msp

# MSDIAL 4.12.
def parser_msdial(file, output):
    # file = 'E:/project/MetaboDIA_data/MS-DIAL Result/PH697338_pos_SWATH.mgf'
    parser = msp.read(file)
    exmz, exrt, exint, mz, intensity = [], [], [], [], []
    for (data, spectrum) in parser:
        if spectrum.shape[0] < 3:
            continue
        keep = np.where(spectrum[:,1] >= 0.02 * np.max(spectrum[:,1]))[0]
        exmz += [data['precursormz']] * len(keep)
        exrt += [data['retentiontime'] * 60] * len(keep)
        exint += [data['intensity']] * len(keep)
        mz += list(spectrum[keep,0])
        intensity += list(spectrum[keep,1])
    res = pd.DataFrame({'precursor_mz':exmz, 'precursor_rt':exrt, 'precursor_intensity':exint,
                           'mz': mz, 'intensity': intensity})
    res.to_csv(output, index=False)
    return res

if __name__ == '__main__':
    
    # MetDIA DIA data
    metdia_dia = 'Comparision/MetDIA_Data/results/30STD_mix 330ppb-1.mgf'
    parser_msdial(metdia_dia, 'Comparision/MetDIA_Data/results/30STD_mix_330ppb_1_MSDIAL.csv')
    
    # MetaboDIA DIA data
    metabodia_dia_pos = 'Comparision/MetaboDIA_Data/results/PH697097_pos_SWATH.mgf'
    metabodia_dia_neg = 'Comparision/MetaboDIA_Data/results/PH697097_neg_SWATH.mgf'
    parser_msdial(metabodia_dia_pos, 'Comparision/MetaboDIA_Data/results/PH697097_pos_MSDIAL.csv')
    parser_msdial(metabodia_dia_neg, 'Comparision/MetaboDIA_Data/results/PH697097_neg_MSDIAL.csv')    
    