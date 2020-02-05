# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:25:49 2020

@author: hcji
"""

import pandas as pd
from libmetgem import msp

def parser_msdial(file, output):
    # file = 'E:/project/MetaboDIA_data/MS-DIAL Result/PH697338_pos_SWATH.mgf'
    parser = msp.read(file)
    exmz, exrt, exint, mz, intensity = [], [], [], [], []
    for (data, spectrum) in parser:
        exmz += [data['precursormz']] * spectrum.shape[0]
        exrt += [data['retentiontime'] * 60] * spectrum.shape[0]
        exint += [data['intensity']] * spectrum.shape[0]
        mz += list(spectrum[:,0])
        intensity += list(spectrum[:,1])
    res = pd.DataFrame({'precursor_mz':exmz, 'precursor_rt':exrt, 'precursor_intensity':exint,
                           'mz': mz, 'intensity': intensity})
    res.to_csv(output, index=False)
    return res

if __name__ == '__main__':
    
    # MetaboDIA data
    metabodia_dda_pos = 'Comparision/MetaboDIA_Data/results/PH697097_pos_IDA.mgf'
    metabodia_dda_neg = 'Comparision/MetaboDIA_Data/results/PH697097_neg_IDA.mgf'
    parser_msdial(metabodia_dda_pos, 'Comparision/MetaboDIA_Data/results/PH697097_pos_IDA.csv')
    parser_msdial(metabodia_dda_neg, 'Comparision/MetaboDIA_Data/results/PH697097_neg_IDA.csv')
    
    # MSDIAL data
    msdial_dda_pos = 'Comparision/MSDIAL_Data/results/Posi_Ida_QC_1_1.mgf'
    msdial_dda_neg = 'Comparision/MSDIAL_Data/results/Nega_Ida_QC_1_1.mgf'
    parser_msdial(metabodia_dda_pos, 'Comparision/MSDIAL_Data/results/Posi_Ida_QC_1_1.csv')
    parser_msdial(metabodia_dda_neg, 'Comparision/MSDIAL_Data/results/Nega_Ida_QC_1_1.csv')
    