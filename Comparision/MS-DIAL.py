# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:25:49 2020

@author: hcji
"""

import pandas as pd
from libmetgem import msp

# MSDIAL 4.12.
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
    
    # MetaboDIA DDA data
    metabodia_dda_pos = 'Comparision/MetaboDIA_Data/results/PH697097_pos_IDA.mgf'
    metabodia_dda_neg = 'Comparision/MetaboDIA_Data/results/PH697097_neg_IDA.mgf'
    parser_msdial(metabodia_dda_pos, 'Comparision/MetaboDIA_Data/results/PH697097_pos_IDA.csv')
    parser_msdial(metabodia_dda_neg, 'Comparision/MetaboDIA_Data/results/PH697097_neg_IDA.csv')
    
    # MSDIAL DDA data
    msdial_dda_pos = 'Comparision/MSDIAL_Data/results/Posi_Ida_QC_1_1.mgf'
    msdial_dda_neg = 'Comparision/MSDIAL_Data/results/Nega_Ida_QC_1_1.mgf'
    parser_msdial(metabodia_dda_pos, 'Comparision/MSDIAL_Data/results/Posi_Ida_QC_1_1.csv')
    parser_msdial(metabodia_dda_neg, 'Comparision/MSDIAL_Data/results/Nega_Ida_QC_1_1.csv')
    
    # MetDIA DIA data
    metdia_dia = 'Comparision/MetDIA_Data/results/30STD_mix 330ppb-1.mgf'
    parser_msdial(metdia_dia, 'Comparision/MetDIA_Data/results/30STD_mix_330ppb_1_MSDIAL.csv')
    
    # MetaboDIA DIA data
    metabodia_dia_pos = 'Comparision/MetaboDIA_Data/results/PH697097_pos_SWATH.mgf'
    metabodia_dia_neg = 'Comparision/MetaboDIA_Data/results/PH697097_neg_SWATH.mgf'
    parser_msdial(metabodia_dia_pos, 'Comparision/MetaboDIA_Data/results/PH697097_pos_MSDIAL.csv')
    parser_msdial(metabodia_dia_neg, 'Comparision/MetaboDIA_Data/results/PH697097_neg_MSDIAL.csv')    
    
    # MSDIAL DIA data
    msdial_dia_pos = 'Comparision/MSDIAL_Data/results/Posi_Swath_QC_1_1.mgf'
    msdial_dia_neg = 'Comparision/MSDIAL_Data/results/Nega_Swath_QC_1_1.mgf'
    parser_msdial(metabodia_dia_pos, 'Comparision/MSDIAL_Data/results/Posi_MSDIAL_QC_1_1.csv')
    parser_msdial(metabodia_dia_neg, 'Comparision/MSDIAL_Data/results/Nega_MSDIAL_QC_1_1.csv')
    