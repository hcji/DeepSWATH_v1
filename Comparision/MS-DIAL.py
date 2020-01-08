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
        exint += [data['intensity'] * spectrum.shape[0]]
        mz += list(spectrum[:,0])
        intensity += list(spectrum[:,1])
    res = pd.DataFrame({'precursor_mz':exmz, 'precursor_rt':exrt, 'precursor_intensity':intensity,
                           'mz': mz, 'intensity': intensity})
    res.to_csv(output, index=False)
    return res

if __name__ == '__main__':
    
    import os
    from tqdm import tqdm
    
    Dir = 'E:/project/MetaboDIA_data/MS-DIAL Result'
    for f in tqdm(os.listdir(Dir)):
        f = 'E:/project/MetaboDIA_data/MS-DIAL Result/' + f
        o = f.replace('.mgf', '.csv')
        parser_msdial(f, o)
