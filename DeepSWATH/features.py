# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 10:11:54 2020

@author: jihon
"""

import sys
import numpy as np
import pandas as pd
from scipy import signal
from tqdm import tqdm
from ICExtract._ICExtract import DIAMS, DIA_PIC, MS2_PIC

def ms2_features(mzfile, tol = 0.025, min_int = 300, missing_scan = 1, min_len = 10):
    # mzfile = "D:/data/DIA-CETSA/K562-STP-48degree-20uM-DIA-1h-Fusion-OT-0621.mzML"
    mzdata = mzfile.encode(sys.getfilesystemencoding())
    diams = DIAMS()
    diams.loadfile(mzdata)
    swaths = diams.Swaths
    
    print('Extracting ion trace')
    # ms1_pic = DIA_PIC(diams, tol, min_int, missing_scan, min_len)
    ms2_pic = []
    for i in tqdm(np.arange(1, len(swaths))):
        p = MS2_PIC(diams, int(i), 0, 10**10, tol, min_int, missing_scan, int(min_len/2))
        ms2_pic.append(p)
    
    print('Detecting peaks')
    '''
    ms1_feat_mz,  ms1_feat_rt, ms1_feat_int = np.array([]), np.array([]), np.array([])
    for p in ms1_pic:
        pk = get_peaks(p)
        ms1_feat_rt = np.concatenate((ms1_feat_rt, pk[0]))
        ms1_feat_mz = np.concatenate((ms1_feat_mz, pk[1]))
        ms1_feat_int = np.concatenate((ms1_feat_int, pk[2]))
    ms1_feat = pd.DataFrame({'mz': ms1_feat_mz, 'rt': ms1_feat_rt, 'intensity': ms1_feat_int})
    '''
    ms2_all = []
    for i in tqdm(range(len(ms2_pic))):
        ms2_feat_mz,  ms2_feat_rt, ms2_feat_int = np.array([]), np.array([]), np.array([])
        for p in ms2_pic[i]:
            pk = get_peaks(p)
            ms2_feat_rt = np.concatenate((ms2_feat_rt, pk[0]))
            ms2_feat_mz = np.concatenate((ms2_feat_mz, pk[1]))
            ms2_feat_int = np.concatenate((ms2_feat_int, pk[2]))
        ms2_feat = pd.DataFrame({'mz': ms2_feat_mz, 'rt': ms2_feat_rt, 'intensity': ms2_feat_int})
        ms2_all.append(ms2_feat)
        
    return swaths, ms2_all


def get_peaks(pic, lamb=0.0, width=1):
    sig = pic[:,2]
    sig = signal.qspline1d(sig, lamb)
    peaks = signal.find_peaks(sig, width=width)[0]
    rt = pic[peaks, 0]
    mz = pic[peaks, 1]
    intensity = pic[peaks, 2]
    return rt, mz, intensity
    
    
