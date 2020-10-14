# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 09:14:11 2019

@author: yn
"""

import pymzml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bisect import bisect_left, bisect_right

def extract_eic(peaks, exmz, exrt, mztol=0.1, rtlength=30):
    rts, abunds = [], []
    mzrange = [exmz - mztol, exmz + mztol]
    rtrange = [exrt - rtlength, exrt + rtlength]
    for p in peaks:
        if p.scan_time[0] > rtrange[1]:
            break
        if (p.scan_time[0] >= rtrange[0]) and (p.scan_time[0] <= rtrange[1]) and (p.ms_level==1):
           rts.append(p.scan_time[0])
           mzs, intensities = p.centroidedPeaks[:,0], p.centroidedPeaks[:,1]
           sel_peaks = np.arange(bisect_left(mzs, mzrange[0]), bisect_right(mzs, mzrange[1]))
           abunds.append(np.sum(intensities[sel_peaks]))
    # plt.plot(rts, abunds)
    return rts, abunds


def fragment_eic(peaks, precursors, exmz, exrt, fragmz, mztol=0.05, rtlength=30):
    precursor = precursors[np.argmin(np.abs(precursors - exmz))]
    rts, abunds = [], []
    mzrange = [fragmz - mztol, fragmz + mztol]
    rtrange = [exrt - rtlength, exrt + rtlength]    
    for p in peaks:
        if p.scan_time[0] > rtrange[1]:
            break
        if (p.scan_time[0] >= rtrange[0]) and (p.scan_time[0] <= rtrange[1]) and (p.ms_level > 1):
            if abs(p.selected_precursors[0]['mz'] - precursor) < 1:
                rts.append(p.scan_time[0])
                mzs, intensities = p.centroidedPeaks[:,0], p.centroidedPeaks[:,1]
                sel_peaks = np.arange(bisect_left(mzs, mzrange[0]), bisect_right(mzs, mzrange[1]))
                abunds.append(np.sum(intensities[sel_peaks])) 
    # plt.plot(rts, abunds)
    return rts, abunds


def get_ms2(peaks, precursors, exmz, exrt):
    precursor = precursors[np.argmin(np.abs(precursors - exmz))]
    nearest = 10**6
    ms2 = None
    for p in peaks:
        if (abs(exrt - p.scan_time[0]) < nearest) and (p.ms_level==2):
            if p.selected_precursors[0]['mz'] == precursor:
                ms2 = (p.centroidedPeaks[:,0], p.centroidedPeaks[:,1])
                nearest = abs(exrt - p.scan_time[0])
    return ms2
  

def dda_dia_compare(dda_result, dia_result, mztol=0.01, rttol=5):
    features = dda_result[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    output = pd.DataFrame(columns=['precursor_mz', 'precursor_rt', 'precursor_intensity', 'mz', 'dda_intensity', 'dia_intensity'])
    for i in range(features.shape[0]):
        f = features.iloc[i,:]
        dda = dda_result[ np.abs(dda_result['precursor_mz'] - f['precursor_mz']) < mztol ]
        dda = dda[ np.abs(dda['precursor_rt'] - f['precursor_rt']) < rttol ]
        dia = dia_result[ np.abs(dia_result['precursor_mz'] - f['precursor_mz']) < mztol ]
        dia = dia[ np.abs(dia['precursor_rt'] - f['precursor_rt']) < rttol ]
        mzs = list(set( list(np.round(dda['mz'], 2)) + list( np.round(dia['mz'], 2)) ))
        dda_int, dia_int = [], []
        for mz in mzs:
            dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.01 ]))
            dia_int.append(np.sum(dia['intensity'][ np.abs(dia['mz'] - mz) < 0.01 ]))
        resi = pd.DataFrame({'precursor_mz': f['precursor_mz'], 'precursor_rt': f['precursor_rt'], 'precursor_intensity': f['precursor_intensity'],
                             'mz': mzs, 'dda_intensity': dda_int, 'dia_intensity': dia_int})
        output = output.append(resi)
    return output
        
    

