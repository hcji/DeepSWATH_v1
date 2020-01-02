# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 09:14:11 2019

@author: yn
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bisect import bisect_left, bisect_right
from pyopenms import MSExperiment, MzMLFile, MzXMLFile
from pyopenms import PeakFileOptions, FeatureFinder
from pyopenms import LogType, FeatureMap
from pyopenms import TheoreticalSpectrumGenerator, AASequence, MSSpectrum

def parser_mzml(raw_file):
    exp = MSExperiment()
    MzMLFile().load(raw_file, exp)
    peaks = exp.getSpectra()
    return peaks


def parser_mzxml(raw_file):
    exp = MSExperiment()
    MzXMLFile().load(raw_file, exp)
    peaks = exp.getSpectra()
    return peaks


def extract_eic(peaks, mz, rt, mztol=0.05, rtlength=30):
    rts = []
    abunds = []
    mzrange = [mz - mztol, mz + mztol]
    rtrange = [rt - rtlength, rt + rtlength]
    for p in peaks:
        if (p.getRT() >= rtrange[0]) and (p.getRT() <= rtrange[1]) and (p.getMSLevel()==1):
            rts.append(p.getRT())
            mzs, intensities = p.get_peaks()
            sel_peaks = np.arange(bisect_left(mzs, mzrange[0]), bisect_right(mzs, mzrange[1]))
            # sel_peaks = np.where( np.logical_and (mzs >= mzrange[0], mzs <= mzrange[1]))[0]
            abunds.append(np.sum(intensities[sel_peaks]))
    return rts, abunds


def fragment_eic(peaks, precursors, exmz, exrt, fragmz, mztol=0.05, rtlength=30):
    precursor = precursors[np.argmin(np.abs(precursors - exmz))]
    rts, abunds = [], []
    mzrange = [fragmz - mztol, fragmz + mztol]
    rtrange = [exrt - rtlength, exrt + rtlength]    
    for p in peaks:
        if (p.getRT() >= rtrange[0]) and (p.getMSLevel()==2):
            if p.getRT() > rtrange[1]:
                break
            if p.getPrecursors()[0].getMZ() == precursor:
                rts.append(p.getRT())
                mzs, intensities = p.get_peaks()
                sel_peaks = np.arange(bisect_left(mzs, mzrange[0]), bisect_right(mzs, mzrange[1]))
                # sel_peaks = np.where( np.logical_and (mzs >= mzrange[0], mzs <= mzrange[1]))[0]
                abunds.append(np.sum(intensities[sel_peaks]))                
    return rts, abunds


def extract_features(raw_file, params=None):
 
    # options: only ms1        
    options = PeakFileOptions()
    options.setMSLevels([1])
    fh = MzMLFile()
    fh.setOptions(options)

    # Load data
    exp = MSExperiment()
    fh.load(raw_file, exp)
    exp.updateRanges()
    
    # set featurefinder
    ff = FeatureFinder()
    ff.setLogType(LogType.CMD)
    
    name = "centroided"
    features = FeatureMap()
    seeds = FeatureMap()
    p = FeatureFinder().getParameters(name)
    
    keys = [str(k, 'utf-8') for k in p.keys()]
    if params is not None:
        for i, val in params.items():
            if i in keys:
                p.setValue(i, val)
    
    # run featurefinder
    ff.run(name, exp, features, p, seeds)
    features.setUniqueIds()
    
    rts, mzs, intensities = [], [], []
    for f in features:
        mzs.append(f.getMZ())
        rts.append(f.getRT())
        intensities.append(f.getIntensity())
    output = pd.DataFrame({'mz': mzs, 'rt':rts, 'intensity':intensities})
    return output


def get_ms2(peaks, precursors, exmz, exrt):
    precursor = precursors[np.argmin(np.abs(precursors - exmz))]
    nearest = 10**6
    ms2 = None
    for p in peaks:
        if (abs(exrt - p.getRT()) < nearest) and (p.getMSLevel()==2):
            if p.getPrecursors()[0].getMZ() == precursor:
                ms2 = p.get_peaks()
                nearest = abs(exrt - p.getRT())
    return ms2
