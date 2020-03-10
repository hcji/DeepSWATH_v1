# -*- coding: utf-8 -*-
"""
Created on Sat Jan 11 08:10:25 2020

@author: hcji
"""

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm


def count_features(ms2):
    features = pd.read_csv(ms2)
    features = features[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    return len(features)


def get_str_for_venn(ms2):
    features = pd.read_csv(ms2)
    features = features[features['intensity'] > 100]
    v = [str(round(features.loc[r,'precursor_mz'],2)) + ' ' +
         str(round(features.loc[r,'precursor_rt'],-1)) + ' ' +
         str(round(features.loc[r,'mz'],1)) for r in features.index] 
    return set(v)


def split_precursor(mslist):
    f1 = list(set(zip(mslist['precursor_mz'], mslist['precursor_rt'])))
    output = []
    for i in f1:
        ms = mslist[mslist['precursor_mz'] == i[0]]
        ms = ms[ms['precursor_rt'] == i[1]]
        output.append(ms)
    return output


def DB_DIA_compare(db, f_deepdia, f_msdial, mztol=0.05):
    db_res = pd.read_csv(db)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    
    features = db_res[['Name', 'PrecursorMz']]
    features = features.drop_duplicates()
    
    output = pd.DataFrame(columns=['Name', 'precursor_mz', 'precursor_intensity', 'DeepDIA_corr', 'MSDIAL_corr'])
    for i in tqdm(range(features.shape[0])):
        f = features.iloc[i,:]
        standard = db_res[db_res['PrecursorMz'] == f['PrecursorMz']]
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['PrecursorMz']) < mztol]
        precursor_intensity = np.mean(deepdia['precursor_intensity'])
        
        deepdia = split_precursor(deepdia)
        msdial = split_precursor(msdial)
        
        mzs = standard['ProductMz']
        std_int = np.array(standard['LibraryIntensity'])
        
        deepdia_corrs, msdial_corrs = [], []
        for j in range(len(deepdia)):
            deepdia_int = []
            for mz in mzs:
                deepdia_int.append(np.sum(deepdia[j]['intensity'][ np.abs(deepdia[j]['mz'] - mz) < 0.1 ]))
            deepdia_corrs.append(np.nanmax([0, pearsonr(std_int, deepdia_int)[0]]))
        
        for j in range(len(msdial)):
            msdial_int = []
            for mz in mzs:
                msdial_int.append(np.sum(msdial[j]['intensity'][ np.abs(msdial[j]['mz'] - mz) < 0.1 ]))
            msdial_corrs.append(np.nanmax([0, pearsonr(std_int, msdial_int)[0]]))

        deepdia_corr = max(deepdia_corrs)
        msdial_corr = max(msdial_corrs)
        output.loc[i] = [f['Name'], f['PrecursorMz'], precursor_intensity, deepdia_corr, msdial_corr]
    return output
        

def DDA_DIA_compare(f_dda, f_deepdia, f_msdial, mztol=0.05, rttol=40):
    dda_res = pd.read_csv(f_dda)
    deepdia_res = pd.read_csv(f_deepdia)
    msdial_res = pd.read_csv(f_msdial)
    # for fair comparesion, exclude low intensity peaks in msdial
    msdial_res = msdial_res[msdial_res['intensity'] > 150]
    
    features = dda_res[['precursor_mz', 'precursor_rt', 'precursor_intensity']]
    features = features.drop_duplicates()
    
    output = pd.DataFrame(columns=['precursor_mz', 'precursor_rt', 'precursor_intensity', 'DeepDIA_corr', 'MSDIAL_corr'])
    for i in tqdm(range(features.shape[0])):
        f = features.iloc[i,:]
        dda = dda_res[dda_res['precursor_mz'] == f['precursor_mz']]
        dda = dda[dda['precursor_rt'] == f['precursor_rt']]
        
        deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - f['precursor_mz']) < mztol ]
        deepdia = deepdia[ np.abs(deepdia['precursor_rt'] - f['precursor_rt']) < rttol ]
        msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - f['precursor_mz']) < mztol ]
        msdial = msdial[ np.abs(msdial['precursor_rt'] - f['precursor_rt']) < rttol ]

        deepdia = split_precursor(deepdia)
        msdial = split_precursor(msdial)

        if dda.shape[0] <= 3:
            continue
        
        deepdia_corrs, msdial_corrs = [], []
        for j in range(len(deepdia)):
            mzs = list(set( list(np.round(dda['mz'], 1)) + list(np.round(deepdia[j]['mz'], 1)) ))
            dda_int, deepdia_int = [], []
            for mz in mzs:
                dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.15 ]))
                deepdia_int.append(np.sum(deepdia[j]['intensity'][ np.abs(deepdia[j]['mz'] - mz) < 0.15 ]))
            deepdia_corrs.append(np.nanmax([0, pearsonr(dda_int, deepdia_int)[0]]))            
        
        for j in range(len(msdial)):
            mzs = list(set( list(np.round(dda['mz'], 1)) + list(np.round(msdial[j]['mz'], 1)) ))
            dda_int, msdial_int = [], []
            for mz in mzs:
                dda_int.append(np.sum(dda['intensity'][ np.abs(dda['mz'] - mz) < 0.15 ]))
                msdial_int.append(np.sum(msdial[j]['intensity'][ np.abs(msdial[j]['mz'] - mz) < 0.15 ]))
            msdial_corrs.append(np.nanmax([0, pearsonr(dda_int, msdial_int)[0]]))            
        
        if len(deepdia_corrs) == 0:
            deepdia_corr = 0
        else:
            deepdia_corr = max(deepdia_corrs)
            
        if len(msdial_corrs) == 0:
            msdial_corr = 0
        else:
            msdial_corr = max(msdial_corrs)
            
        if (deepdia_corr < 0.2) or (msdial_corr < 0.2):
            continue
            
        output.loc[i] = [f['precursor_mz'], f['precursor_rt'], f['precursor_intensity'], deepdia_corr, msdial_corr]
    return output

        
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    # MetDIA Comparision
    db = 'Comparision/MetDIA_data/results/30STD_Database.csv'
    f_deepdia = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_DeepDIA.csv'
    f_msdial = 'Comparision/MetDIA_data/results/30STD_mix_330ppb_1_MSDIAL.csv'
    metdia = DB_DIA_compare(db, f_deepdia, f_msdial, mztol=0.05)
    metdia.to_csv('Comparision/MetDIA_Data/results/comparision.csv')
    
    # MetaboDIA Comparision
    p_dda = 'Comparision/MetaboDIA_data/results/PH697097_pos_IDA.csv'
    n_dda = 'Comparision/MetaboDIA_data/results/PH697097_neg_IDA.csv'
    p_deepdia = 'Comparision/MetaboDIA_data/results/PH697097_pos_DeepDIA.csv'
    n_deepdia = 'Comparision/MetaboDIA_data/results/PH697097_neg_DeepDIA.csv'
    p_msdial = 'Comparision/MetaboDIA_data/results/PH697097_pos_MSDIAL.csv'
    n_msdial = 'Comparision/MetaboDIA_data/results/PH697097_neg_MSDIAL.csv'
    
    p_metabodia = DDA_DIA_compare(p_dda, p_deepdia, p_msdial)
    n_metabodia = DDA_DIA_compare(n_dda, n_deepdia, n_msdial)
    p_metabodia.to_csv('Comparision/MetaboDIA_Data/results/comparision_pos.csv')
    n_metabodia.to_csv('Comparision/MetaboDIA_Data/results/comparision_neg.csv')

    # Correlation Violin Plot
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
    axes[0,0].violinplot([list(metdia['DeepDIA_corr'])], [1], showmeans=False, showmedians=True)
    axes[0,0].violinplot([list(metdia['MSDIAL_corr'])], [2], showmeans=False, showmedians=True)
    axes[0,0].set_xticks(range(4))
    axes[0,0].set_xticklabels(['', 'DeepMetDIA', 'MSDIAL', ''])
    axes[0,0].set_ylabel('Correlation')
      
    axes[0,1].violinplot([list(p_metabodia['DeepDIA_corr']), list(n_metabodia['DeepDIA_corr'])], [1,5], showmeans=False, showmedians=True)
    axes[0,1].violinplot([list(p_metabodia['MSDIAL_corr']), list(n_metabodia['MSDIAL_corr'])], [3,7], showmeans=False, showmedians=True)
    axes[0,1].set_xticks(range(8))
    axes[0,1].set_xticklabels(['', 'DeepSWATH', '\nPositive', 'MSDIAL', '', 'DeepSWATH', '\nNegative', 'MSDIAL'])
    axes[0,1].set_ylabel('Correlation')
    
    exp_mz = 165.054
    exp_rt = 123.067
    mztol = 0.05
    rttol = 30
    dda_res = pd.read_csv(p_dda)
    deepdia_res = pd.read_csv(p_deepdia)
    msdial_res = pd.read_csv(p_msdial) 
    
    dda = dda_res[ np.abs(dda_res['precursor_mz'] - exp_mz) < mztol ]
    dda = dda[ np.abs(dda['precursor_rt'] - exp_rt) < rttol ]
    deepdia = deepdia_res[ np.abs(deepdia_res['precursor_mz'] - exp_mz) < mztol ]
    deepdia = deepdia[ np.abs(deepdia['precursor_rt'] - exp_rt) < rttol ]
    msdial = msdial_res[ np.abs(msdial_res['precursor_mz'] - exp_mz) < mztol ]
    msdial = msdial[ np.abs(msdial['precursor_rt'] - exp_rt) < rttol ]
    # msdial = msdial[msdial['intensity'] > 150]
    
    axes[1,0].vlines(dda['mz'], 0, dda['intensity'] / np.max(dda['intensity']), color='red', alpha=0.7, label='DDA')
    axes[1,0].vlines(deepdia['mz'], 0, -deepdia['intensity'] / np.max(deepdia['intensity']), color='blue', alpha=0.7, label='DeepSWATH')
    axes[1,0].text(400, -0.8, 'MS/MS \n precursor:'+str(exp_mz))
    axes[1,0].axhline(0, color='black')
    axes[1,0].set_xlabel('m/z')
    axes[1,0].set_ylabel('Abundance')
    axes[1,0].legend()
    
    axes[1,1].vlines(dda['mz'], 0, dda['intensity']/ np.max(dda['intensity']), color='red', alpha=0.8, label='DDA')
    axes[1,1].vlines(msdial['mz'], 0, -msdial['intensity']/ np.max(msdial['intensity']), color='orange', alpha=0.8, label='MS-DIAL')
    axes[1,1].text(400, -0.8, 'MS/MS \nprecursor:'+str(exp_mz))
    axes[1,1].axhline(0, color='black')
    axes[1,1].set_xlabel('m/z')
    axes[1,1].set_ylabel('Abundance')
    axes[1,1].legend()
    
    plt.savefig('Figure/DDA_DIA_compare.png')
    

    # coelution handling
    # example 1
    from DeepDIA.utils import parser_mzml, extract_eic, fragment_eic, get_ms2
    from DeepDIA.plot import plot_ms
    
    file = 'Comparision/MetaboDIA_Data/data/PH697097_pos_SWATH.mzML'
    peaks = parser_mzml(file)
    peaks1 = [p for p in peaks if p.getMSLevel()==1]
    peaks2 = [p for p in peaks if p.getMSLevel()==2]
    precursors = np.unique([p.getPrecursors()[0].getMZ() for p in peaks2])
    del(peaks)

    exmz1 = 266.094
    exrt1 = 391.224
    exmz2 = 266.122
    exrt2 = 398.333
    
    ms2 = get_ms2(peaks2, precursors, exmz1, exrt1)
    exeic = extract_eic(peaks1, exmz1, exrt1, rtlength=30)
     
    dia1 = deepdia_res[np.round(deepdia_res['precursor_mz'],3) == exmz1]
    dia1 = dia1[np.round(dia1['precursor_rt'],3) == exrt1]
    dia2 = deepdia_res[np.round(deepdia_res['precursor_mz'],3) == exmz2]
    dia2 = dia2[np.round(dia2['precursor_rt'],3) == exrt2]
    
    fragmz1 = dia1[dia1['intensity'] > 0.1* max(dia1['intensity'])]['mz']
    fragmz2 = dia2[dia2['intensity'] > 0.1* max(dia2['intensity'])]['mz']
    
    plt.plot(exeic[0], exeic[1], color = 'black', linewidth=2, alpha = 0.5, label='precursor') 
    for fragmz in list(fragmz1)+list(fragmz2)+list(ms2[0]):
        frageic = fragment_eic(peaks2, precursors, exmz1, exrt1, fragmz, mztol=0.01, rtlength=30)
        plt.plot(frageic[0], frageic[1], color='blue', linewidth=2, alpha = 0.8)
        plt.xlim(380, 420)
        plt.legend(['precursor', 'candidates'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')
    
    plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, linewidth=2, label='precursor') 
    for fragmz in fragmz1:
        frageic = fragment_eic(peaks2, precursors, exmz1, exrt1, fragmz, mztol=0.01, rtlength=30)
        plt.plot(frageic[0], frageic[1], color='red', linewidth=2)
        plt.xlim(380, 420)
        # plt.xlim(385,395)
        # plt.ylim(0, 8000)
        plt.legend(['precursor', 'fragments'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')
    
    plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, linewidth=2, label='precursor') 
    for fragmz in fragmz2:
        frageic = fragment_eic(peaks2, precursors, exmz1, exrt1, fragmz, mztol=0.01, rtlength=30)
        plt.plot(frageic[0], frageic[1], color='green', linewidth=2)
        plt.xlim(380, 420)
        plt.legend(['precursor', 'fragments'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')
    
    
    # example 2
    exmz1 = 136.075
    exrt1 = 101.729
    
    ms2 = get_ms2(peaks2, precursors, exmz1, exrt1)
    exeic = extract_eic(peaks1, exmz1, exrt1, mztol=0.05, rtlength=30)
     
    dia1 = deepdia_res[np.round(deepdia_res['precursor_mz'],3) == exmz1]
    dia1 = dia1[np.round(dia1['precursor_rt'],3) == exrt1]
    
    fragmz1 = dia1[dia1['intensity'] > 0.1* max(dia1['intensity'])]['mz']
    
    plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, label='precursor') 
    for fragmz in ms2[0]:
        frageic = fragment_eic(peaks2, precursors, exmz1, exrt1, fragmz, mztol=0.01, rtlength=30)
        plt.plot(frageic[0], frageic[1], color='blue', alpha = 0.8)
        plt.legend(['precursor', 'candidates'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')
    
    fragmat = []
    plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, label='precursor') 
    for fragmz in fragmz1:
        frageic = fragment_eic(peaks2, precursors, exmz1, exrt1, fragmz, mztol=0.01, rtlength=30)
        fragmat.append(frageic[1])
        plt.plot(frageic[0], frageic[1], linewidth=2)
        plt.legend(['precursor', 'fragments'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')
    
    # combine with mcr
    from mcr import MCR
    from mcr.efa import efa
    
    M = np.array(fragmat).T
    C_ = efa(M, ncomp=2, plot=False)
    mcr = MCR(maxiter=1000, chkpnt=10, tol=1e-10)
    mcr.fit(M, C=C_, debug=True)
    B = mcr.B
    plt.plot(exeic[0], exeic[1], color = 'black', alpha = 0.5, label='precursor', linewidth=2) 
    for i in range(M.shape[1]):
        plt.plot(exeic[0], M[:,i] * B[:,1]/ max(B[:,1]), color = 'red', linewidth=2)
        plt.legend(['precursor', 'fragments'])
    plt.xlabel('RT')
    plt.ylabel('Abundance')    
    