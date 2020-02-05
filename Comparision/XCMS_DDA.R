# xcms version 3.8.1
library(xcms)

get_DDA_MS <- function(file, snthresh = 10, prefilter = c(1, 100), ppm = 100, peakwidth = c(5, 60)){
  dda_data <- readMSData(file, mode = "onDisk")
  
  cwp <- CentWaveParam(snthresh = snthresh, prefilter = prefilter, ppm = ppm, peakwidth = peakwidth)
  dda_data <- findChromPeaks(dda_data, param = cwp)
  dda_spectra <- chromPeakSpectra(dda_data)
  
  peaks <- chromPeaks(dda_data)
  ex_id <- rownames(chromPeaks(dda_data, ppm = 30))
  all_ms2 <- NULL
  for (id in ex_id){
    ex_mz <- peaks[id, 'mz']
    ex_rt <- peaks[id, 'rt']
    ex_int <- peaks[id, 'maxo']
    ex_spectra <- dda_spectra[mcols(dda_spectra)$peak_id == id]
    ex_spectrum <- combineSpectra(ex_spectra, method = consensusSpectrum, mzd = 0,
                                  ppm = ppm, minProp = 0.8, weighted = FALSE,
                                  intensityFun = median, mzFun = median)
    if (length(ex_spectrum) == 0 || length(mz(ex_spectrum)[[1]]) == 0){
      next
    }
    frag_mz = mz(ex_spectrum)[[1]]
    frag_abund = intensity(ex_spectrum)[[1]]
    frag_mz = frag_mz[frag_abund >= 100]
    frag_abund = frag_abund[frag_abund >= 100]
    if (length(frag_mz) == 0){
      next
    }
    ms2 <- cbind(ex_mz, ex_rt, ex_int, frag_mz, frag_abund)
    all_ms2 <- rbind(all_ms2, ms2)
  }
  colnames(all_ms2) <- c('precursor_mz', 'precursor_rt', 'precursor_intensity', 'mz', 'intensity')
  return(all_ms2)
}


# MetaboDIA data
metabodia_pos_file <- 'Comparision/MetaboDIA_Data/data/PH697097_pos_IDA.mzML'
metabodia_neg_file <- 'Comparision/MetaboDIA_Data/data/PH697097_neg_IDA.mzML'
metabodia_pos_fet <- get_DDA_MS(metabodia_pos_file, prefilter = c(1, 500))
metabodia_neg_fet <- get_DDA_MS(metabodia_neg_file, prefilter = c(1, 500))
write.csv(metabodia_pos_fet, 'Comparision/MetaboDIA_Data/results/DDA_pos.csv')
write.csv(metabodia_neg_fet, 'Comparision/MetaboDIA_Data/results/DDA_neg.csv')

# MS-DIAL data
msdial_pos_file <- 'Comparision/MSDIAL_Data/data/Posi_Ida_QC_1_1.mzML'
msdial_neg_file <- 'Comparision/MSDIAL_Data/data/Nega_Ida_QC_1_1.mzML'
msdial_pos_fet <- get_DDA_MS(msdial_pos_file, prefilter = c(1, 500))
msdial_neg_fet <- get_DDA_MS(msdial_neg_file, prefilter = c(1, 500))
write.csv(msdial_pos_fet, 'Comparision/MSDIAL_Data/results/DDA_pos.csv')
write.csv(msdial_neg_fet, 'Comparision/MSDIAL_Data/results/DDA_neg.csv')
