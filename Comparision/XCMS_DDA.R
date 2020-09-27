# xcms version 3.8.1
library(xcms)

get_DDA_MS <- function(file, snthresh = 10, prefilter = c(1, 100), ppm = 100, peakwidth = c(5, 60), noise = 100){
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
    frag_mz = frag_mz[frag_abund >= noise]
    frag_abund = frag_abund[frag_abund >= noise]
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
write.csv(metabodia_pos_fet, 'Comparision/MetaboDIA_Data/results/PH697097_pos_IDA.csv')
write.csv(metabodia_neg_fet, 'Comparision/MetaboDIA_Data/results/PH697097_neg_IDA.csv')

# MetaboKit data
metabodia_pos_file <- 'D:/MetaboKit_data/IDA/POS/SWATH-POS-MODE-PathwayMetabolites_PosIDA-1-Metabolite-Mix-_PosIDA-1.mzXML'
metabodia_neg_file <- 'D:/MetaboKit_data/IDA/NEG/SWATH-NEG-MODE-PathwayMetabolites_NEGIDA-2-Metabolite-Mix-_NEGIDA-1.mzML'
metabodia_pos_fet <- get_DDA_MS(metabodia_pos_file, snthresh = 3, prefilter = c(1, 50), noise = 30)
metabodia_neg_fet <- get_DDA_MS(metabodia_neg_file, snthresh = 3, prefilter = c(1, 50), noise = 30)

compounds <- read.csv('D:/MetaboKit_data/IDA/compounds.csv')
compounds.mzpos <- round(compounds$NEUTRAL.MONOISOTOPIC.MOLECULAR.MASS + 1.0034, 1)
compounds.mzneg <- round(compounds$NEUTRAL.MONOISOTOPIC.MOLECULAR.MASS - 1.0034, 1)

metabodia_pos_fet <- metabodia_pos_fet[round(metabodia_pos_fet[,'mz'], 1) %in% compounds.mzpos, ]
metabodia_neg_fet <- metabodia_neg_fet[round(metabodia_neg_fet[,'mz'], 1) %in% compounds.mzneg, ]

write.csv(metabodia_pos_fet, 'Comparision/MetaboKit_Data/results/MetaboKit_pos_IDA.csv')
write.csv(metabodia_neg_fet, 'Comparision/MetaboKit_Data/results/MetaboKit_neg_IDA.csv')
