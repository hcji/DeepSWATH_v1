# xcms version 3.8.1
library(xcms)

cs_pos <- list.files('D:/MetaboDIA_data/CS', pattern = '*_pos_IDA.mzXML')
for (f in cs_pos){
  f <- paste('D:/MetaboDIA_data/CS/', f, sep='')
  o <- stringr::str_replace(f, '.mzXML', '.ms2.csv')
  dda_data <- readMSData(f, mode = "onDisk")
  
  cwp <- CentWaveParam(snthresh = 5, noise = 50, ppm = 30, peakwidth = c(5, 60))
  dda_data <- findChromPeaks(dda_data, param = cwp)
  dda_spectra <- chromPeakSpectra(dda_data)
  
  peaks <- chromPeaks(dda_data)
  ex_id <- rownames(chromPeaks(dda_data, ppm = 30))
  all_ms2 <- NULL
  for (id in ex_id){
    ex_mz <- peaks[id, 'mz']
    ex_rt <- peaks[id, 'rt']
    ex_spectra <- dda_spectra[mcols(dda_spectra)$peak_id == id]
    ex_spectrum <- combineSpectra(ex_spectra, method = consensusSpectrum, mzd = 0,
                                  ppm = 30, minProp = 0.8, weighted = FALSE,
                                  intensityFun = median, mzFun = median)
    if (length(ex_spectrum) == 0 || length(mz(ex_spectrum)[[1]]) == 0){
      next
    }
    ms2 <- cbind(ex_mz, ex_rt, mz(ex_spectrum)[[1]], intensity(ex_spectrum)[[1]])
    all_ms2 <- rbind(all_ms2, ms2)
  }
  colnames(all_ms2) <- c('precursor_mz', 'precursor_rt', 'mz', 'intensity')
  write.csv(all_ms2, o)
  print(paste(f, 'finished'))
}