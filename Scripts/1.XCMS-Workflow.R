# xcms version 3.8.1
library(xcms)

xcms_process <- function(dir, replace=TRUE){
  files <- list.files(dir, pattern = '*_IDA.mzXML', full.names = TRUE)
  exist <- list.files(dir, pattern = '*.ms2.csv', full.names = TRUE)
  
  for (f in files){
    o <- stringr::str_replace(f, '.mzXML', '.ms2.csv')
    if (!replace && o %in% exist){
      next
    }
    dda_data <- readMSData(f, mode = "onDisk")
    
    cwp <- CentWaveParam(snthresh = 10, noise = 200, ppm = 30, peakwidth = c(5, 60))
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
                                    ppm = 30, minProp = 0.8, weighted = FALSE,
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
    write.csv(all_ms2, o)
    print(paste(f, 'finished'))
  }
}

xcms_process('D:/MetaboDIA_data/CS')
xcms_process('D:/MetaboDIA_data/PH')
