# xcms version 3.8.1
library(xcms)

getisotarmz <- function(win_file, swath_data){
  exp <- read.delim(win_file)
  tarmz <- rowMeans(exp[,3:4])
  tarmz[1] <- NA
  return(rep(tarmz, length(swath_data)/length(tarmz)))
}

xcms_swath <- function(file){
  cwp <- CentWaveParam(snthresh = 10, noise = 200, ppm = 30,
                       peakwidth = c(5, 60))
  swath_data <- readMSData(file, mode = "onDisk")
  swath_data@featureData@data$isolationWindowTargetMZ <- swath_data@featureData@data$precursorMZ
  swath_data@featureData@data$isolationWindowLowerOffset <- swath_data@featureData@data$precursorMZ - 12.5
  swath_data@featureData@data$isolationWindowUpperOffset <- swath_data@featureData@data$precursorMZ + 12.5
  
  swath_data <- findChromPeaks(swath_data, param = cwp)

  cwp2 <- CentWaveParam(snthresh = 5, noise = 200, ppm = 30,
                       peakwidth = c(5, 60))
  swath_data <- findChromPeaksIsolationWindow(swath_data, param = cwp2)
  swath_spectra <- reconstructChromPeakSpectra(swath_data, minCor = 0.7)
  

  
}

file <- 'Example/CS52684_neg_SWATH.mzXML'
win_file <- 'Example/windows.txt'
