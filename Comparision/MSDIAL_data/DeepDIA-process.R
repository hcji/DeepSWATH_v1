# xcms version 3.8.1
library(xcms)

get_ms1_features <- function(file, snthresh = 10, prefilter = c(1, 100), noise = 100, ppm = 30, peakwidth = c(5, 60)){
  data <- readMSData(file, mode = "onDisk") 
  cwp <- CentWaveParam(snthresh = snthresh, noise = noise, prefilter = prefilter, ppm = ppm, peakwidth = peakwidth)
  data <- findChromPeaks(data, param = cwp)
  peaks <- chromPeaks(data)
  return(peaks)
}

file <- 'E:/project/MSDIAL_data/20140809_MSDIAL_DemoFiles_Swath (wiff.wiffscan)/Nega_Swath_QC_1_1.mzML'
features = get_ms1_features(file, snthresh = 5, prefilter = c(1, 1000))
write.csv(features, stringr::str_replace(file, '.mzML', '.feature.csv'))