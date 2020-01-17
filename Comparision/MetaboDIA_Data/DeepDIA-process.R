# xcms version 3.8.1
library(xcms)

get_ms1_features <- function(file, snthresh = 10, noise = 200, ppm = 30, peakwidth = c(5, 60)){
  data <- readMSData(file, mode = "onDisk") 
  cwp <- CentWaveParam(snthresh = snthresh, noise = noise, ppm = ppm, peakwidth = peakwidth)
  data <- findChromPeaks(data, param = cwp)
  peaks <- chromPeaks(data)
  return(peaks)
}

files <- list.files('D:\\MetaboDIA_data\\PH', pattern = '*.mzXML', full.names = TRUE)
for (f in files) {
  features = get_ms1_features(f)
  write.csv(features, stringr::str_replace(f, '.mzXML', '.feature.csv'))
}