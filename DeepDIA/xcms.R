# check xcms
if (!requireNamespace("xcms", quietly = TRUE)){
  if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
  } 
  BiocManager::install("xcms")
}

# xcms version 3.8.1
library(xcms)

get_ms1_features <- function(file, snthresh = 10, noise = 200, ppm = 30, peakwidth = c(5, 60)){
  data <- readMSData(file, mode = "onDisk") 
  cwp <- CentWaveParam(snthresh = snthresh, noise = noise, ppm = ppm, peakwidth = peakwidth)
  data <- findChromPeaks(data, param = cwp)
  peaks <- chromPeaks(data)
  return(peaks)
}

# file = 'Example/CS52684_neg_SWATH.mzXML'
# features = get_ms1_features(file)
# write.csv(features, 'Example/CS52684_neg_SWATH.features.csv')