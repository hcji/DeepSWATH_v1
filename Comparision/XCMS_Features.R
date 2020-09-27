# xcms version 3.8.1
library(xcms)

get_ms1_features <- function(file, snthresh = 10, prefilter = c(1, 100), ppm = 100, peakwidth = c(5, 60)){
  data <- readMSData(file, mode = "onDisk") 
  cwp <- CentWaveParam(snthresh = snthresh, prefilter = prefilter, ppm = ppm, peakwidth = peakwidth)
  data <- findChromPeaks(data, param = cwp)
  peaks <- chromPeaks(data)
  return(peaks)
}

# MetDIA data
metdia_file <- 'Comparision/MetDIA_Data/data/30STD_mix 330ppb-1.mzML'
metdia_fet <- get_ms1_features(metdia_file, prefilter = c(1, 500))
write.csv(metdia_fet, 'Comparision/MetDIA_Data/results/xcms_ms1_feature.csv')

# MetaboDIA data
metabodia_pos_file <- 'Comparision/MetaboDIA_Data/data/PH697097_pos_SWATH.mzML'
metabodia_neg_file <- 'Comparision/MetaboDIA_Data/data/PH697097_neg_SWATH.mzML'
metabodia_pos_fet <- get_ms1_features(metabodia_pos_file, prefilter = c(1, 500))
metabodia_neg_fet <- get_ms1_features(metabodia_neg_file, prefilter = c(1, 500))
write.csv(metabodia_pos_fet, 'Comparision/MetaboDIA_Data/results/xcms_ms1_feature_pos.csv')
write.csv(metabodia_neg_fet, 'Comparision/MetaboDIA_Data/results/xcms_ms1_feature_neg.csv')

# MetaboKit data
metabodia_pos_file <- 'Comparision/MetaboKit_Data/data/SWATH-POS-Metab-Mix-10ng-ml_r2-SWATH-POS-PathwayMetab-Mix-10ng-ml_r2.mzML'
metabodia_neg_file <- 'Comparision/MetaboKit_Data/data/SWATH-NEG-Metab-Mix-10ng-ml_r2-SWATH-NEG-PathwayMetab-Mix-10ng-ml_r2.mzML'
metabodia_pos_fet <- get_ms1_features(metabodia_pos_file, snthresh = 3, prefilter = c(1, 100))
metabodia_neg_fet <- get_ms1_features(metabodia_neg_file, snthresh = 3, prefilter = c(1, 100))
write.csv(metabodia_pos_fet, 'Comparision/MetaboKit_Data/results/xcms_ms1_feature_pos.csv')
write.csv(metabodia_neg_fet, 'Comparision/MetaboKit_Data/results/xcms_ms1_feature_neg.csv')
