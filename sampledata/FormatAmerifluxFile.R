library(tidyverse)
library(stringr)
# format Ameriflux data to create a sample file:
setwd('/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/CarbonExchangeParameters')
tower.data <- read.csv('sampledata/AMF_US-Skr_BASE_HH_2-5.csv', skip=2) %>%
  mutate(Year = str_sub(TIMESTAMP_START, 1, 4),
Month = str_sub(TIMESTAMP_START, 5,6),
Day = str_sub(TIMESTAMP_START, 7,8),
Hour = str_sub(TIMESTAMP_START, 9,10),
Minutes = str_sub(TIMESTAMP_START, 11,12),
Date = paste( Year, Month, Day, sep ="-" ),
Time = paste( Hour, Minutes, sep =":" ),
YearMon = paste(Year, Month, sep="-"))

write.csv(tower.data, 'AMF_US-Skr_BASE_HH_2-5_Formatted.csv' )
