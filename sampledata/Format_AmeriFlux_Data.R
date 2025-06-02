# Format Ameriflux data : https://ameriflux.lbl.gov/
library(tidyverse)
library( brms)

setwd('sampledata')
af.data <- read.csv('AMF_US-Skr_BASE_HH_2-5.csv', skip=2 )

af.data[af.data == -9999]<- NA # replace -9999 with NA

# Format the date and time
af.data.formatted <- af.data %>%  mutate( Year = TIMESTAMP_START %>% substr( 1,4), # Year
                                          Month = TIMESTAMP_START %>% substr( 5,6), # month
                                          Day = TIMESTAMP_START %>% substr( 7,8), # day
                                          Hour = TIMESTAMP_START %>% substr( 9,10), # Hour
                                          Mintutes = TIMESTAMP_START %>% substr( 11,12), # minutes
                                          Date = paste(Year,Month,Day, sep = "-" ) %>% as.Date,
                                          Time = paste(Hour,  Mintutes, sep=":"),
                                          YearMon = paste(Year, Month, sep="-")) %>% filter( Month  == '07')


write.csv(af.data.formatted, 'AMF_US-Skr_BASE_HH_2-5_Formatted.csv')
