# Example Workflow:

rm(list=ls())

library(tidyverse)

source("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/CarbonExchangeParameters/SampleData/05_Flow_SRS6Fix.R" )

srs6 <- srs6.sites %>% mutate(PAR= SW_IN, nee=NEE, YearMon = format( TIMESTAMP, format="%Y-%m")) 

source('LRC_PARMS.R' )

srs6.sub <- srs6 %>% filter(YearMon =="2009-07") # Subset the datafile 

test.LRC.parms <- LRC_PARMS(data.frame=srs6.sub , 
                            idx.colname='YearMon', 
                            priors = priors.lrc, 
                            iterations = 4000,
                            NEE.colname ='NEE', 
                            PAR.colname ='PAR')

write.csv(test.LRC.parms , 'SRS6-Yearmon-LRC-Parms.csv' )

test.LRC.parms <-read.csv('SRS6-Yearmon-LRC-Parms.csv')


priors.trc <- prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.09)

source('TRC_PARMS.R' )

test.TRC.parms <- TRC_PARMS(data.frame=srs6.sub , 
                            idx.colname='YearMon', 
                            priors = priors.trc, 
                            iterations = 4000,
                            NEE.colname='NEE', 
                            PAR.colname='PAR',
                            TA.colname = 'TA')

write.csv(test.TRC.parms , 'SRS6-Yearmon-TRC-Parms.csv' )



