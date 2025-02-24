# Example Workflow:
library(tidyverse)
library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
#
rm(list=ls())
#
source("SampleData/05_Flow_SRS6Fix.R" )
#
srs6 <- srs6.sites %>% mutate(PAR= SW_IN, nee=NEE, YearMon = format( TIMESTAMP, format="%Y-%m")) %>% filter(year(TIMESTAMP) < 2007)
#
####
source('LRC_PARMS_06.R')

test.LRC.parms.06 <- LRC_PARMS_06(data.frame=srs6 ,
                            priors = priors.lrc, 
                            iterations = 4000,
                            idx.colname="YearMon", 
                            NEE.colname='nee', 
                            PAR.colname='PAR')

write.csv(test.LRC.parms.06 , 'SRS6-Yearmon-LRC-Parms_06.csv' )

test.LRC.parms.06 <-read.csv('SRS6-Yearmon-LRC-Parms_06.csv')

#### test the temperature response curve
source('TRC_PARMS_02.R')
test.TRC.parms.02 <- TRC_PARMS_02(data.frame=srs6,
                                  priors.trc = priors.trc, 
                                  iterations = 4000,
                                  idx.colname="YearMon", 
                                  NEE.colname='nee', 
                                  PAR.colname='PAR',
                                  TA.colname="TA")

write.csv(test.TRC.parms.02 , 'SRS6-Yearmon-TRC-Parms_02.csv')

test.TRC.parms.02 <-read.csv('SRS6-Yearmon-TRC-Parms_02.csv')


