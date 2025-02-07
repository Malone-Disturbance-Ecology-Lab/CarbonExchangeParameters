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
srs6 <- srs6.sites %>% mutate(PAR= SW_IN, nee=NEE, YearMon = format( TIMESTAMP, format="%Y-%m")) 
#
####
source('LRC_PARMS_Junna.R')

test.LRC.parms.Junna <- LRC_PARMS_Junna(data.frame=srs6 , 
                            idx=srs6$YearMon, 
                            priors = priors.lrc, 
                            iterations = 4000,
                            nee='nee', 
                            PAR='PAR')

write.csv(test.LRC.parms.Junna , 'SRS6-Yearmon-LRC-Parms_Junna.csv' )

test.LRC.parms <-read.csv('SRS6-Yearmon-LRC-Parms_Junna.csv')



