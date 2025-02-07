# Example Workflow:
library(tidyverse)
library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)

rm(list=ls())
source("/Users/bz294/Documents/Proposals/GPP/Code/CarbonExchangeParameters/SampleData/05_Flow_SRS6Fix.R" )

srs6 <- srs6.sites %>% mutate(PAR= SW_IN, nee=NEE, YearMon = format( TIMESTAMP, format="%Y-%m")) 

parms <- data.frame(idx=as.character(),
                    a.mean = as.numeric(),
                    a.se = as.numeric(),
                    a.Bulk_ESS= as.numeric() ,
                    a.Tail_ESS= as.numeric() ,
                    a.Rhat= as.numeric(),
                    
                    B.mean = as.numeric(),
                    B.se = as.numeric(),
                    B.Bulk_ESS= as.numeric() ,
                    B.Tail_ESS= as.numeric() ,
                    B.Rhat= as.numeric(),
                    
                    Y.mean = as.numeric(),
                    Y.se= as.numeric(),
                    Y.Bulk_ESS= as.numeric() ,
                    Y.Tail_ESS= as.numeric() ,
                    Y.Rhat= as.numeric(),
                    
                    Z.mean = as.numeric(),
                    Z.se= as.numeric(),
                    Z.Bulk_ESS= as.numeric() ,
                    Z.Tail_ESS= as.numeric() ,
                    Z.Rhat= as.numeric(),
                    samples= as.numeric())

priors <-  
  prior(normal(-0.01, 0.1), nlpar = "a", lb=-0.2, ub=0) +   
  prior(normal(0, 10), nlpar = "B") +   
  prior(normal(0, 10), nlpar = "Y") +  
  prior(normal(0, 10), nlpar = "Z")

#equation <- nee ~ (a1 * PAR * ax)/(a1 * PAR + ax) - r
equation <- nee ~ a * exp(-B * PAR) - Y * exp(-Z * PAR)

model.brms <- brm( bf( equation, a+B+Y+Z ~ 1, nl=TRUE),
                   prior = priors , data = srs6 %>% filter(YearMon == '2004-12') , 
                   iter = 4000, cores =6,  backend = "cmdstanr")

model.brms.df <- summary(model.brms)$fixed 

model.brms.df.a <- model.brms.df %>% filter( row.names(model.brms.df) == 'a_Intercept')
model.brms.df.B <- model.brms.df %>% filter( row.names(model.brms.df) == 'B_Intercept')
model.brms.df.Z <- model.brms.df %>% filter( row.names(model.brms.df) == 'Z_Intercept')
model.brms.df.Y <- model.brms.df %>% filter( row.names(model.brms.df) == 'Y_Intercept')

samples <- srs6 %>% filter(YearMon == '2004-12')%>% select(nee)  %>% na.omit %>% nrow
baseline <- as.Date(paste('2004-12', '-01', sep="")) %>% lubridate::days_in_month()*48 %>% as.numeric

results <- data.frame( idx = '2004-12', 
                       a.mean = model.brms.df.a$Estimate ,
                       a.se = model.brms.df.a$Est.Error ,
                       a.Bulk_ESS = model.brms.df.a$Bulk_ESS , 
                       a.Tail_ESS = model.brms.df.a$Tail_ESS ,
                       a.Rhat = model.brms.df.a$Rhat ,
                       
                       B.mean = model.brms.df.B$Estimate ,
                       B.se = model.brms.df.B$Est.Error ,
                       B.Bulk_ESS = model.brms.df.B$Bulk_ESS , 
                       B.Tail_ESS = model.brms.df.B$Tail_ESS ,
                       B.Rhat = model.brms.df.B$Rhat ,
                       
                       Y.mean = model.brms.df.Y$Estimate ,
                       Y.se = model.brms.df.Y$Est.Error ,
                       Y.Bulk_ESS = model.brms.df.Y$Bulk_ESS , 
                       Y.Tail_ESS = model.brms.df.Y$Tail_ESS ,
                       Y.Rhat = model.brms.df.Y$Rhat,
                       
                       Z.mean = model.brms.df.Z$Estimate ,
                       Z.se = model.brms.df.Z$Est.Error ,
                       Z.Bulk_ESS = model.brms.df.Z$Bulk_ESS , 
                       Z.Tail_ESS = model.brms.df.Z$Tail_ESS ,
                       Z.Rhat = model.brms.df.Z$Rhat,
                       samples= samples/baseline *100)

parms <- parms %>% rbind(results)
# needs nee and PAR

source('LRC_PARMS_Mod_Exponential.R')

test.LRC.parms <- LRC_PARMS_Mod_Exponential(data.frame=srs6 , 
                            idx=srs6$YearMon, 
                            priors = priors, 
                            iterations = 4000,
                            nee='nee', 
                            PAR='PAR')

write.csv(test.LRC.parms , 'SRS6-Yearmon-LRC-Parms.csv' )

test.LRC.parms <-read.csv('SRS6-Yearmon-LRC-Parms.csv')


priors.trc <- prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.09)

source('TRC_PARMS.R' )

test.TRC.parms <- TRC_PARMS(data.frame=srs6 , 
                            idx=srs6$YearMon, 
                            priors = priors.trc, 
                            iterations = 4000,
                            nee='nee', 
                            PAR='PAR',
                            TA = 'TA')

write.csv(test.TRC.parms , 'SRS6-Yearmon-TRC-Parms.csv' )



