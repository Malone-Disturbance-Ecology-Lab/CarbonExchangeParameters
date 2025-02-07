1# Example Workflow:
library(tidyverse)
library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)

source("/Users/sm3466/YSE Dropbox/Sparkle Malone/Research/CarbonExchangeParameters/SampleData/05_Flow_SRS6Fix.R" )

srs6 <- srs6.sites %>% mutate(PAR= SW_IN, nee=NEE, YearMon = format( TIMESTAMP, format="%Y-%m")) 

parms <- data.frame(idx=as.character(), 
                    a1.mean = as.numeric(), 
                    a1.se = as.numeric(),
                    a1.Bulk_ESS= as.numeric() ,
                    a1.Tail_ESS= as.numeric() ,
                    a1.Rhat= as.numeric(),
                    
                    ax.mean = as.numeric(),
                    ax.se = as.numeric(),
                    ax.Bulk_ESS= as.numeric() ,
                    ax.Tail_ESS= as.numeric() ,
                    ax.Rhat= as.numeric(),
                    
                    r.mean = as.numeric(),
                    r.se= as.numeric(),
                    r.Bulk_ESS= as.numeric() ,
                    r.Tail_ESS= as.numeric() ,
                    r.Rhat= as.numeric(),
                    samples= as.numeric(),
                    
                    theta.mean = as.numeric(),
                    theta.se= as.numeric(),
                    theta.Bulk_ESS= as.numeric() ,
                    theta.Tail_ESS= as.numeric() ,
                    theta.Rhat= as.numeric()
)

priors <-  prior(normal(-0.01, 0.1), nlpar = "a1", lb=-0.2, ub= 0) +
  prior(normal( -7.65 ,  0.33), nlpar = "ax", lb=-30, ub= -5) +
  prior(normal(2.10, 0.11), nlpar = "r", lb=1.9, ub= 2.2)+
  prior(normal(1.9, .1), nlpar = "theta", lb=0, ub= 3)


equation <- nee ~ ((a1 * PAR) + ax - sqrt((a1 * PAR + ax)^2 - 4 * (a1 * PAR * theta * ax))/(2 * theta))-r


model.brms <- brm( bf( equation, a1+ax+r+theta ~ 1, nl=TRUE),
                   prior = priors , data = srs6 %>% filter(YearMon == '2004-12') , 
                   iter = 4000, cores =6,  backend = "cmdstanr")

model.brms.df <- summary(model.brms)$fixed 

model.brms.df.a1 <- model.brms.df %>% filter( row.names(model.brms.df) == 'a1_Intercept')
model.brms.df.ax <- model.brms.df %>% filter( row.names(model.brms.df) == 'ax_Intercept')
model.brms.df.r <- model.brms.df %>% filter( row.names(model.brms.df) == 'r_Intercept')
model.brms.df.theta <- model.brms.df %>% filter( row.names(model.brms.df) == 'theta_Intercept')

samples <- srs6 %>% filter(YearMon == '2004-12')%>% select(nee)  %>% na.omit %>% nrow
baseline <- as.Date(paste('2004-12', '-01', sep="")) %>% lubridate::days_in_month()*48 %>% as.numeric

results <- data.frame( idx = '2004-12', 
                       a1.mean = model.brms.df.a1$Estimate ,
                       a1.se = model.brms.df.a1$Est.Error ,
                       a1.Bulk_ESS = model.brms.df.a1$Bulk_ESS , 
                       a1.Tail_ESS = model.brms.df.a1$Tail_ESS ,
                       a1.Rhat = model.brms.df.a1$Rhat ,
                       
                       ax.mean = model.brms.df.ax$Estimate ,
                       ax.se = model.brms.df.ax$Est.Error ,
                       ax.Bulk_ESS = model.brms.df.ax$Bulk_ESS , 
                       ax.Tail_ESS = model.brms.df.ax$Tail_ESS ,
                       ax.Rhat = model.brms.df.ax$Rhat ,
                       
                       r.mean = model.brms.df.r$Estimate ,
                       r.se = model.brms.df.r$Est.Error ,
                       r.Bulk_ESS = model.brms.df.r$Bulk_ESS , 
                       r.Tail_ESS = model.brms.df.r$Tail_ESS ,
                       r.Rhat = model.brms.df.r$Rhat,
                       
                       theta.mean = model.brms.df.r$Estimate ,
                       theta.se = model.brms.df.r$Est.Error ,
                       theta.Bulk_ESS = model.brms.df.r$Bulk_ESS , 
                       theta.Tail_ESS = model.brms.df.r$Tail_ESS ,
                       theta.Rhat = model.brms.df.r$Rhat,
                       
                       samples= samples/baseline *100)

parms <- parms %>% rbind(results)
# needs nee and PAR









source('LRC_PARMS_NRH.R')

test.LRC.parms <- LRC_PARMS(data.frame=srs6 , 
                            idx=srs6$YearMon, 
                            priors = priors.lrc, 
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



