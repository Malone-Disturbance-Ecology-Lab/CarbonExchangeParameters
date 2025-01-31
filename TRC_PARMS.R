#' fit a Temperature Response Curves (TRC) by an index to get a parameter file
#'
#' `TRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `TRC_PARMS` requires a few arguments @data.frame,  @iterations, and @priors.
#' `dataframe` is a dataframe that containes nee, an index, and TA.
#' `iterations` is the number of iterations to run brms. 
#' `priors ` is the priors for the brms to use.

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)


TRC_PARMS <- function( data.frame, iterations, priors, idx){
  
  data.frame <- data.frame %>% mutate(idx= idx)
  
  try(equation <- nee ~ a * exp(b*TA) , silent =T)
  
  # PARM Dataframe:
  try(parms <- data.frame(idx=as.character(), 
                      a.mean = as.numeric(), 
                      a.se = as.numeric(),
                      a.Bulk_ESS= as.numeric() ,
                      a.Tail_ESS= as.numeric() ,
                      a.Rhat= as.numeric(),
                      
                      b.mean = as.numeric(),
                      b.se = as.numeric(),
                      b.Bulk_ESS= as.numeric() ,
                      b.Tail_ESS= as.numeric() ,
                      b.Rhat= as.numeric(),
                      samples= as.numeric()), silent = T)
  
  for ( i in unique(idx)){
    print(i)
    
    # Subset the file:
    try( df <- data.frame %>% filter(idx == i), silent= T)
    # get priors:
    
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    try(model.brms <-  brm( bf(nee ~ a * exp(b*TA),a+b ~ 1, nl=TRUE),
                        prior = priors , data = data.frame, 
                        backend = "cmdstanr", iter = iterations, cores =4, seed=101), silent= T)
    
    
    try(model.brms.df <- summary(model.brms)$fixed , silent= T)
    
    try(model.brms.df.a <- model.brms.df %>% filter( row.names(model.brms.df) == 'a_Intercept'), silent= T)
    try(model.brms.df.b <- model.brms.df %>% filter( row.names(model.brms.df) == 'b_Intercept'), silent= T)
    
    try(samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow , silent= T)
    try(baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric, silent= T)
    
    
    try(results <- data.frame( idx = i, 
                           a.mean = model.brms.df.a$Estimate ,
                           a.se = model.brms.df.a$Est.Error ,
                           a.Bulk_ESS = model.brms.df.a$Bulk_ESS , 
                           a.Tail_ESS = model.brms.df.a$Tail_ESS ,
                           a.Rhat = model.brms.df.a$Rhat ,
                           
                           b.mean = model.brms.df.b$Estimate ,
                           b.se = model.brms.df.b$Est.Error ,
                           b.Bulk_ESS = model.brms.df.b$Bulk_ESS , 
                           b.Tail_ESS = model.brms.df.b$Tail_ESS ,
                           b.Rhat = model.brms.df.b$Rhat,
                           samples= samples/baseline *100 ), silent= T)
    try( parms <- parms %>% rbind(results), silent= T)
  }
  
  return(parms ) 
}

# Example of priors: 

priors <- prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.09)
