#' fit a Temperature Response Curves (TRC) by an index to get a parameter file
#'
#' `TRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `TRC_PARMS` requires a few arguments @data.frame,  @iterations, @priors, @idx.colname, @NEE.colname, @PAR.colname, and @TA.colname.
#' `dataframe` is a dataframe that containes NEE, idx, TA, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors.trc ` is the priors for the brms to use.
#' `idx.colname` is the column containing the index.
#' `NEE.colname` is the name of the column containing NEE.
#' `PAR.colname` is the name of the column containing PAR.
#' `TA.colname` is the name of the column containing air temperature.

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)

message("Using the TRC_PARMS function requires the following libraries: brms, 
        cmdstanr, ggplot2, beepr, tidybayes, tidyverse, and ggpubr.")

message("This function uses the equation:
        
        nee ~ a * exp(b*TA)
  
        
        The equation requires NEE (𝜇mol m-2 s-1) andair temperature")



# Example of priors: 
priors.ts1.trc <-  prior(normal(0.2 , 1), nlpar = "a", lb=0.1, ub= 1) +
  prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.9)

message("To see the default priors run: 'priors.trc' ")

TRC_PARMS <- function( data.frame, iterations, priors.trc, idx.colname, NEE.colname, PAR.colname, TA.colname){
  
  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  data.frame$TA <- data.frame[,TA.colname]
  
  df <- data.frame %>% select(idx, nee, PAR, TA)
  
  if( c("idx") %in% names(df) ) {
    print("GREAT JOB! your dataframe contains idx")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("nee") %in% names(df) ) {
    print("YIPEE! your dataframe contains nee")
  }  else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("PAR") %in% names(df) ) {
    print("Hooray! your dataframe contains PAR")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("TA") %in% names(df) ) {
    print("Hooray! your dataframe contains TA")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
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
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    try( df.sub <- df %>% filter(idx == i, PAR < 10), silent= T)
    # get priors:
    
    try(model.brms <-  brm( bf(nee ~ a * exp(b*TA),a+b ~ 1, nl=TRUE),
                        prior = priors.trc , data = df.sub, 
                        backend = "cmdstanr", iter = iterations, cores =4, seed=101), silent= F)
    
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent= F)
    
    try(model.brms.df.a <- model.brms.df %>% filter( row.names(model.brms.df) == 'a_Intercept'), silent= F)
    try(model.brms.df.b <- model.brms.df %>% filter( row.names(model.brms.df) == 'b_Intercept'), silent= F)
    
    try(samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow , silent= F)
    try(baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric, silent= F)
    
    
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
    
    message( 'YOUR DID IT!')
    print( results) 
    try( parms <- parms %>% rbind(results), silent= T)
    
  }
   
  return(parms ) 
}

message("Model parameters Are fit using the R package brms.
        
        Rhat (Potential Scale Reduction Factor):
        
        Indicates how well the different Markov chains in your analysis 
        have converged to the same posterior distribution. 
        
        Ideally, Rhat should be close to 1 for all parameters.
        
        A high Rhat value suggests potential convergence issues and 
        the need to run the chains longer. 
        
        Bulk ESS (Effective Sample Size - Bulk):
        
        Estimates the effective number of independent samples from the 
        central part of the posterior distribution. 
        
        Tail ESS (Effective Sample Size - Tail):
        Estimates the effective number of independent samples from the tails of the posterior distribution. 
        
        Important for assessing the reliability of quantile estimates (e.g., 95% confidence intervals). 
        
        Key points to remember:
        Aim for Rhat close to 1 and high values for both Bulk ESS and Tail ESS." )

#EOF