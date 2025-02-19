#' fit a Temperature Response Curves (TRC) by an index to get a parameter file
#'
#' `TRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `TRC_PARMS` requires a few arguments @data.frame,  @iterations, @priors, @idx.colname, @NEE.colname, @PAR.colname, and @TA.colname.
#' `dataframe` is a dataframe that containes NEE, idx, TA, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors.trc ` is the priors for the brms to use.
#' `idx.colname` is the column containing the index.
#' `nee.colname` is the name of the column containing NEE.
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
        
        Reco ~ Rref * exp(E0(1/(Tref - T0) - 1/(Tair - T0)))")



# Example of priors: 
priors.trc <- prior(normal( 0.5 ,  0.3), nlpar = "E0", lb=0.01, ub= 1)+
  prior(normal( 1.0 ,  0.3), nlpar = "Rref", lb=0.01, ub= 1)

message("To see the default priors run: 'priors.trc' ")

TRC_PARMS_01 <- function(data.frame, iterations, priors, idx.colname, nee.colname, TA.colname, PAR.colname){
  
  data.frame$nee <- data.frame[,nee.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$TA <- data.frame[,TA.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  
  df <- data.frame %>% select(idx, nee, TA, PAR)
  
  if( c("idx") %in% names(df) ) {
    print("GREAT JOB! your dataframe contains idx")
  } else{ print("The dataframe must include: idx,  Reco, TA")}
  
  if( c("Reco") %in% names(df) ) {
    print("YIPEE! your dataframe contains Reco")
  }  else{ print("The dataframe must include: idx,  Reco, TA")}
  
  if( c("PAR") %in% names(df) ) {
    print("Hooray! your dataframe contains PAR")
  } else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  if( c("TA") %in% names(df) ) {
    print("Hooray! your dataframe contains TA")
  } else{ print("The dataframe must include: idx,  Reco, TA")}
  
  try(equation <- NEE ~ Rref * exp(E0 * (1 / (10 + 46.02) - 1 / (TA + 46.02))) , silent =T)
  
  # PARM Dataframe:
  try(parms <- data.frame(idx=as.character(), 
                          E0.mean = as.numeric(), 
                          E0.se = as.numeric(),
                          E0.Bulk_ESS= as.numeric() ,
                          E0.Tail_ESS= as.numeric() ,
                          E0.Rhat= as.numeric(),
                          
                          Rref.mean = as.numeric(),
                          Rref.se = as.numeric(),
                          Rref.Bulk_ESS= as.numeric() ,
                          Rref.Tail_ESS= as.numeric() ,
                          Rref.Rhat= as.numeric(),
                          samples= as.numeric()), silent = T)
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- try(df %>% filter(idx == i, PAR < 10), silent = TRUE) ## to get night time data
    
    # get priors:
    print(df.sub)
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    try(model.brms <-  brm( bf(nee ~ Rref * exp(E0 * (1 / (10 + 46.02) - 1 / (TA + 46.02))),E0+Rref ~ 1, nl=TRUE),
                            prior = priors.trc , data = df.sub, 
                            backend = "cmdstanr", iter = iterations, cores =4, seed=101), silent= F)
    
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent= F)
    
    # Extract parameter estimates correctly
    model.brms.df.E0 <- model.brms.df %>% filter(rownames(model.brms.df) == "E0_Intercept")
    model.brms.df.Rref <- model.brms.df %>% filter(rownames(model.brms.df) == "Rref_Intercept")
    
    try(samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow , silent= F)
    try(baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric, silent= F)
    
    
    try(results <- data.frame( idx = i, 
                               E0.mean = model.brms.df.E0$Estimate ,
                               E0.se = model.brms.df.E0$Est.Error ,
                               E0.Bulk_ESS = model.brms.df.E0$Bulk_ESS , 
                               E0.Tail_ESS = model.brms.df.E0$Tail_ESS ,
                               E0.Rhat = model.brms.df.E0$Rhat ,
                               
                               Rref.mean = model.brms.df.Rref$Estimate ,
                               Rref.se = model.brms.df.Rref$Est.Error ,
                               Rref.Bulk_ESS = model.brms.df.Rref$Bulk_ESS , 
                               Rref.Tail_ESS = model.brms.df.Rref$Tail_ESS ,
                               Rref.Rhat = model.brms.df.Rref$Rhat,
                               samples= samples/baseline *100 ), silent= T)
    
    message( 'YOUR DID IT!')
    print(summary(model.brms))
    #print(results) 
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