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

message("Using the LRC_PARMS function requires the following libraries: brms, 
        cmdstanr, ggplot2, beepr, tidybayes, tidyverse, and ggpubr.")

message("This function uses the equation:
        
        nee ~  r0 * exp(alpha*TA + beta*TA*TA)
        
        Where nee is nighttime nee, that is, ecosystem respiration (𝜇mol CO2 m-2 s-1), 
        TA is air temperature (°C), 
        r0 is ecosystem respiration at TA = 0 °C,
        alpha and beta are coefficients, 
        The equations require TA (°C),  PAR (𝜇mol m-2 s-1) and NEE (𝜇mol m-2 s-1), 
        ")


# Example of priors: 
priors.trc <-  prior(normal(1, 5), nlpar = "r0", lb=0, ub= 5) +
  prior(normal(0.07, 0.1), nlpar = "alpha", lb=0, ub= 0.2) +
  prior(normal(0, 0.05), nlpar = "beta", lb=-0.1, ub= 0.1)

message("To see the default priors run: 'priors.trc' ")

TRC_PARMS_02 <- function( data.frame, iterations, priors.trc, idx.colname, NEE.colname, PAR.colname, TA.colname){
  
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
  
  try(equation <- nee ~ r0 * exp(alpha*TA + beta*TA*TA) , silent =T)
  
  # PARM Dataframe:
  parms <- data.frame(idx=as.character(), 
                      r0.mean = as.numeric(), 
                      r0.se = as.numeric(),
                      r0.Bulk_ESS= as.numeric() ,
                      r0.Tail_ESS= as.numeric() ,
                      r0.Rhat= as.numeric(),
                      
                      alpha.mean = as.numeric(),
                      alpha.se = as.numeric(),
                      alpha.Bulk_ESS= as.numeric() ,
                      alpha.Tail_ESS= as.numeric() ,
                      alpha.Rhat= as.numeric(),
                      
                      beta.mean = as.numeric(),
                      beta.se= as.numeric(),
                      beta.Bulk_ESS= as.numeric() ,
                      beta.Tail_ESS= as.numeric() ,
                      beta.Rhat= as.numeric(),
                      samples= as.numeric())
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    try( df.sub <- df %>% filter(idx == i, PAR < 10), silent= T)
    # get priors:
    
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    try(model.brms <- brm( bf( equation, r0+alpha+beta ~ 1, nl=TRUE), 
                           prior = priors.trc, data = df.sub, 
                           backend = "cmdstanr", iter = iterations, cores =4, seed=101), silent= F)
    # Here it should be df.sub
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent = T)
    try(model.brms.df.r0 <- model.brms.df %>% filter( row.names(model.brms.df) == 'r0_Intercept'), silent = F)
    try(model.brms.df.alpha <- model.brms.df %>% filter( row.names(model.brms.df) == 'alpha_Intercept'), silent = F)
    try(model.brms.df.beta <- model.brms.df %>% filter( row.names(model.brms.df) == 'beta_Intercept'), silent = F)
    
    samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow
    
    baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric
    
    try(results <- data.frame( idx = i, 
                               r0.mean = model.brms.df.r0$Estimate ,
                               r0.se = model.brms.df.r0$Est.Error ,
                               r0.Bulk_ESS = model.brms.df.r0$Bulk_ESS , 
                               r0.Tail_ESS = model.brms.df.r0$Tail_ESS ,
                               r0.Rhat = model.brms.df.r0$Rhat ,
                               
                               alpha.mean = model.brms.df.alpha$Estimate ,
                               alpha.se = model.brms.df.alpha$Est.Error ,
                               alpha.Bulk_ESS = model.brms.df.alpha$Bulk_ESS , 
                               alpha.Tail_ESS = model.brms.df.alpha$Tail_ESS ,
                               alpha.Rhat = model.brms.df.alpha$Rhat ,
                               
                               beta.mean = model.brms.df.beta$Estimate ,
                               beta.se = model.brms.df.beta$Est.Error ,
                               beta.Bulk_ESS = model.brms.df.beta$Bulk_ESS , 
                               beta.Tail_ESS = model.brms.df.beta$Tail_ESS ,
                               beta.Rhat = model.brms.df.beta$Rhat,
                               samples= samples/baseline *100 ), silent = T)
    
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

# EOF
