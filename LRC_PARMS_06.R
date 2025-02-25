#' fit a Light Response Curves (LRC) by an index to get a parameter file
#'
#' `LRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `LRC_PARMS` requires a few arguments @data.frame,  @iterations, @priors, @idx.colname, @NEE.colname, @PAR.colname, and @TA.colname.
#' `dataframe` is a dataframe that contains nee, an index, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors.LRC ` is the priors for the brms to use.
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
        
        NEE ~ Pmax * (1 - exp(alpha * (PAR - Icomp)))
        
        Where nee is daytime net ecosystem exchange (𝜇mol CO2 m-2 s-1), 
        Pmax: the maximum photosynthesis rate (𝜇mol CO2 m-2 s-1), 
        alpha is a coefficient,
        Icomp is the compensation light level with the same unit of PAR.
        
        The equations require PAR (𝜇mol m-2 s-1) and NEE (𝜇mol m-2 s-1)")


priors.lrc <-  prior(normal(-5, 5), nlpar = "Pmax", lb=-30, ub= 0) +
  prior(normal(0.1, 1), nlpar = "alpha", lb=0, ub= 1) +
  prior(normal(300, 100), nlpar = "Icomp", lb=0, ub= 500)

message("To see the default priors run: 'priors.lrc' ")


LRC_PARMS_06 <- function( data.frame, iterations, priors.lrc, idx.colname, NEE.colname , PAR.colname ){
  
  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  
  df <- data.frame %>% select(idx,
                              nee,
                              PAR)
  
  equation <- nee ~ Pmax * (1 - exp(alpha * (PAR - Icomp)))
  
  
  if( c("idx") %in% names(df) ) {
    print("GREAT JOB! your dataframe contains idx")
  } else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  if( c("nee") %in% names(df) ) {
    print("YIPEE! your dataframe contains nee")
  }  else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  if( c("PAR") %in% names(df) ) {
    print("Hooray! your dataframe contains PAR")
  } else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  # PARM Dataframe:
  parms <- data.frame(idx=as.character(), 
                      Pmax.mean = as.numeric(), 
                      Pmax.se = as.numeric(),
                      Pmax.Bulk_ESS= as.numeric() ,
                      Pmax.Tail_ESS= as.numeric() ,
                      Pmax.Rhat= as.numeric(),
                      
                      alpha.mean = as.numeric(),
                      alpha.se = as.numeric(),
                      alpha.Bulk_ESS= as.numeric() ,
                      alpha.Tail_ESS= as.numeric() ,
                      alpha.Rhat= as.numeric(),
                      
                      Icomp.mean = as.numeric(),
                      Icomp.se= as.numeric(),
                      Icomp.Bulk_ESS= as.numeric() ,
                      Icomp.Tail_ESS= as.numeric() ,
                      Icomp.Rhat= as.numeric(),
                      samples= as.numeric())
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- df %>% filter(idx == i, PAR > 0) %>% na.omit
    
    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())
    
    try(model.brms <- brm( bf( equation, Pmax+alpha+Icomp ~ 1, nl=TRUE), prior = priors.lrc, 
                           data = df.sub, iter = iterations, cores =3, chains = 1, backend = "cmdstanr"), silent= F)
    # Here it should be df.sub
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent = T)
    try(model.brms.df.Pmax <- model.brms.df %>% filter( row.names(model.brms.df) == 'Pmax_Intercept'), silent = F)
    try(model.brms.df.alpha <- model.brms.df %>% filter( row.names(model.brms.df) == 'alpha_Intercept'), silent = F)
    try(model.brms.df.Icomp <- model.brms.df %>% filter( row.names(model.brms.df) == 'Icomp_Intercept'), silent = F)
    
    samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow
    
    baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric
    
    try(results <- data.frame( idx = i, 
                               Pmax.mean = model.brms.df.Pmax$Estimate ,
                               Pmax.se = model.brms.df.Pmax$Est.Error ,
                               Pmax.Bulk_ESS = model.brms.df.Pmax$Bulk_ESS , 
                               Pmax.Tail_ESS = model.brms.df.Pmax$Tail_ESS ,
                               Pmax.Rhat = model.brms.df.Pmax$Rhat ,
                               
                               alpha.mean = model.brms.df.alpha$Estimate ,
                               alpha.se = model.brms.df.alpha$Est.Error ,
                               alpha.Bulk_ESS = model.brms.df.alpha$Bulk_ESS , 
                               alpha.Tail_ESS = model.brms.df.alpha$Tail_ESS ,
                               alpha.Rhat = model.brms.df.alpha$Rhat ,
                               
                               Icomp.mean = model.brms.df.Icomp$Estimate ,
                               Icomp.se = model.brms.df.Icomp$Est.Error ,
                               Icomp.Bulk_ESS = model.brms.df.Icomp$Bulk_ESS , 
                               Icomp.Tail_ESS = model.brms.df.Icomp$Tail_ESS ,
                               Icomp.Rhat = model.brms.df.Icomp$Rhat,
                               samples= samples/baseline *100 ), silent = T)
    
    message( 'YOUR DID IT!')
    
    print(results)
    
    try(parms <- parms %>% rbind(results), silent = T)
    
    message( 'Just keep swimming')
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
