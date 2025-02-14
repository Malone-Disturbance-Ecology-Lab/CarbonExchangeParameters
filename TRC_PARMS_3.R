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
        
        Reco ~ a * ((T - Topt)^2) + ERmax
        Where a, Topt (optimum temperature) and ERmax (maximum respiration) are parameters.
        The equations require T (celcius) and Respiration (Reco) (𝜇mol m-2 s-1)
        quadratic temperature response function")

# paper 
#https://www.nature.com/articles/s41559-023-02121-w#Equ1

# Example of priors: 
#priors.trc <- prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.09) ## sample
priors.trc <- prior(normal( 0.1 ,  0.03), nlpar = "a", lb=0.5, ub= 0.09)+
  prior(normal( 20 ,  4), nlpar = "Topt", lb=5, ub= 40) +
  prior(normal(7, 5), nlpar = "ERmax", lb=1, ub= 20)

message("To see the default priors run: 'priors.trc' ")

TRC_PARMS_3 <- function( data.frame, iterations, priors.trc, idx.colname, Reco.colname, TA.colname){
  
  data.frame$Reco <- data.frame[,Reco.colname]
  #data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  #data.frame$PAR <- data.frame[,PAR.colname]
  data.frame$TA <- data.frame[,TA.colname]
  
  df <- data.frame %>% select(idx, Reco,TA)
  
  if( c("idx") %in% names(df) ) {
    print("GREAT JOB! your dataframe contains idx")
  } else{ print("The dataframe must include: idx,  Reco, TA")}
  
  if( c("Reco") %in% names(df) ) {
    print("YIPEE! your dataframe contains nee")
  }  else{ print("The dataframe must include: idx,  Reco, TA")}
  
 # if( c("PAR") %in% names(df) ) {
  #  print("Hooray! your dataframe contains PAR")
  #} else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("TA") %in% names(df) ) {
    print("Hooray! your dataframe contains TA")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  #try(equation <- nee ~ a * exp(b*TA) , silent =T)
  try(equation <- Reco ~ a * ((TA - Topt)^2) + ERmax , silent =T)
  
  
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
    df.sub <- try(df %>% filter(idx == i, TA < 40), silent = TRUE)
    
    # get priors:
    print(df.sub)
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
    try(model.brms <-  brm( bf(Reco ~ a * ((TA - Topt)^2) + ERmax, a+Topt+ERmax ~ 1, nl=TRUE),
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