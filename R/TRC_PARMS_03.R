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
        
        NEE ~ a * ((TA - Topt)^2) + ERmax   (quadratic temperature response function)
        
        Where 
        a    = a parameter estimate of the fitted quadratic function 
        T    = Temperature in celcius 
        Topt = optimum temperature, a parameter, vertex of the parabola
        ERmax =maximum respiration, a parameter, vertex of the parabola
        
        The equations require T (celcius) PAR (𝜇mol m-2 s-1), and NEE (𝜇mol m-2 s-1)
       
    
        Reference: https://www.nature.com/articles/s41559-023-02121-w#Equ1
        
        ")

# Example of priors: 
priors.trc <- prior(normal( 0.08 ,  0.03), nlpar = "a", lb=0,ub= 0.9 )+
  prior(normal( 24 ,  4), nlpar = "Topt", lb=10, ub= 40) +
  prior(normal(3.8, 5), nlpar = "ERmax", lb=0,ub= 10)


 

#priors.trc <- prior(prior(default(), nlpar = "a") + 
 #                     prior(default(), nlpar = "Topt")
  #                  +prior(default(), nlpar = "ERmax"))

message("To see the default priors run: 'priors.trc' ")


#data.frame=srs6.sub 
#idx.colname='YearMon' 
#priors = priors.trc 
#iterations = 4000
#NEE.colname='NEE' 
#TA.colname = 'TA'
#PAR.colname = 'PAR'



TRC_PARMS_03 <- function( data.frame, iterations, priors.trc, idx.colname, NEE.colname, PAR.colname, TA.colname){
  
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
    print("Hooray! your dataframe contains TA, this is a test")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  try(equation <- nee ~ a * ((TA - Topt)^2) + ERmax , silent =T)
  
  
  # PARM Dataframe:
  try(parms <- data.frame(idx=as.character(), 
                          
                          a.mean = as.numeric(), 
                          a.se = as.numeric(),
                          a.Bulk_ESS= as.numeric() ,
                          a.Tail_ESS= as.numeric() ,
                          a.Rhat= as.numeric(),
                          
                          Topt.mean = as.numeric(), 
                          Topt.se = as.numeric(),
                          Topt.Bulk_ESS= as.numeric() ,
                          Topt.Tail_ESS= as.numeric() ,
                          Topt.Rhat= as.numeric(),
                          
                          ERmax.mean = as.numeric(),
                          ERmax.se = as.numeric(),
                          ERmax.Bulk_ESS= as.numeric() ,
                          ERmax.Tail_ESS= as.numeric() ,
                          ERmax.Rhat= as.numeric(),
                          samples= as.numeric()), silent = T)
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- try(df %>% filter(idx == i, PAR < 10), silent = TRUE)
    
    # get priors:
    print(df.sub)

    try(model.brms <-  brm( bf(nee ~ a * ((TA - Topt)^2) + ERmax, a+Topt+ERmax ~ 1, nl=TRUE),
                            prior = priors.trc , data = df.sub, 
                            backend = "cmdstanr", iter = iterations, cores =4, seed=101), silent= F)
    
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent= F)
    
    try(model.brms.df.a <- model.brms.df %>% filter( row.names(model.brms.df) == 'a_Intercept'), silent= F)
    try(model.brms.df.Topt <- model.brms.df %>% filter( row.names(model.brms.df) == 'Topt_Intercept'), silent= F)
    try(model.brms.df.ERmax <- model.brms.df %>% filter( row.names(model.brms.df) == 'ERmax_Intercept'), silent= F)
    
    try(samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow , silent= F)
    try(baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric, silent= F)
    
    
    try(results <- data.frame( idx = i, 
 
                               
                               a.mean = model.brms.df.a$Estimate ,
                               a.se = model.brms.df.a$Est.Error ,
                               a.Bulk_ESS = model.brms.df.a$Bulk_ESS , 
                               a.Tail_ESS = model.brms.df.a$Tail_ESS ,
                               a.Rhat = model.brms.df.a$Rhat ,                              
                               
                               Topt.mean = model.brms.df.Topt$Estimate ,
                               Topt.se = model.brms.df.Topt$Est.Error ,
                               Topt.Bulk_ESS = model.brms.df.Topt$Bulk_ESS , 
                               Topt.Tail_ESS = model.brms.df.Topt$Tail_ESS ,
                               Topt.Rhat = model.brms.df.Topt$Rhat ,
                               
                               ERmax.mean = model.brms.df.ERmax$Estimate ,
                               ERmax.se = model.brms.df.ERmax$Est.Error ,
                               ERmax.Bulk_ESS = model.brms.df.ERmax$Bulk_ESS , 
                               ERmax.Tail_ESS = model.brms.df.ERmax$Tail_ESS ,
                               ERmax.Rhat = model.brms.df.ERmax$Rhat,
                               samples= samples/baseline *100 ), silent= T)
    
    message( 'YOUR DID IT!')
    print(results) 
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
