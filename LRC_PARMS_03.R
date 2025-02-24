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
        
             nee ~ ax * (1 - exp(-a1 * PAR / ax)) - r  (exponential equation)
  
          or nee ~ Pmax*(1-exp(-a*PAR/Pmax))-Rd  (equation in paper format)
          where NEE= net ecosystem exchange in 𝜇mol m-2 s-1
                ax=Pmax is light-saturated net photosynthetic rate, a parameter
                a1=a is the initial quantum efficiency or intrinsic quantum yield , a parameter
                PAR = Photosynthetically active radiation in mol m-2 s-1
                r=Rd is dark respiration rate, a parameter

  
          The equations require PAR (𝜇mol m-2 s-1) and NEE (𝜇mol m-2 s-1)
          
          reference https://link.springer.com/article/10.1007/s12010-015-1543-0
          Equation # 8 in the paper (reference)
                ")


# Example of priors: 

priors.lrc <-  prior(normal(-9.9, 0.4), nlpar = "ax", lb=-15, ub= 1) +
  prior(normal( -0.05 ,  0.1), nlpar = "a1", lb=-0.2, ub= 0.5) +
  prior(normal(1.2, 0.11), nlpar = "r", lb=1.0, ub= 3.2)



message("To see the default priors run: 'priors.lrc' ")


LRC_PARM_03 <- function( data.frame, iterations, priors.lrc, idx.colname, NEE.colname , PAR.colname ){
  
  data.frame$nee <- data.frame[,NEE.colname]
  data.frame$idx <- data.frame[,idx.colname]
  data.frame$PAR <- data.frame[,PAR.colname]
  
  df <- data.frame %>% select(idx,
                              nee,
                              PAR)
  
  equation <- nee ~ ax * (1 - exp(-a1 * PAR / ax)) - r

  
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
  # PARM Dataframe:
  
  parms <- data.frame(idx=as.character(), 
                      ax.mean = as.numeric(), 
                      ax.se = as.numeric(),
                      ax.Bulk_ESS= as.numeric() ,
                      ax.Tail_ESS= as.numeric() ,
                      ax.Rhat= as.numeric(),
                      
                      a1.mean = as.numeric(),
                      a1.se = as.numeric(),
                      a1.Bulk_ESS= as.numeric() ,
                      a1.Tail_ESS= as.numeric() ,
                      a1.Rhat= as.numeric(),
                      
                      r.mean = as.numeric(),
                      r.se= as.numeric(),
                      r.Bulk_ESS= as.numeric() ,
                      r.Tail_ESS= as.numeric() ,
                      r.Rhat= as.numeric(),
                      samples= as.numeric()
                      
                      
  )
  
  
  message(" Your dataframe looks good and you are now ready to start fitting models")
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- df %>% filter(idx == i, PAR > 0) %>% na.omit
    
    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())
    
    try(model.brms <- brm( bf( equation, ax+a1+r ~ 1, nl=TRUE),
                           prior = priors.lrc , data = df, iter = iterations, cores =3, chains = 1, backend = "cmdstanr"), silent= F)
    
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent = T)
    try(model.brms.df.ax <- model.brms.df %>% filter( row.names(model.brms.df) == 'ax_Intercept'), silent = F)
    try( model.brms.df.a1 <- model.brms.df %>% filter( row.names(model.brms.df) == 'a1_Intercept'), silent = F)
    try(model.brms.df.r <- model.brms.df %>% filter( row.names(model.brms.df) == 'r_Intercept'), silent = F)

    samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow
    
    baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric
    
    
    try(results <- data.frame( idx = i, 
                               ax.mean = model.brms.df.ax$Estimate ,
                               ax.se = model.brms.df.ax$Est.Error ,
                               ax.Bulk_ESS = model.brms.df.ax$Bulk_ESS , 
                               ax.Tail_ESS = model.brms.df.ax$Tail_ESS ,
                               ax.Rhat = model.brms.df.ax$Rhat ,
                               
                               a1.mean = model.brms.df.a1$Estimate ,
                               a1.se = model.brms.df.a1$Est.Error ,
                               a1.Bulk_ESS = model.brms.df.a1$Bulk_ESS , 
                               a1.Tail_ESS = model.brms.df.a1$Tail_ESS ,
                               a1.Rhat = model.brms.df.a1$Rhat ,
                               
                               r.mean = model.brms.df.r$Estimate ,
                               r.se = model.brms.df.r$Est.Error ,
                               r.Bulk_ESS = model.brms.df.r$Bulk_ESS , 
                               r.Tail_ESS = model.brms.df.r$Tail_ESS ,
                               r.Rhat = model.brms.df.r$Rhat,
                               
                               samples= samples/baseline *100), silent = T)
    
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
