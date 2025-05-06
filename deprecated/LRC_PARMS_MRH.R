#' fit a Light Response Curves (LRC) by an index to get a parameter file
#'
#' `LRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `LRC_PARMS` requires a few arguments @data.frame,  @iterations, and @priors.
#' `dataframe` is a dataframe that containes nee, an index, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors.LRC ` is the priors for the brms to use.
#' `idx` is the column containing the index.
#' `nee` is the name of the column containing NEE.
#' `PAR` is the name of the column containing PAR.
#' `TA` is the name of the column containing air temperature.


# Add notes to indicate successful or unsuccessful fit. 

rm(list=ls())

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)


# Example of priors: 
priors.lrc <-  prior(normal(-0.01, 0.1), nlpar = "a1", lb=-0.2, ub= 0) +
  prior(normal( -7.65 ,  0.33), nlpar = "ax", lb=-30, ub= -5) +
  prior(normal(2.10, 0.11), nlpar = "r", lb=1.9, ub= 2.2)+
  prior(normal(.3, .2), nlpar = "b", lb=0, ub= 1)

LRC_PARMS_MRH <- function( data.frame, iterations, priors.lrc, idx, nee, PAR){
  
  df <- data.frame %>% mutate(idx= idx,
                              nee= nee,
                              PAR= PAR ) %>% select(idx, nee, PAR)
  
  equation <- nee ~ a1 * ( (1-b*PAR) / (1+(a1/ax)*PAR) ) - r
  
  
  if( c("idx") %in% names(df) ) {
    print("idx present")
  } else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  if( c("nee") %in% names(df) ) {
    print("nee data present")
  }  else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  if( c("PAR") %in% names(df) ) {
    print("PAR data present")
  } else{ print("The dataframe must include: idx,  nee, and PAR")}
  
  # PARM Dataframe:
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
                      
                      b.mean = as.numeric(),
                      b.se= as.numeric(),
                      b.Bulk_ESS= as.numeric() ,
                      b.Tail_ESS= as.numeric() ,
                      b.Rhat= as.numeric()
  )
  
  
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- df %>% filter(idx == i, PAR > 0) %>% na.omit
    
    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())
    
    try(model.brms <- brm( bf( equation, a1+ax+r+b ~ 1, nl=TRUE),
                           prior = priors.lrc , data = df, iter = iterations, cores =3, chains = 1, backend = "cmdstanr"), silent= F)
    
    print(model.brms)
    
    try(model.brms.df <- summary(model.brms)$fixed , silent = T)
    try(model.brms.df.a1 <- model.brms.df %>% filter( row.names(model.brms.df) == 'a1_Intercept'), silent = F)
    try( model.brms.df.ax <- model.brms.df %>% filter( row.names(model.brms.df) == 'ax_Intercept'), silent = F)
    try(model.brms.df.r <- model.brms.df %>% filter( row.names(model.brms.df) == 'r_Intercept'), silent = F)
    try(model.brms.df.b <- model.brms.df %>% filter( row.names(model.brms.df) == 'b_Intercept'), silent = F)
    
    samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow
    
    baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric
    
    
    try(results <- data.frame( idx = i, 
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
                               
                               b.mean = model.brms.df.r$Estimate ,
                               b.se = model.brms.df.r$Est.Error ,
                               b.Bulk_ESS = model.brms.df.r$Bulk_ESS , 
                               b.Tail_ESS = model.brms.df.r$Tail_ESS ,
                               b.Rhat = model.brms.df.r$Rhat,
                               
                               samples= samples/baseline *100), silent = T)
    
    print(results)
    
    try(parms <- parms %>% rbind(results), silent = T)
  }
  
  return(parms ) 
}

# EOF

