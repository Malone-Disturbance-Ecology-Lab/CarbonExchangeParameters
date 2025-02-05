#' fit a Temperature Response Curves (TRC) by an index to get a parameter file
#'
#' `TRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `TRC_PARMS` requires a few arguments @data.frame,  @iterations, and @priors.
#' `dataframe` is a dataframe that containes nee, idx, TA, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors.trc ` is the priors for the brms to use.
#' `idx` is the column containing the index.
#' `nee` is the name of the column containing NEE.
#' `PAR` is the name of the column containing PAR.
#' `TA` is the name of the column containing air temperature.

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)

# Example of priors: 
priors.trc <- prior(normal( 0.5 ,  0.03), nlpar = "b", lb=0.001, ub= 0.09)

TRC_PARMS <- function( data.frame, iterations, priors.trc, idx, nee, PAR, TA){
  
  df <- data.frame %>% mutate(idx= idx,
                              nee= nee,
                              PAR= PAR,
                              TA = TA) %>% select(idx, nee, PAR, TA)
  
  if( c("idx") %in% names(df) ) {
    print("idx present")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("nee") %in% names(df) ) {
    print("nee data present")
  }  else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  
  if( c("PAR") %in% names(df) ) {
    print("PAR data present")
  } else{ print("The dataframe must include: idx,  nee, TA, and PAR")}
  if( c("TA") %in% names(df) ) {
    print("PAR data present")
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
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    try( df.sub <- df %>% filter(idx == i, PAR < 10), silent= T)
    # get priors:
    
    # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
    #                    data = df %>% filter(PAR > 0),
    #                    family = poisson())
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
    
    print( results) 
    try( parms <- parms %>% rbind(results), silent= T)
    
  }
   
  return(parms ) 
}

#EOF