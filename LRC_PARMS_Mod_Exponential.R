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

#rm(list=ls())

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)


# Example of priors: 

priors.lrc <- prior(normal(-0.01, 0.1), nlpar = "a", lb=-0.2, ub= 0) + 
  prior(default(), nlpar = "B") + 
  prior(default(), nlpar = "Z") + 
  prior(default(), nlpar = "Y")

LRC_PARMS_Mod_Exponential <- function( data.frame, iterations, priors.lrc, idx, nee, PAR){
  
  df <- data.frame %>% mutate(idx= idx,
                                      nee= nee,
                                      PAR= PAR ) %>% select(idx, nee, PAR)
  
  equation <- nee ~ a * exp(-B * PAR) - Y * exp(-Z * PAR)
  
  
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
                    a.mean = as.numeric(), 
                    a.se = as.numeric(),
                    a.Bulk_ESS= as.numeric() ,
                    a.Tail_ESS= as.numeric() ,
                    a.Rhat= as.numeric(),
                    
                    B.mean = as.numeric(),
                    B.se = as.numeric(),
                    B.Bulk_ESS= as.numeric() ,
                    B.Tail_ESS= as.numeric() ,
                    B.Rhat= as.numeric(),
                    
                    Y.mean = as.numeric(),
                    Y.se= as.numeric(),
                    Y.Bulk_ESS= as.numeric() ,
                    Y.Tail_ESS= as.numeric() ,
                    Y.Rhat= as.numeric(),
                    
                    Z.mean = as.numeric(),
                    Z.se= as.numeric(),
                    Z.Bulk_ESS= as.numeric() ,
                    Z.Tail_ESS= as.numeric() ,
                    Z.Rhat= as.numeric(),
                    samples= as.numeric())



for ( i in unique(df$idx)){
  print(i)
   
  # Subset the file:
  df.sub <- df %>% filter(idx == i, PAR > 0) %>% na.omit

  # get priors:
   #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())
  
  try(model.brms <- brm( bf( equation, a+B+Y+Z ~ 1, nl=TRUE),
                 prior = priors , data = df, iter = iterations, cores =3, chains = 1, backend = "cmdstanr"), silent= F)
 
  print(model.brms)
  
  try(model.brms.df <- summary(model.brms)$fixed , silent = T)
  try(model.brms.df.a <- model.brms.df %>% filter( row.names(model.brms.df) == 'a_Intercept'), silent = F)
  try(model.brms.df.B <- model.brms.df %>% filter( row.names(model.brms.df) == 'B_Intercept'), silent = F)
  try(model.brms.df.Y <- model.brms.df %>% filter( row.names(model.brms.df) == 'Y_Intercept'), silent = F)
  try(model.brms.df.Z <- model.brms.df %>% filter( row.names(model.brms.df) == 'Z_Intercept'), silent = F)
  
  samples <- data.frame %>% filter(YearMon == i)%>% select(nee)  %>% na.omit %>% nrow
  
  baseline <- as.Date(paste(i, '-01', sep="")) %>% lubridate::days_in_month() *48 %>% as.numeric
  
  
  try(results <- data.frame( idx = i, 
                         a.mean = model.brms.df.a$Estimate ,
                         a.se = model.brms.df.a$Est.Error ,
                         a.Bulk_ESS = model.brms.df.a$Bulk_ESS , 
                         a.Tail_ESS = model.brms.df.a$Tail_ESS ,
                         a.Rhat = model.brms.df.a$Rhat ,
                         
                         B.mean = model.brms.df.B$Estimate ,
                         B.se = model.brms.df.B$Est.Error ,
                         B.Bulk_ESS = model.brms.df.B$Bulk_ESS , 
                         B.Tail_ESS = model.brms.df.B$Tail_ESS ,
                         B.Rhat = model.brms.df.B$Rhat ,
                         
                         Y.mean = model.brms.df.Y$Estimate ,
                         Y.se = model.brms.df.Y$Est.Error ,
                         Y.Bulk_ESS = model.brms.df.Y$Bulk_ESS , 
                         Y.Tail_ESS = model.brms.df.Y$Tail_ESS ,
                         Y.Rhat = model.brms.df.Y$Rhat,
                         
                         Z.mean = model.brms.df.Z$Estimate ,
                         Z.se = model.brms.df.Z$Est.Error ,
                         Z.Bulk_ESS = model.brms.df.Z$Bulk_ESS , 
                         Z.Tail_ESS = model.brms.df.Z$Tail_ESS ,
                         Z.Rhat = model.brms.df.Z$Rhat,
                         samples= samples/baseline *100 ), silent = T)
  
 print(results)
 
  try(parms <- parms %>% rbind(results), silent = T)
}
 
return(parms ) 
}

# EOF

