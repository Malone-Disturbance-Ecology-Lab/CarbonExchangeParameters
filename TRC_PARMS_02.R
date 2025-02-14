library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)

# this function is very sensitive to priors.
# Example of priors: 
priors.trc <-  prior(normal(1, 5), nlpar = "r0", lb=0, ub= 5) +
  prior(normal(0.07, 0.1), nlpar = "alpha", lb=0, ub= 0.2) +
  prior(normal(0, 0.05), nlpar = "beta", lb=-0.1, ub= 0.1)

TRC_PARMS_02 <- function( data.frame, iterations, priors.trc, idx, nee, PAR, TA){
  
  df <- data.frame %>% mutate(idx= idx,
                              nee= nee,
                              PAR=PAR,
                              TA= TA ) %>% select(idx, nee, PAR, TA)
  
  equation <- nee ~ r0 * exp(alpha * TA + beta * TA * TA)
  
  if( c("idx") %in% names(df) ) {
    print("idx present")
  } else{ print("The dataframe must include: idx,  nee, PAR, and TA")}
  
  if( c("nee") %in% names(df) ) {
    print("nee data present")
  }  else{ print("The dataframe must include: idx,  nee, PAR, and TA")}
  
  if( c("PAR") %in% names(df) ) {
    print("PAR data present")
  } else{ print("The dataframe must include: idx,  nee, PAR, and TA")}

  if( c("TA") %in% names(df) ) {
    print("TA data present")
  } else{ print("The dataframe must include: idx,  nee, PAR, and TA")}
  
  
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
  
  
  
  for ( i in unique(df$idx)){
    print(i)
    
    # Subset the file:
    df.sub <- df %>% filter(idx == i, PAR <= 0) %>% na.omit

    # get priors:
    #priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),data = df %>% filter(PAR > 0), family = poisson())
    
    try(model.brms <- brm( bf( equation, r0+alpha+beta ~ 1, nl=TRUE), prior = priors.trc, 
                           data = df.sub, iter = iterations, cores =3, chains = 1, backend = "cmdstanr"), silent= F)
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
    
    print(results)
    
    try(parms <- parms %>% rbind(results), silent = T)
  }
  
  return(parms ) 
}
