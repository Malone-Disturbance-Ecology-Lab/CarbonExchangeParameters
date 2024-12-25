#' fit a Light Response Curves (LRC) by an index to get a parameter file
#'
#' `LRC_PARMS` returns a dataframe of parameter values by the index used to fit them.
#' `LRC_PARMS` requires a few arguments @data.frame,  @iterations, and @priors.
#' `dataframe` is a dataframe that containes nee, an index, and PAR.
#' `iterations` is the number of iterations to run brms. 
#' `priors ` is the priors for the brms to use.

rm(list=ls())

library(brms) 
library(cmdstanr)
library(ggplot2)
library(beepr)
library(tidybayes)
library(tidyverse)
library(ggpubr)

LRC_PARMS <- function( data.frame, iterations, priors){
  
equation <- nee ~ (a1 * PAR * ax)/(a1 * PAR + ax) + r

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
                    r.Rhat= as.numeric())

for ( i in unique(data.frame$idx)){
  print(i)
   
  # Subset the file:
  df <- data.frame %>% filter(idx == i)
  # get priors:
 
  # priors <- get_prior(bf(equation, a1+ax+r ~ 1, nl=TRUE),
  #                    data = df %>% filter(PAR > 0),
  #                    family = poisson())
  model.brms <- brm( bf( equation, a1+ax+r ~ 1, nl=TRUE),
                 prior = priors , data = df, 
                 backend = "cmdstanr", iter = iterations, cores =6)
  
  model.brms.df <- summary(model.brms)$fixed 
  
  model.brms.df.a1 <- model.brms.df %>% filter( row.names(model.brms.df) == 'a1_Intercept')
  model.brms.df.ax <- model.brms.df %>% filter( row.names(model.brms.df) == 'ax_Intercept')
  model.brms.df.r <- model.brms.df %>% filter( row.names(model.brms.df) == 'r_Intercept')
  
  results <- data.frame( idx = i, 
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
                         r.Rhat = model.brms.df.r$Rhat )
  parms <- parms %>% rbind(results)
}
 
return(parms ) 
}

# Example of priors: 

priors <-  prior(normal(-0.01, 0.1), nlpar = "a1", lb=-0.2, ub= 0) +
prior(normal( -7.65 ,  0.33), nlpar = "ax", lb=-8, ub= -5) +
  prior(normal(2.10, 0.11), nlpar = "r", lb=1.9, ub= 2.2)

# EOF

