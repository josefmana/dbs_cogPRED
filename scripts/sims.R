# This is a script to check false error rates in the best case scenario for choosing significant predictors of a longitudinal
# outcome from a number of potential predictors via two-stage procedure that first pre-selects potential predictors via
# univariate regression and then continues to multiple regressions compared to Bayesian Lasso.

rm( list = ls() ) # clear environment

s = 87542 # seed for reproducibility

library(here) # directory management
library(tidyverse) # data wrangling
library(MASS) # multivariate normal distribution
library(conflicted) # for conflicted functions
library(lmerTest) # two-step procedure mixed-effects models
library(brms) # Bayesian lasso mixed-effects models
library(bayestestR) # probability of direction

conflict_prefer( name = "select", winner = "dplyr" )

options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores for computing
struct <- read.csv( here("_data","struct.csv"), sep = "," ) # read the structure of our data set


# FUNCTIONS ----

## DATA-GENERATING ----

# function on varying number of subjects, number of predictors, varying pre-selection threshold and effect sizes
# that returns synthetic dato be used for ensuing analysis; takes the following inputs:
# b0 = population-level pre-surgery intercept
# b1 = population-level annual cognitive decline
# eff = list of effect sizes (vector)
# cor = covariance matrix (identical with correlation in case of standard normals)
# struct = longitudinal data structure
sim <-
  
  function( b0 = 0, b1 = -.3, eff = rep(0,7), cor = diag( rep(1,7) ), struct = struct ) {
    
    # extract number of patients and their IDs
    ID <- unique(struct$id)
    N <- length(ID)
    
    # extract number of predictors
    n_prds <- length(eff)
    
    # simulate patient-specific predictors/variation
    t <- mvrnorm( N, rep(0,n_prds), cor ) %>% `rownames<-`( ID ) # "n_prds" predictors from multivariate standard normal distribution
    u <- mvrnorm( N, c(0,0), diag( c(.05^2,.05^2) ) ) %>% `row.names<-`( ID ) # patient-specific intercepts and slopes
    
    # prepare a data set combining pre-defined data structure (struct) and generative regression parameters
    d <-
      struct %>%
      mutate(
        b0 = b0, b1 = b1, # population-level slope
        c0 = sapply( 1:nrow(struct), function(i) u[ struct$id[i], 1 ] ), # patient-specific intercept
        c1 = sapply( 1:nrow(struct), function(i) u[ struct$id[i], 2 ] ) # patient-specific slope
      )
    
    # add patients' cognitive profile
    for ( i in 1:ncol(t) ) d[ , paste0("p",i) ] <- sapply( 1:nrow(d), function(j) t[ d$id[j], i ] )
    
    # generate the outcome
    if ( length(eff) == 1 ) d$o <- with( d, rnorm( N, ( b0 + c0 )  + ( b1 + c1 ) * time + d[ , "p1"] * eff[1] * time, 1 ) ) # exactly one non-zero effect
    else d$o <- with( d, rnorm( N, ( b0 + c0 )  + ( b1 + c1 ) * time + rowSums( sapply( 1:n_prds, function(i) d[ , paste0("p",i) ] * eff[i] * d$time ) ), 1 ) ) # none or more than one non-zero effects
    
    # return the data set
    return(d)
    
  }


## TWO-STEP PROCEDURE ----

# d = data set to be used
# thres = p-value used as a threshold for pre-selection in the two-step method (numeric)
twostep <-
  
  function( d = d, thres = .2 ) {
    
    # extract the number of predictors
    k <- names(d) %>% grepl("p", . ) %>% sum()
    
    # pre-select predictors
    sel <-
      sapply(
        # loop through predictors individually
        paste0( "p", 1:k ),
        function(i)
          # assign 1 for significant and 0 for non-significant predictors
          ifelse( summary( lmer( as.formula( paste0( "o ~ 1 + time * ", i, "+ (1 + time | id)") ), data = d ) )[["coefficients"]][3,5] < thres, 1, 0 )
      )
    
    # compute final multiple linear regression for the two-step procedure
    # use only pre-selected variables
    if( sum(sel) == 0 ) reg <- lmer( o ~ 1 + time + (1 + time | id ), data = d ) # intercept only if no variable was selected
    else reg <- lmer( as.formula( paste0( "o ~ 1 + ", paste("time", names(sel)[sel == 1], sep = " * ", collapse = " + "), " + (1 + time | id )" ) ), data = d )
    
    # extract p-values
    if ( sum(sel) == 0 ) p <- NULL # there ain't no p-value if no variable was pre-selected
    else p <- as.vector( summary(reg)$coefficients[ paste0("time:",names(sel)[sel == 1]) , 5 ] )
    
    # sum-up the results and print them
    res <- data.frame( var = paste0("p", 1:k), p = rep(1, k) ) # prepare ones (for those that were not tested, technically p = 1, easier to work with than NAs)
    res$p[ which(sel == 1) ] <- p
    return(res) # return the results
    
  }

## BAYESIAN LASSO ----

# d = data set to be used
lasso <-
  
  function( d = d ) {
    
    # extract the number of predictors
    k <- names(d) %>% grepl("p", . ) %>% sum()
    
    # fit the Bayesian Lasso
    reg <-
      brm(
        as.formula( paste0( "o ~ 1 + ", paste( "time", paste0("p",1:k), sep = " * ", collapse = " + " ), " + (1 + time | id)" ) ), # linear model
        family = gaussian(), prior = prior( lasso(1), class = b ), # likelihood and priors
        data = d, # observed variables
        seed = 87542 # seed for reproducibility
      )
    
    # extract probability of direction and recalculate to a p-value equivalent
    pd <- pd_to_p( with( p_direction(reg), pd[ Parameter %in% paste0("b_time:p",1:k) ] ), direction = "two-sided" )
    
    # sum-up the results and print them
    res <- data.frame( var = paste0("p", 1:k), p = pd )
    return(res) # return the results
    
  }


# COVARIANCE MATRIXES ----

# ---- set-up covariance matrixes ----

# first prepare a function for no-covariance cases
nocov_fun <- function(K = 7) return( diag( rep(1,K) ) ) # a function for no correlation with varying number of measures

# prepare no covariance matrixes for 7 and 23 predictors (equivalent to number of latent factors and single test scores in our study)
for ( i in c(7,23) ) assign( paste0("nocov",i), nocov_fun(K = i) )

# putative covariation (correlation in this case) structure of predictors in our study
# based on grouping of tests to common clusters with ~50 % shared variability (i.e., r = .7, R2 = .49) 
yocov23 <-
  
  matrix(
    
    c(
      1, .7, rep(0,21), .7, 1, rep(0,21), # State-Trait Anxiety Inventory
      rep(0,2), 1, .7, rep(0,19), rep(0,2), .7, 1, rep(0,19), # Trail Making Test
      rep(0,4), 1, rep(.7,2), rep(0,16), rep(0,4), .7, 1, .7, rep(0,16), rep(0,4), rep(.7,2), 1, rep(0,16), # Digit Span
      rep(0,7), 1, .7, rep(0,14), rep(0,7), .7, 1, rep(0,14), # Spatial Span
      rep(0,9), 1, rep(0,13), # Tower of London
      rep(0,10), 1, rep(.7,2), rep(0,10), rep(0,10), .7, 1, .7, rep(0,10), rep(0,10), rep(.7,2), 1, rep(0,10), # Stroop test
      rep(0,13), 1, .7, rep(0,8), rep(0,13), .7, 1, rep(0,8), # verbal fluency
      rep(0,15), 1, rep(0,7), # Similarities (abstraction)
      rep(0,16), 1, rep(.7,4), rep(0,2), rep(0,16), .7, 1, rep(.7,3), rep(0,2), rep(0,16), rep(.7,2), 1, rep(.7,2), rep(0,2), rep(0,16), rep(.7,3), 1, .7, rep(0,2), rep(0,16), rep(.7,4), 1, rep(0,2), # RAVLT 
      rep(0,21), 1, .7, rep(0,21), .7, 1 # Family Pictures test
    ),
    
    ncol = 23,
    byrow = T
  )

# save it
write.table( yocov23, file = here("_data","sim_covstruct.csv"), sep = ",", row.names = F, col.names = F )


# SIMULATE DATA ----

# as the code is fitting few hundreds of Stan models, it crashes every so often, so we will put some insurances in place
if ( !file.exists( here("mods","sims.rds") ) ) {
  
  # use the "sim" function to generate thousand of data sets that can be used to replicate our findings
  d0 <-
    
    lapply(
      
      # start by assigning average rate of decline
      setNames( -c(0,.3,.5), c("none","mild","moderate") ),
      function(i)
        
        # set-up number of predictors and covariance structure
        lapply(
          setNames( c("nocov7","nocov23","yocov23"), c("nocov7", "nocov23", "yocov23") ),
          function(j)
            # loop through 1000 iterations
            lapply( 1:1e2, function(k) sim( b0 = 0, b1 = i, eff = rep( 0, nrow( get(j) ) ), cor = get(j), struct = struct ) )
        )
      
    )
  
  # save it
  saveRDS( list( d0 = d0 ), here("mods","sims.rds")  )
  
} else {
  
  for ( i in names( readRDS( here("mods","sims.rds") ) ) ) assign( i, readRDS( here("mods","sims.rds") )[[i]] )
  
}


# FIT MODELS ----

## TWO-STEP PROCEDURE ----

# loop through all the data sets simulated above in this case
if ( !exists("t0") ) {
  
  t0 <-
    
    lapply(
      
      # loop through the average rates of decline
      setNames( names(d0), names(d0) ),
      function(i)
        
        # loop through covariance structures
        lapply(
          
          setNames( names(d0[[i]]), names(d0[[i]]) ),
          function(j)
            # loop through all one hundred data sets in each decline/covariance combination
            lapply( 1:length(d0[[i]][[j]]), function(k) twostep( d = d0[[i]][[j]][[k]], thres = .2 ) )
          
        )
    )
  
  # save after each iteration
  saveRDS( list( d0 = d0, l0 = l0, t0 = t0 ), here("mods","sims.rds") )

}


## BAYESIAN LASSO ----

# because R can crash due to hundreds of Stan models being fit, put in place some safeguards
# start from a scratch if no models were fitted yet
if ( !exists("l0") ) {
  
  # prepare a list for the Bayesian Lasso model
  l0 <- list()
  
  # loop through all decline slopes simulated above
  for ( i in names(d0) ) {
    
    # fit the Bayesian Lasso regression for each decline/covariance combination for the first one hundred data sets
    l0[[i]] <- lapply( setNames( names(d0[[i]]), names(d0[[i]]) ), function(j) lapply( 1:1e2, function(k) lasso( d = d0[[i]][[j]][[k]] ) ) )
    
    # save after each iteration
    saveRDS( list( d0 = d0, l0 = l0 ), here("mods","sims.rds") )
    
  }
  
  # if there are some models already, continue with the next decline slope in line
} else {
  
  # loop through all decline slopes simulated above
  for ( i in names(d0) ) {
    
    # compute only decline slopes that are not computed yet
    if ( !exists( i, where = l0 ) ) {
      
      # fit the Bayesian Lasso regression for each decline/covariance combination for the first one hundred data sets
      l0[[i]] <- lapply( setNames( names(d0[[i]]), names(d0[[i]]) ), function(j) lapply( 1:1e2, function(k) lasso( d = d0[[i]][[j]][[k]] ) ) )
      
      # save after each iteration
      saveRDS( list( d0 = d0, l0 = l0 ), here("mods","sims.rds") )
      
    }
  }
}


# EXTRACT SUMMARIES ----

s <- mget( c("l0","t0") ) # prepare a list

# collapse simulation results for Lasso and Two-step procedures
for( i in names(s) ) {
  
  # average decline
  for ( j in names( s[[i]] ) ) {
    
    # covariance structure
    for ( k in names( s[[i]][[j]] ) ) {
      
      # data set
      for ( l in 1:length( s[[i]][[j]][[k]] ) ) s[[i]][[j]][[k]][[l]] <- s[[i]][[j]][[k]][[l]] %>% mutate( meth = i, decl = j, covariance = k, dataset = l ) # collapse all data sets within covariance structures across methods and decline
      s[[i]][[j]][[k]] <- do.call( rbind.data.frame, s[[i]][[j]][[k]] )
    }
    
    # collapse data sets within decline levels across methods
    s[[i]][[j]] <- do.call( rbind.data.frame, s[[i]][[j]] )
  }
  
  # collapse data sets within methods
  s[[i]] <- do.call( rbind.data.frame, s[[i]] )
}

# summarise of false positives conditional on null hypothesis
falsepos <-
  
  do.call( rbind.data.frame, s ) %>%
  mutate( sig = ifelse(p < .05, 1, 0 ) ) %>% # flag findings significant on 5% level
  pivot_wider( id_cols = c("dataset","meth","covariance","decl"), values_from = sig, names_from = var ) %>% # pivot again such that each potential predictor has its own column with indicator of a "hit" per simulation/method pairs
  mutate( fp = rowSums( across( starts_with("p") ), na.rm = T ) ) %>% # calculate number of "hits"/false positives across variables for each simulation/method pair
  select( meth, covariance, decl, fp ) %>% # keep only variables of interest
  
  # prepare a table of frequencies of false positives per method
  table() %>%
  as.data.frame() %>%
  
  # prepare variables for plotting
  mutate(
    `Method:` = case_when( meth == "l0" ~ "Bayesian Lasso", meth == "t0" ~ "Two-step procedure"  ),
    Decline =
      factor(
        decl,
        levels = c("none","mild","moderate"),
        labels = c("No~average~decline~(b[1]==0.0)","Mild~average~decline~(b[1]==-0.3)","Moderate~average~decline~(b[1]==-0.5)"),
        ordered = T
      ),
    Covariance =
      factor(
        covariance,
        levels = c("nocov7","nocov23","yocov23"),
        labels = c("Independent~predictors~(k==7)","Independent~predictors~(k==23)","Covaried~predictors~(k==23)"),
        ordered = T
      ),
  )

# save it
write.table( falsepos, here("mods","sims_fp.csv"), sep = ",", row.names = F, quote = F )

