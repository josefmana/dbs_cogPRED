# This is a script that calculates exploratory factor analysis (EFA) of pre-surgery battery for the longitudinal cognition
# in DBS study

rm( list = ls() ) # clear environment

imp = 100 # number of multiple imputations to account for missing pre-surgery data

library(here) # directory management
library(tidyverse) # data wrangling
library(psych) # factor analysis

# create folders to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( "mods", function(i) if( !dir.exists(i) ) dir.create(i) )


# READ ----

d0 <- lapply( 1:imp, function(i) read.csv( here( "_data", "imputed", paste0("imputed_df_",i,".csv") ), sep = "," ) )


# CALCULATE EFA ----

# loop through all imputed data sets
efa <-
  
  lapply(
    
    1:imp,
    function(i)
      
      # loop through three to eight latent factors for each imputation
      lapply(
        
        3:8,
        function(j)
          
          fa( d0[[i]][ ,-1], # one-by-one use each imputed data 
              nfactors = j, # fit 3-8 factor solutions
              rotate = "varimax", # rotate varimax to enforce orthogonality and for interpretation purposes
              scores = "regression" # compute regression scores for each patient
              )
        
      )
      
  )


# POST-PROCESSING ----

# one-by-one inspect all six-factor and seven-factor solutions' loading matrices
# create a convenience function so that I don't go crazy immediately
printload <- function( i, c = .4, f = 7 ) print( efa[[i]][[f-2]]$loadings, cutoff = c, sort = T )

# choosing 7-factor solution due to good performance indexes,
# and theoretically sound loading patterns across imputed data sets
nf = 7

# prepare an array for labels of the seven-factor solution factors
domsum <- array( data = NA, dim = c(2, nf, imp), dimnames = list( c("nms","sgn"), paste0("F", 1:nf), 1:imp ) )

# read the table with seven-factor labels
domsum["nms", , ] <- t( read.csv( here("_data","efa_labels.csv"), sep = "," , row.names = 1, header = T) )

# fill-in signs of each factor in each imputation to know which scores should be reversed
domsum["sgn", , ] <-
  
  apply( domsum["nms", , ] , 2 , function(x) startsWith( x , "-") ) %>%
  t() %>%
  as.data.frame() %>%
  mutate( across( everything() , ~ ifelse( .x, -1 , 1 ) ) ) %>%
  t()

# get rid of the minus sign in labels table
domsum["nms", , ] <-
  
  domsum["nms", , ] %>%
  t() %>%
  as.data.frame() %>%
  mutate( across( everything() , ~ gsub( "-" , "" , .x ) ) ) %>%
  t()

# list all the domains
doms <-
  
  c("exec_fun", # loaded on primarily by PST, the first factor in 82% data sets
    "epis_mem", # loaded on primarily by RAVLT, the second factor in 79% data sets
    "verb_wm", # loaded on primarily by DS, the third factor in 62% data sets
    "visp_mem", # loaded on primarily by FP, the fourth factor in 45% data sets
    "set_shift", # loaded on primarily by TMT and RAVLT-B, the fifth factor in 28% data sets
    "anxiety", # loaded on primarily by STAI, the sixth factor in 60% data sets
    "visp_wm" # loaded on primarily by SS, the seventh factor in 49% data sets
    )

# switch signs where appropriate in EFA loadings and scores, and rename and sort columns
for ( i in 1:imp ) {
  
  for ( j in c("loadings","scores","Vaccounted") ) {
    
    # multiply by a diagonal matrix of 1 and -1
    if ( j %in% c("loadings","scores") ) efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]] %*% diag(domsum["sgn", , i ] )
    
    # rename the columns
    colnames( efa[[i]][[5]][[j]] ) <- domsum["nms", , i ]
    
    # reorder the columns such that they are in the same order for each imputation
    efa[[i]][[nf-2]][[j]] <- efa[[i]][[nf-2]][[j]][, doms]
    
  }
}


# SAVE RESULTS ----

saveRDS( object = list( efa = efa, imp = imp, nf = nf, doms = doms ), file = here("mods","factanal.rds") )


# RENV UPDATE ----

renv::snapshot()
