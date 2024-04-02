# This is a script that computes generalized mixed models to describe and predict post-surgery cognitive decline for
# the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

s = 87542 # seed for reproducibility

# rstan options
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores
ch = 4 # number of chains
it = 2500 # iterations per chain
wu = 500 # warm-up iterations
ad = .99 # adapt_delta parameter

library(here) # directory management
library(tidyverse) # data wrangling
library(brms) # model fitting
library(cmdstanr) # model fitting
library(tidybayes) # posterior manipulation

# create folders to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("mods","figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )


# DATA READ ----

# original data set
d0 <-
  read.csv( here("_data","20220508_dbs_longCOG_data.csv") , sep = "," ) %>%
  filter( included == 1 ) %>% # only STN-DBS treated patients with pre- and post-surgery data
  filter( complete.cases(drs_tot) ) # get rid of three dummy rows due to more than one stimulation parameter (no DRS-2)

# EFA models including regression-based factor scores
for( i in names( readRDS( here("mods","factanal.rds") ) ) ) assign( i , readRDS( here("mods","factanal.rds") )[[i]] )

# imputed data sets
d.imp <- lapply( 1:imp, function(i) read.csv( here( "_data", "imputed", paste0("imputed_df_",i,".csv") ), sep = "," ) )


# PRE-PROCESSING ----

# extract means and SDs for scaling
scl <-
  
  lapply(
    setNames( c("mean","sd"), c("M","SD") ),
    function(i)
      
      lapply(
        setNames( c("drs_tot","bdi","ledd_mg","age_ass_y"), c("drs","bdi","led","age") ),
        function(j)
          do.call( i, list( d0[ ,j], na.rm = T ) )

      )
  )

# add median time of assessment before surgery
scl$Md <- list( time = -median( d0[d0$ass_type == "pre", ]$time_y , na.rm = T ) )

# scale the data for description
df <-
  
  d0 %>%
  # scale all DRS-2 and time already
  mutate(
    time = time_y + scl$Md$time,
    drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
    cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
  ) %>%
  # keep only variables of interest
  select( id, time, drs, cens_drs )

# merge longitudinal df with baseline factor scores in efa and test scores in d.imp for prediction
d1 <-
  
  lapply(
    
    1:imp,
    function(i)
      
      d0 %>%
      left_join( cbind.data.frame( id = d.imp[[i]]$id, efa[[i]][[nf-2]]$scores ) , by = "id" ) %>% # add pre-surgery factor scores
      mutate(
        time = time_y + scl$Md$time,
        drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
        bdi = ( bdi - scl$M$bdi ) / scl$SD$bdi,
        led = ( ledd_mg - scl$M$led ) / scl$SD$led,
        age = ( age_ass_y - scl$M$age ) / scl$SD$age,
        sex = as.factor( sex ), # for better estimation of BDI in the second (covariate) model
        cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ) # right censoring for DRS == 144
      ) %>%
      # keep only variables of interest
      select(
        id, time, drs, cens_drs, bdi, led, age, sex, # outcomes, demographics, clinics
        exec_fun, epis_mem, verb_wm, visp_mem, set_shift, anxiety, visp_wm # pre-surgery cognition
      ) %>%
      # add raw test scores for model comparisons
      left_join( cbind.data.frame( d.imp[[i]] ), by = "id" )

  )

# extract test column names
tests <- names(d.imp[[1]])[-1]

# loop through all imputations to get means and SDs of pre-surgery cognitive domains and tests scores,
# then transform the pre-surgery cognition in each data set to (pre-surgery) zero mean, unit SD variables
for ( i in 1:imp ) {
  
  # start with pre-processing the cognitive domains
  for ( j in doms ) {
    
    # calculate scaling values
    scl$M[[j]][[i]] <- efa[[i]][[nf-2]]$scores[ ,j] %>% mean()
    scl$SD[[j]][[i]] <- efa[[i]][[nf-2]]$scores[ ,j] %>% sd()
    
    # scale in the jth imputed data set
    d1[[i]][[j]] <- case_when(
      # all but anxiety measures will be inverse such that parameters
      # can be interpreted as effect of deficit in said measure
      j == "anxiety" ~ ( d1[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]],
      j != "anxiety" ~ ( scl$M[[j]][[i]] - d1[[i]][[j]] ) / scl$SD[[j]][[i]]
    )
  }
  
  # next pre-process single cognitive tests
  for ( j in tests ) {
    
    # calculate scaling values
    scl$M[[j]][[i]] <- d.imp[[i]][,j] %>% mean()
    scl$SD[[j]][[i]] <- d.imp[[i]][,j] %>% sd()
    
    # scale in the jth imputed data set
    if ( j %in% c( paste0("sc_tmt_",c("a","b")), paste0("sc_pst_",c("d","w","c")), paste0("sc_staix",1:2) ) ) {

      # all but reaction speed and anxiety measures will be inversed such that parameters
      # can be interpreted as effect of deficit in said measure
      d1[[i]][[j]] <- ( d1[[i]][[j]] - scl$M[[j]][[i]] ) / scl$SD[[j]][[i]]
    } else d1[[i]][[j]] <- ( scl$M[[j]][[i]] - d1[[i]][[j]] ) / scl$SD[[j]][[i]]

  }
}

# save the data and scaling values
saveRDS( list(d0 = d0, df = d1, scl = scl, tests = tests, doms = doms, imp = imp ), here("_data","longitudinal_df.rds") )


# DESCRIPTIVE MODELS ----

# set-up the linear model
f0 <-
  
  list(
    m0_linear = bf( drs | cens(cens_drs) ~ 1 + time + (1 + time | id) ),
    m0_spline = bf( drs | cens(cens_drs) ~ t2(time) + (1 + time | id) )
  )

# set-up priors (using brms default non-informative to allow for information from data to prevail)
p0 <- NULL

# model fitting
m0 <-

  lapply(

    setNames( names(f0), names(f0) ),
    function(i)
      brm( formula = f0[[i]], family = student(), prior = p0,
           data = df, sample_prior = T, seed = s, chains = ch,
           iter = it, warmup = wu, control = list( adapt_delta = ad ),
           file = here( "mods",paste0(i,".rds") ), save_model = here( "mods", paste0(i,".stan") )
           )

  )


## POSTERIOR PREDICTIONS ----

if ( !file.exists( here("_data","ppred.rds") ) ) {
  
  # simulate values for each patient each six months from 2 years before to 12 years after surgery
  n_seq = 24
  
  # re-format id in d to a factor
  d0$id <- factor( d0$id, levels = unique(d0$id) )
  
  # prepare data sets for prediction for each subject in the data set
  d_seq <-
    
    # prepare all the columns
    expand.grid( seq( from = -2, to = 12, length.out = n_seq ), levels(d0$id) ) %>%
    `colnames<-` ( c("time_y","id") ) %>%
    mutate( time = time_y + scl$Md$time ) %>%
    mutate( !!!setNames(rep( NA, length( c(doms,tests) ) ), c(doms,tests) ) ) %>% # add empty columns for predictors
    
    # add subjects' median (w.r.t. imputations) cognitive profile
    mutate(
      across(
        all_of( c(doms,tests) ),
        ~ sapply(
          1:length(.x),
          function(i)
            sapply( 1:imp, function(j) d1[[j]][ d1[[j]]$id == id[i] & d1[[j]]$time < ( 0 + scl$Md$time ) , cur_column() ] ) %>%
            median()
        )
      )
    )
  
  # calculate posterior predictions
  ppred <- list()
  
  # loop through ids and add predictions
  for ( j in unique(d_seq$id) ) {
    
    ppred[[names(m0)[1]]][[j]] <-
      
      d_seq[ d_seq$id == j, ] %>%
      add_predicted_draws( m0[[1]], seed = s ) %>%
      mutate(drs = scl$M$drs + scl$SD$drs * .prediction) %>%
      median_hdci( .width = .95 ) %>% # need HDCI instead of HDI to deal with multimodality
      mutate(drs.upper = ifelse(drs.upper > 144, 144, drs.upper) )
    
  }
  
} else {
  
  ppred <- readRDS( here("_data","ppred.rds") )
  
}


# PREDICTIVE MODELS ----

# set-up the linear models
f1 <-
  
  list(
    m1_lasso_doms = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf(),
    m2_lasso_tests = paste0( "drs | cens(cens_drs) ~ 1 + ", paste("time", tests, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf()
  )

# set-up priors (using the same priors for both models)
p1 <-
  
  c(
    # fixed effects
    prior( normal(0.3, .1), class = Intercept ),
    prior( lasso(1), class = b ),
    # random effects
    prior( normal(0, .1), class = sd, coef = Intercept, group = id ),
    prior( normal(0, .1), class = sd, coef = time, group = id ),
    prior( lkj(2), class = cor ),
    # other distributional parameters
    prior( exponential(1), class = sigma ),
    prior( gamma(2, 0.1), class = nu )
  )

# model fitting
m1 <-
  
  lapply(
    
    setNames( names(f1), names(f1) ),
    function(i)
      brm_multiple( formula = f1[[i]], family = student(), prior = p1,
                    data = d1, sample_prior = T, seed = s, chains = ch,
                    iter = it, warmup = wu, control = list( adapt_delta = ad ),
                    file = here( "mods",paste0(i,".rds") ), save_model = here( "mods", paste0(i,".stan") )
                    )
    
  )

# clean the environment
rm( list = ls()[ !( ls() %in% c("d0","d1","d_seq","imp","m1","scl","doms","tests","ch","it","wu","ad","s","ppred") ) ] )
gc()

# compute PSIS-LOO for each imputation in the primary models
# first read loo if there is already something computed
if ( file.exists( here("mods","lasso_psis_loo.rds") ) ) l <- readRDS( here("mods","lasso_psis_loo.rds" ) )

# if no PSIS-LOO was already computed, start from a scratch
if ( !exists("l") ) {
  
  l <- list()
  for ( i in names(m1) ) {
    for ( j in 1:imp ) {

      l[[i]][[j]] <- loo( m1[[i]] , newdata = d1[[j]] )
      saveRDS( l, here("mods","lasso_psis_loo.rds") ) # save after each iteration
      print( paste0(i,", dataset #",j) ) # print the number of the last data set with PSIS-LOO

    }
  }
  
  # otherwise continue from the next data set after the last one with already computed PSIS-LOO
} else {
  
  for ( i in names(m1) ) {
    for ( j in 1:(imp-1) ) {
      
      j = length(l[[i]])
      if ( j < imp ) {
        
        l[[i]][[j+1]] <- loo( m1[[i]] , newdata = d1[[j+1]] )
        saveRDS( l, here("mods","lasso_psis_loo.rds") )
        print(paste0(i,", dataset #",j+1) ) # print the number of the last data set with PSIS-LOO

      }
    }
  }
}


## POSTERIOR PREDICTIONS ----

# run the code only if they were not predicted yet to save time
if( !file.exists( here("_data","ppred.csv") ) ) {
  
  if( !file.exists( here("_data","ppred.rds") ) ) {
    
    # prepare a list for predictions stratified by subjects chunks
    # calculating prediction of single data points based on fixed-effects, random-effects and remaining distributional (residual) parameters
    #ppred <- list()
    gc()
    
    # add predictions to ppred
    for (i in names(m1) ) {
      for ( j in unique(d_seq$id) ) {
        
        ppred[[i]][[j]] <-
          
          d_seq[ d_seq$id == j, ] %>%
          add_predicted_draws( m1[[i]], seed = s ) %>%
          mutate(drs = scl$M$drs + scl$SD$drs * .prediction) %>%
          median_hdi( .width = .95 ) %>%
          mutate(drs.upper = ifelse(drs.upper > 144, 144, drs.upper) ) # manual censoring
        
      }
    } # took 6558.316 sec in total
    
    # save the original prediction file as .rds
    saveRDS( ppred, here("_data","ppred.rds") )
    
  } else {
    
    ppred <- readRDS( here("_data","ppred.rds") )
    
  }
  
  # prepare a table of posterior predictions
  ppred <-
    
    lapply(
      
      names(ppred),
      function(i)
        do.call( rbind.data.frame, ppred[[i]] ) %>%
        add_column(
          `Predicted by:` =
            case_when(
              i == "m0_linear" ~ "time only",
              i == "m1_lasso_doms" ~ "factor scores",
              i == "m2_lasso_tests" ~ "test scores"
            )
        )
      
    ) %>%
    
    do.call( rbind.data.frame, . ) # collapse to a single file
  
  # save as .csv
  write.table( ppred, here( "_data","ppred.csv"), sep = ",", row.names = F, quote = F )
  
} else {
  
  ppred <- read.csv( here("_data","ppred.csv"), sep = "," )

}


# ROBUSTNESS CHECK MODELS ----

# clean the environment
rm( list = ls()[ !( ls() %in% c("d0","d1","scl","doms","tests","ch","it","wu","ad","s") ) ] )
gc()

# set-up the linear models for drs
f.drs <-
  
  list(
    m3_doms_cov = ( paste0( "drs | cens(cens_drs) ~ 1 + ", paste( paste0( "time * ", c("age","mi(bdi)","mi(led)") ), collapse = " + " ), " + ", paste("time", doms, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf() ) + student(),
    m4_tests_cov = ( paste0( "drs | cens(cens_drs) ~ 1 + ", paste( paste0( "time * ", c("age","mi(bdi)","mi(led)") ), collapse = " + " ), " + ", paste("time", tests, sep = " * ", collapse = " + "), " + (1 + time | id)" ) %>% as.formula() %>% bf() ) + student()
  )

# set-up linear models for covariates with missing values
f.bdi <- bf( bdi | mi() ~ 1 + time + sex + age + mi(led) + (1 + time | id) ) + gaussian()
f.led <- bf( led | mi() ~ t2(time) + (1 | id) ) + gaussian()

# set-up priors (using the same priors for both models)
p <-
  
  c(
    # DRS-2
    prior( normal(0.3, .1), class = Intercept, resp = drs),
    prior( lasso(1), class = b, resp = drs ),
    prior( normal(0, .1), class = sd, coef = Intercept, group = id, resp = drs ),
    prior( normal(0, .1), class = sd, coef = time, group = id, resp = drs ),
    prior( exponential(1), class = sigma, resp = drs ),
    prior( gamma(2, 0.1), class = nu, resp = drs ),
    # BDI-II
    prior( normal(.6, .5), class = Intercept, resp = bdi ),
    prior( normal(0, .5), class = b, coef = time, resp = bdi ),
    prior( normal(0, .5), class = b, coef = sexmale, resp = bdi ),
    prior( normal(0, .5), class = b, coef = miled, resp = bdi ),
    prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = bdi ),
    prior( normal(0, .5), class = sd, coef = time, group = id, resp = bdi ),
    prior( exponential(1), class = sigma, resp = bdi ),
    # LEDD
    prior( normal(0, 100), class = Intercept, resp = led ),
    prior( normal(0, 100), class = b, coef = t2time_1, resp = led ),
    prior( normal(0, .5), class = sd, coef = Intercept, group = id, resp = led ),
    prior( exponential(1), class = sigma, resp = led ),
    # covariance structure
    prior( lkj(2), class = cor)
  )

# model fitting
m <-
  
  lapply(
    setNames( names(f.drs), names(f.drs) ),
    function(i)
      brm_multiple( formula = f.drs[[i]] + f.bdi + f.led, prior = p,
                    data = df, sample_prior = T, seed = s, chains = ch,
                    iter = it, warmup = wu, control = list( adapt_delta = ad ),
                    file = here( "mods", paste0(i,".rds") ),
                    save_model = here("mods", pastee0(i,".stan") )
                    )
  )

# print the highest Rhat to get an idea whether chains converged
sapply(
  names(m),
  function(i)
    sapply(
      paste0("_",c("drs","bdi","led"),"_" ),
      function(j)
        max( m[[i]]$rhats %>% select( contains(j) ), na.rm = T ) )
)


# STIMULATION/STN INTERSECTIONS ----

# clean the environment
rm( list = ls()[ !( ls() %in% c("scl","s") ) ] )
gc()

# new rstan options
ch = 4 # number of chains
it = 2000 # iterations per chain
wu = 1000 # warm-up iterations
ad = .90 # adapt_delta parameter

# read data & variable names
v <- read.csv( here("_data","var_nms.csv"), sep = ";" )[1:35, ] # variables to be inspected
d0 <- read.csv( here("_data","20240330_data_vat_plus.csv"), sep = "," ) # full data-set
d1 <- d0 %>% filter( ass_type == "pre" ) # baseline data

# pre-process the data
df <-
  
  d0 %>%
  mutate(
    time = time_y + scl$Md$time,
    drs = ( drs_tot - scl$M$drs ) / scl$SD$drs,
    cens_drs = ifelse( drs == max(drs, na.rm = T) , "right" , "none" ),
    elloc = as.factor(elloc),
    post = as.numeric(post)
  ) %>%
  select( id, time, drs, cens_drs, elloc, post, starts_with("STN") )


## PRE-SURGERY PROFILE ----

# printing rounded number
rprint <- function( x, dec=2 ) sprintf( paste0("%.",dec,"f"), round( x , dec) )

# variables to be log-transformed
logv <- c( paste0("tmt_", c("a","b") ), paste0("pst_", c("d","w","c") ) )

# fit models and summarise
tab <-
  
  sapply(
    
    v$variable[-4], # do not use sex in this model for continuous responses
    function(i)
      
      with(
        
        d1, {
          
          if ( i %in% logv ) y <- log( get(i) ) else y <- get(i) # log transform if applicable
          
          # extract outcomes
          y0 <- c( na.omit( y[elloc == 0] ) )
          y1 <- c( na.omit( y[elloc == 1] ) )
          
          # extract scaling values
          M <- mean( y, na.rm = T )
          SD <- sd( y, na.rm = T )
          
          # data for stan model
          dlist <- list( y0 = ( (y0-M)/SD ), y1 = ( (y1-M)/SD ), N0 = length(y0), N1 = length(y1) )
          
          # fit the model
          mod <- cmdstan_model( here("mods","m5_ttest.stan") )
          fit <- mod$sample( data = dlist, seed = s, chains = ch, iter_warmup = wu , iter_sampling = it, adapt_delta = ad ) # fit it
          
          # posteriors draws
          post <- fit$draws( format = "data.frame" )
          post <- apply( post[ ,2:5], 2, function(x) x * SD + M )
          if( i %in% logv ) post <- exp(post)
          
          # posterior differences
          diff <- cbind( mu = post[ ,"mu_1"] - post[ ,"mu_0"], sigma = post[ ,"sigma_1"] - post[ ,"sigma_0"] )
          
          # summaries  of posterior differences
          sum <- apply( diff, 2, function(x) paste0( rprint( median(x) ), " [", rprint( hdi(x)[1] ),", ", rprint( hdi(x)[2]), "]" ) )
          prob <- apply( diff, 2, function(x) sum( x > 0 ) / length(x) )
          
          # sample description
          n <- with( dlist, paste( N0, N1, sep = "/" ) )
          desc <- sapply( 0:1, function(j) paste0( rprint( mean( get(i)[elloc == j], na.rm = T ) ), " ± ", rprint( sd( get(i)[elloc == j], na.rm = T ) ) ) )
          
          # print results
          return(
            c( n = n,
               `0` = desc[1], `1` = desc[2],
               mu_sum = sum[["mu"]], mu_prob = prob[["mu"]],
               sigma_sum = sum[["sigma"]], sigma_prob = prob[["sigma"]]
               )
          )
          
        }
      )

  ) %>%
  
  t() %>%
  as.data.frame() %>%
  rownames_to_column("var") %>%
  
  # add row for sex
  add_row(
    var = "sex",
    n = paste( nrow( d1[ d1$elloc == 0, ] ), nrow( d1[ d1$elloc == 1, ] ) ) ,
    `0` = paste0( nrow( d1[ d1$elloc == 0 & d1$sex == "male", ] ), " (", rprint( 100 * nrow( d1[ d1$elloc == 0 & d1$sex == "male", ] ) / nrow(d1), 1 ) , "%)" ),
    `1` = paste0( nrow( d1[ d1$elloc == 1 & d1$sex == "male", ] ), " (", rprint( 100 * nrow( d1[ d1$elloc == 1 & d1$sex == "male", ] ) / nrow(d1), 1 ) , "%)" ),
    mu_sum = NA,
    mu_prob = NA,
    sigma_sum = NA,
    sigma_prob = NA,
    .after = 3
  ) %>%
  
  # finishing touches
  mutate( var = sapply( var, function(i) v[ v$variable == i, "name" ] ) ) %>%
  mutate( across( all_of( ends_with("prob") ), ~ as.numeric(.x) ) )

# save it
write.table( tab, here("tabs","elloc_comaprsions.csv"), sep = ";", row.names = F, quote = F )


## DESCRIPTIVE MODEL ----

f6 <- bf( drs | cens(cens_drs) ~ 1 + time * elloc + (1 + time | id) ) # linear model
p6 <- NULL # prior
contrasts(df$elloc) <- -contr.sum(2)/2

# model fitting
m6 <- brm( formula = f6, family = student(), prior = p6,
           data = df, sample_prior = T, seed = s, chains = ch,
           iter = it, warmup = wu, control = list( adapt_delta = ad ),
           file = here("mods","m6_desc_check.rds"), save_model = here("mods","m6_desc_check.stan")
           )


## VAT/STN INTERSECTIONS MODEL ----

# extract structure mapping
struct <-

  expand.grid(
    atlas = paste0( "Atlas", c("Volume","Intersection") ),
    side = c("left","right"),
    struct = paste0( "STN", c("","_motor","_associative","_limbic") )
  ) %>%

  mutate( out = paste( struct, side, atlas, sep = "_" ) ) %>%
  pivot_wider( names_from = "atlas", values_from = "out" ) %>%
  mutate( proportion = paste( struct, side, "proportion", sep = "_" ) )

# compute proportions of VAT/STN intersections
for ( i in 1:nrow(struct) ) df[ , struct$proportion[i] ] <- df[ , struct$AtlasIntersection[i] ] / df[ , struct$AtlasVolume[i] ]

# add scaling for VAT/STN intersection proportions
for ( i in names(df)[grepl( "proportion", names(df) )] ) {
  
  scl$M[[i]] <- mean( df[ df$post == 0, i ], na.rm = T )
  scl$SD[[i]] <- sd( df[ df$post == 0, i ], na.rm = T )
  
}

# scale the predictors
df <-

  df %>%
  filter( elloc == 1 ) %>%
  select( id, time, drs, cens_drs, post, ends_with("proportion") ) %>%
  mutate( across( all_of( ends_with("proportion") ), ~ ( .x - scl$M[[cur_column()]] ) / scl$SD[[cur_column()]] ) )

# linear model
f7 <-
  
  paste0( "drs | cens(cens_drs) ~ 1 + ", paste( paste0( "time * ", struct$proportion[-c(1:2)] ), collapse = "+" ), " + (1 + time | id)" ) %>%
  as.formula() %>%
  bf()

# prior
p7 <-

  c(
    # fixed effects
    prior( normal(0, .5), class = Intercept ),
    prior( normal(0, .5), class = b ),
    # random effects
    prior( exponential(2), class = sd, coef = Intercept, group = id ),
    prior( exponential(2), class = sd, coef = time, group = id ),
    prior( lkj(2), class = cor ),
    # other distributional parameters
    prior( exponential(1), class = sigma ),
    prior( gamma(2, 0.1), class = nu )
  )

# model fitting
m7 <-
  
  brm(
    formula = f7, family = student(), prior = p7,
    data = df, sample_prior = T, seed = s, chains = ch, iter = it, warmup = wu, control = list( adapt_delta = ad ),
    file = here( "mods","m7_stn_intersect.rds"), save_model = here( "mods", "m7_stn_intersect.stan")
  )


# plot 'Bayesian p-values' to summarise differences ----
#
#tab %>%
#  
#  # extract variables needed
#  select( var, ends_with("prob") ) %>%
#  mutate(
#    var =
#      case_when(
#        var == "MDS-UPDRS III (ON medication)" ~ "MDS-UPDRS III ON",
#        var == "MDS-UPDRS III (OFF medication)" ~ "MDS-UPDRS III OFF",
#        var == "Disease duration at surgery (years)" ~ "Dis. dur. at surgery",
#        .default = var
#      )
#  ) %>%
#  
#  # pre-process the table
#  mutate( var = sapply( var, function(i) strsplit( i, " (", fixed = T )[[1]][1] ) ) %>%
#  mutate( var = factor( var, levels = reorder(var,mu_prob), ordered = T ) ) %>%
#  filter( var != "Sex" ) %>%
#  pivot_longer( cols = -var, names_to = "Parameter", values_to = "prob" ) %>%
#  mutate(
#    pdir = ifelse( prob < .5, 1-prob, prob ),
#    `Directionality: ` = factor( ifelse( prob < .5, "included < excluded", "excluded < included" ) ),
#    p = 2 * (1 - pdir)
#  ) %>%
#  
#  # plot it
#  ggplot() +
#  aes( x = p, y = reorder( var, p, decreasing = T ), fill = `Directionality: ` ) +
#  geom_bar( stat = "identity" ) +
#  scale_fill_manual( values = c("navy","skyblue") ) +
#  geom_vline( xintercept = .05, linetype = "dashed", linewidth = 1.2, colour = "red" ) +
#  facet_wrap( ~ Parameter, ncol = 2, labeller = as_labeller( c( mu_prob = "mu", sigma_prob = "sigma" ), label_parsed ) ) +
#  labs( y = NULL, x = "Bayesian p-value" ) +
#  theme_minimal( base_size = 14 ) +
#  theme( legend.position = "bottom" )
