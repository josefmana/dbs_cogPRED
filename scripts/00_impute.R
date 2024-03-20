# This is a script computes imputed data sets of pre-surgery battery for the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

imp = 100 # number of multiple imputations to account for missing pre-surgery data
s = 87542 # seed for reproducibility

library(here) # directory management
library(tidyverse) # data wrangling
library(missMDA) # imputations

# create folders to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( "_data/imputed", function(i) if( !dir.exists(i) ) dir.create(i) )


# READ ----

d <-

  read.csv( here("_data","20220508_dbs_longCOG_data.csv"), sep = "," ) %>%
  filter( included == 1 & ass_type == "pre" ) %>%

  # transform data as required
  select( c( 2, which(names(.) == "tmt_a"):which(names(.) == "fp_dr"), which(names(.) %in% paste0("staix",1:2)) ) ) %>% # keep test only
  rename_with( ~paste0("sc_",.), -!!"id" ) %>% # change names such that they get correct label in post-processing
  mutate( across( c( paste0("sc_tmt_", c("a","b")), paste0("sc_pst_", c("d","w","c")) ), log ) ) %>% # log-transform reaction times
  mutate( across( !starts_with("id"), ~ as.vector( scale( .x, center = T, scale = T ) ) ) ) # standardise (center and scale)


# IMPUTE ----

nb <- estim_ncpPCA( d[ ,-1] , ncp.min = 0, ncp.max = 10 , nbsim = imp ) # find out the optimal number of components

set.seed(s) # seed for reproducibility
d.imp <- MIPCA( d[ ,-1] , ncp = nb$ncp , nboot = imp ) # impute via PCA-based multiple imputation


# SAVE ----

for ( i in 1:imp ) write.table( x = cbind( id = d$id, d.imp$res.MI[[i]] ),
                                file = here( "_data", "imputed", paste0("imputed_df_",i,".csv") ),
                                sep = ",",
                                row.names = F,
                                quote = F
                                )


# RENV UPDATE ----

renv::snapshot()
