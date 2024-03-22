# This is a script that prepares graphical abstract figure for the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

s = 87542 # seed for reproducibility

library(here) # directory management
library(tidyverse) # data wrangling
library(tidybayes) # posterior manipulation
library(patchwork) # figures


# LOAD MATERIALS ----

for ( i in names( readRDS( here("_data","longitudinal_df.rds") ) ) ) assign( i, readRDS( here("_data","longitudinal_df.rds") )[[i]] ) # load all the data
m <- readRDS( here("mods","m1_lasso_doms.rds") ) # load the model
prds <- c("exec_fun","epis_mem") # list of predictors effect which we're going to visualize


# PRE-PROCESSING ----

# prepare data set with patients' median pre-surgery cognitive profile across imputations
d <-
  
  d0 %>%
  select( id, time_y , drs_tot ) %>%
  filter( complete.cases(drs_tot) ) %>%
  rename( "drs" = "drs_tot" ) %>%
  mutate( !!!setNames( rep( NA, length(doms) ), doms ) ) %>% # add empty columns for predictors
  mutate( across( all_of(doms), ~ apply( sapply( 1:imp, function(i) df[[i]][ , cur_column() ] ), 1, median ) ) ) %>%
  
  # add the quantile groups to the longitudinal data set
  left_join(
    # extract quantile groups for each patient according to each latent cognitive factor
    .[ .$time_y < 0 , ] %>%
      mutate( across( all_of(doms), ~ -ntile(.x,5)+6, .names = "{.col}_pent" ) ) %>%
      select( id, ends_with("pent") ),
    # joining by id
    by = "id"
  )


# POSTERIOR PREDICTION ----

# takes a lot of time to compute so load if the predictions already exist
if ( !file.exists( here("_data","graphstract_epred.rds") ) ) {
  
  # prepare a prediction dummy data set
  d_seq <-
    list(
      epred_fix = list(), # prediction of the expectation (epred) based on fixed-effects only
      epred_all = list() # prediction of the expectation based on both fixed- and random-effects
    )
  
  # number of predictions (time-points) per predictor group/quantile to calculate
  n_seq = 30
  
  # loop through predictors
  for ( i in names(d_seq) ) {
    for ( j in prds ) {
      
      d_seq[[i]][[j]] <-
        
        expand.grid(
          as.vector( by( d[[j]], d[[ paste0(j,"_pent") ]], median ) ),
          seq( from = -2, to = 12, length.out = n_seq )
        ) %>%
        
        `colnames<-` ( c(j, "time_y") ) %>%
        
        # add all variables needed
        mutate(
          id = "sub000",
          time = time_y + scl$Md$time,
          grp = rep( 1:5, n_seq ), # in the expand.grid above, need to write predictor first, time second, otherwise this is incorrect
          across( any_of(prds), ~ if( j == cur_column() ) {.x} else 0 ),
          !!!setNames( rep( NA, length( doms[ !(doms %in% j) ] ) ), doms[ !(doms %in% j) ] ),
          across( all_of( doms[ !(doms %in% j) ] ), ~ 0 )
        ) %>%
        
        # add predictions of the expectation (e_pred) to d_seq
        add_epred_draws( m , allow_new_levels = T, seed = s, if (i == "epred_fix") { re_formula = NA } ) %>%
        mutate(.epred = scl$M$drs + scl$SD$drs * .epred) %>%
        median_hdi(.width = .95)
      
    } 
  }
  
  # save it
  saveRDS( d_seq, here("_data","graphstract_epred.rds") )
  
} else {
  
  # if it was calculated in past, just load it
  d_seq <- readRDS( here("_data","graphstract_epred.rds") )
  
}


# PLOTTING ----

# plot facets' names
nms <-
  list(
    exec_fun = bquote("Longitudinal cognition predicted by pre-surgery" ~bold("executive functions") ),
    epis_mem = bquote("Longitudinal cognition predicted by pre-surgery" ~bold("episodic memory") )
  )

# prepare figure for the graphical abstract
f <-
  
  lapply(
    
    setNames(prds,prds),
    function(i)
      
      d[ , c( "time_y","id","drs") ] %>%
      mutate( grp = d[ , paste0(i, "_pent") ] ) %>%
      ggplot( aes(x = time_y, y = drs, group = id) ) +
      geom_hline( yintercept = 139, linetype = "dotted", linewidth = 1.5, alpha = 1 ) +
      geom_line( linewidth = .5, alpha = .33 ) +
      geom_point( size = 4, alpha = .33 ) +
      geom_ribbon(data = d_seq$epred_all[[i]], alpha = .1, aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = grp ) ) +
      geom_ribbon(data = d_seq$epred_fix[[i]], alpha = .3, aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = grp ) ) +
      geom_line(data = d_seq$epred_fix[[i]], linewidth = 2.5, aes( x = time_y, y = .epred, color = grp ) ) +
      scale_color_gradient( low = "red", high = "grey66" ) +
      scale_y_continuous(name = "DRS-2", limits = c(80,153), breaks = seq(80,140,20), labels = seq(80,140,20) ) +
      scale_x_continuous(name = "Time from surgery (years)", limits = c(-3,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
      facet_wrap(
        ~  grp,
        nrow = 1,
        labeller = as_labeller( c(`1` = "first pentile",
                                  `2` = "second pentile",
                                  `3` = "third pentile",
                                  `4` = "fourth pentile",
                                  `5` = "fifth pentile"
                                  )
                                )
      ) +
      ggtitle( nms[[i]] ) +
      theme_minimal(base_size = 22) +
      theme( legend.position = "none", plot.title = element_text(hjust = .5) )
      
  )


# arrange figure's subplots for saving
( f$exec_fun / f$epis_mem )

# save the figure
ggsave( here("figs","graphstract_epred.jpg"), dpi = 300, width = 1.5 * 13.1, height = 14.4 )

