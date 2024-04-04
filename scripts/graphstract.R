# This is a script that prepares graphical abstract figure for the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

s = 87542 # seed for reproducibility

library(here) # directory management
library(tidyverse) # data wrangling
library(tidybayes) # posterior manipulation


# LOAD MATERIALS ----

for ( i in names( readRDS( here("_data","longitudinal_df.rds") ) ) ) assign( i, readRDS( here("_data","longitudinal_df.rds") )[[i]] ) # load all the data
m <- readRDS( here("mods","m1_lasso_doms.rds") ) # load the model
prds <- doms
#prds <- c("exec_fun","epis_mem") # list of predictors effect which we're going to visualize


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
if ( !file.exists( here("_data","graphstract_epred.csv") ) ) {
  
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
        mutate(.epred = scl$M$drs + scl$SD$drs * .epred ) %>%
        median_hdi(.width = .95) %>%
        mutate( predictor = j, type = i )
      
    } 
  }
  
  # pull it all together
  d_seq <-
    
    do.call(
      rbind.data.frame,
      list(
        do.call( rbind.data.frame, d_seq$epred_fix ),
        do.call( rbind.data.frame, d_seq$epred_all )
      )
    )
  
  # save it
  #saveRDS( d_seq, here("_data","graphstract_epred.rds") )
  write.table( d_seq, file = here("_data","graphstract_epred.csv"), sep = ",", row.names = F, quote = F )

} else {
  
  # if it was calculated in past, just load it
  #d_seq <- readRDS( here("_data","graphstract_epred.rds") )
  d_seq <- read.csv( here("_data","graphstract_epred.csv"), sep = "," )
  
}


# PLOTTING ----

# set-up colors (grouped according to cognitive functions)
clr <- c( "#78A641", "#8A60B0", "#2CA030", "#12A2A8", "#FF7F0E", "#6F63BB", "#1F83B4" )

# extract order of variables from the most to the least predictive
ord <-
  
  m %>%
  spread_draws( `b_.*`, regex = T ) %>%
  select( contains("time:") ) %>%
  rename_with( ~ sub("b_time:","",.x) ) %>%
  apply( 2, median ) %>%
  sort() %>%
  names()

# make predictor an ordered factor in posterior predictive dataset
d_seq <- d_seq %>% mutate( predictor = factor( predictor, levels = ord, ordered = T ) )

# prepare figure for the graphical abstract
d %>%
  
  # prepare data
  pivot_longer(
    cols = ends_with("pent"),
    names_to = "predictor",
    names_transform = function(x) sub("_pent","",x), values_to = "grp"
  ) %>%
  mutate( predictor = factor( predictor, levels = ord, ordered = T ) ) %>%
  
  # plot raw data
  ggplot( aes(x = time_y, y = drs, group = id) ) +
  geom_hline( yintercept = 139, linetype = "dotted", linewidth = 1, alpha = 1 ) +
  geom_line( linewidth = .5, alpha = .33 ) +
  geom_point( size = 2.5, alpha = .33 ) +
  
  # plot posterior predictions
  geom_ribbon( data = subset( d_seq, type == "epred_all" ), alpha = .2, aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = predictor ) ) +
  geom_ribbon( data = subset( d_seq, type == "epred_fix" ), alpha = .3, aes(x = time_y, y = .epred, ymin = .lower, ymax = .upper, fill = predictor ) ) +
  geom_line( data = subset( d_seq, type == "epred_fix" ), linewidth = 1, aes( x = time_y, y = .epred, colour = predictor ) ) +
  
  # tidy-up axes and colours
  scale_y_continuous( name = "DRS-2", limits = c(80,153), breaks = seq(80,140,20), labels = seq(80,140,20) ) +
  scale_x_continuous( name = "Time from surgery (years)", limits = c(-2,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
  scale_colour_manual( values = clr ) +
  scale_fill_manual( values = clr ) +
  
  # final touches
  theme_minimal( base_size = 16 ) +
  theme( legend.position = "none" ) +
  facet_grid(
    rows = vars(predictor),
    cols = vars(grp),
    labeller =
      as_labeller(
        c( `1` = "1^th~pentile",
           `2` = "2^nd~pentile",
           `3` = "3^rd~pentile",
           `4` = "4^th~pentile",
           `5` = "5^th~pentile",
           exec_fun = "Executive~functions",
           visp_mem = "Visuospatial~memory",
           set_shift = "Set~shifting",
           epis_mem = "Episodic~memory",
           anxiety = "Anxiety",
           visp_wm = "Spatial~working~memory",
           verb_wm = "Verbal~working~memory"
        ),
        label_parsed
      )
  )

# save the figure
ggsave( here("figs","graphstract_epred.jpg"), dpi = 300, width = 12.3, height = 14.7 )
