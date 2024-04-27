# This is a script that prepares graphical abstract figure for the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

s = 87542 # seed for reproducibility

library(here) # directory management
library(patchwork) # putting graphs together
library(tidyverse) # data wrangling
library(tidybayes) # posterior manipulation


# POSTERIOR PREDICTION ----

## LOAD MATERIALS ----

for ( i in names( readRDS( here("_data","longitudinal_df.rds") ) ) ) assign( i, readRDS( here("_data","longitudinal_df.rds") )[[i]] ) # load all the data
m <- readRDS( here("mods","m1_lasso_doms.rds") ) # load the model
prds <- doms


## PRE-PROCESSING ----

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


## POSTERIOR PREDICTION ----

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
  write.table( d_seq, file = here("_data","graphstract_epred.csv"), sep = ",", row.names = F, quote = F )

} else {
  
  # if it was calculated in past, just load it
  d_seq <- read.csv( here("_data","graphstract_epred.csv"), sep = "," )
  
}


## PLOTTING ----

theme_set( theme_minimal( base_size = 14 ) ) # set-up theme
clr <- c( "#78A641", "#8A60B0", "#2CA030", "#FF7F0E", "#12A2A8", "#6F63BB", "#1F83B4" ) # set-up colors (grouped according to cognitive functions)

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

f3 <-
  
  # prepare data
  d %>%
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
  scale_y_continuous( name = "DRS-2", limits = c(80,155), breaks = seq(80,140,20), labels = seq(80,140,20) ) +
  scale_x_continuous( name = "Time post-surgery (years)", limits = c(-2,12), breaks = seq(-2,12,2), labels = seq(-2,12,2) ) +
  scale_colour_manual( values = clr ) +
  scale_fill_manual( values = clr ) +
  
  # final touches
  theme( legend.position = "none", panel.grid.minor = element_blank() ) +
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
ggsave( plot = f3, here("figs","pred_stratified.jpg"), dpi = 300, width = 11.1, height = 13.3 )


# DESCRIPTION ----

# read the data set and prepare subsets for individual analyses
d1 <- read.csv( here("_data","20220508_dbs_longCOG_data.csv") , sep = "," ) %>% filter( included == 1 )

# plot it
f1 <-
  
  d1[ complete.cases(d1$drs_tot) , ] %>% 
  ggplot( aes(x = time_y) ) +
  stat_bin( breaks = seq(-2,12,.5) ) + # creates bars
  stat_bin( breaks = seq(-2,12,.5), geom = "text", aes(label = after_stat(count) ), vjust = -0.8, size = 3.5 ) + # add numbers
  labs( x = "Time post-surgery (years)", y = "Number of Assessments" ) +
  scale_y_continuous( expand = c(0, 0), limits = c(0, 107), breaks = seq(0, 100, 10), labels = seq(0, 100, 10) ) +
  scale_x_continuous( limits = c(-2, 12), breaks = seq(-2, 12, 1), labels = seq(-2, 12, 1) ) +
  theme( panel.grid.minor = element_blank() )


# CONTRASTS ----

# prepare colors to use in graphs (a colorblind-friendly palette)
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )

## PREPARE DATA ----
d_seq <-
  
  data.frame( id = "sub000", time_y = c( -.3, seq(.5,5,1/12) ) ) %>%
  mutate( time = time_y + scl$Md$time )

# contrast for the true score
contr <-
  
  lapply(
    
    c("fixed","random"),
    function(i)
      
      # extract predictions
      d_seq %>%
      add_epred_draws(
        object = readRDS( here("mods","m0_linear.rds") ),
        newdata = . ,
        allow_new_levels = T,
        if ( i == "fixed" ) { re_formula = NA },
        seed = 87542
      ) %>%
      
      # post-process
      mutate( pred = scl$M$drs + scl$SD$drs * .epred ) %>%
      mutate( pred = case_when( pred > 144 ~ 144, pred < 0 ~ 0, .default = pred) ) %>%
      ungroup() %>%
      select( id, time_y, .draw, pred ) %>%
      pivot_wider( names_from = time_y, values_from = pred ) %>%
      mutate( across( all_of( as.character( d_seq$time_y ) ), ~ .x-`-0.3` ) ) %>%
      pivot_longer( cols = -c("id",".draw"), names_to = "contrast" ) %>%
      group_by(contrast) %>%
      summarise( md = median(value), eti_low = qi(value,.width=.9)[1], eti_hig = qi(value,.width=.9)[2] ) %>%
      mutate( model = i, .before = 1 )
    
  )

# extract predictive contrast of one year post-minus-pre assessments
contr$all <-
  
  d_seq %>%
  add_predicted_draws( object = readRDS( here("mods","m0_linear.rds") ), newdata = . , allow_new_levels = T, re_formula = NA, seed = 87542 ) %>%
  mutate( pred = scl$M$drs + scl$SD$drs * .prediction ) %>%
  mutate( pred = case_when( pred > 144 ~ 144, pred < 0 ~ 0, .default = pred ) ) %>%
  ungroup() %>%
  select( id, time_y, .draw, pred ) %>%
  pivot_wider( names_from = time_y, values_from = pred ) %>%
  mutate( across( all_of( as.character( d_seq$time_y ) ), ~ .x - `-0.3` ) ) %>%
  pivot_longer( cols = -c("id",".draw"), names_to = "contrast" ) %>%
  group_by(contrast) %>%
  summarise( md = median(value), eti_low = qi(value,.width=.9)[1], eti_hig = qi(value,.width=.9)[2] ) %>%
  mutate( model = "all", .before = 1 )

# put it together
contr <-
  
  contr %>%
  do.call( rbind.data.frame, . ) %>%
  filter( !grepl("-",contrast) ) %>%
  pivot_wider( names_from = model, values_from = c("md","eti_low","eti_hig") ) %>%
  mutate( contrast = as.numeric(contrast) )


## PLOT IT ----
f2 <-
  
  contr %>%
  ggplot() +
  aes( y = contrast, x = md_fixed ) +
  
  # add lines and shades
  geom_ribbon( data = contr, aes(xmin = eti_low_all, xmax = eti_hig_all), fill = cbPal[8], alpha = .4 ) +
  geom_ribbon( data = contr, aes(xmin = eti_low_random, xmax = eti_hig_random), fill = cbPal[8], alpha = .6 ) +
  geom_ribbon( data = contr, aes(xmin = eti_low_fixed, xmax = eti_hig_fixed), fill = cbPal[8], alpha = 1 ) +
  geom_line( size = 2, colour = "black" ) +
  
  # finish it
  geom_vline( xintercept = 0, linewidth = 1, linetype = "dashed", colour = "navyblue" ) +
  labs( y = "Time post-surgery (years)", x = parse( text = "Delta~`DRS-2`" ) ) +
  scale_x_continuous( breaks = seq(-15,5,1), labels = seq(-15,5,1) ) +
  coord_flip( xlim = c(-14,3) )

# save it
ggsave( plot = f2, file = here("figs","desc_contrasts.jpg"), dpi = 300, width = 9.66, height = 5.5 )


# MERGING ----

( (f1/f2) | f3 ) + plot_layout(widths = c(1,1.25) )
ggsave( plot = last_plot(), file = here("figs","graphstract.jpg"), dpi = 300, width = 18.3, height = 13.3 )
