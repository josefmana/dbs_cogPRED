# This is a script that computes generalized mixed models to describe and predict post-surgery cognitive decline for
# the longitudinal cognition in DBS study

rm( list = ls() ) # clear environment

# list packages to be used
pkgs <- c("here","tidyverse","brms","tidybayes")

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders to store results in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present.
sapply( c("mods","figs","tabs"), function(i) if( !dir.exists(i) ) dir.create(i) )
