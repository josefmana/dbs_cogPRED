# Extract & describe electrode locations of (a subset of) included patients

rm( list = ls() ) # clear environment

library(here) # directory management
library(tidyverse) # data wrangling
library(R.matlab) # read matlab files


# IN-HOUSE FUNCTIONS ----

# collapse electrode contact paramaters from original LEAD-DBS outcome to a nice table
parcol <-
  
  function( d1, d2, colnms, side = c("right","left"), contnm = "contact" ) {
    
    rbind.data.frame(
      d1 %>% `colnames<-`(colnms) %>% as.data.frame() %>% add_column( side = side[1] ) %>% rownames_to_column(contnm),
      d2 %>% `colnames<-`(colnms) %>% as.data.frame() %>% add_column( side = side[2] ) %>% rownames_to_column(contnm)
    )
    
  }


# PATIENTS PRESENCE ----

d0 <- readRDS( here("_data","longitudinal_df.rds") )$d0 # data set for analysis
elloc_meta <- read.csv( here("_data","elloc","elloc_subjects.csv"), sep = ";" ) # read metadata
elloc_pats <- list.dirs( here("_data","elloc","subs"), recursive = F ) %>% substr( nchar(.)-5, nchar(.) ) # patients from elloc files

# manually list reasons for missing localisation data for each patient
elloc_miss <-
  
  data.frame(
    id = sort( unique(d0$id)[ !(unique(d0$id) %in% elloc_pats ) ] ), # patients with missing elloc
    elloc_miss =
      c("no mri","no mri","post only","post only","post only","no mri","post only","no mri","post only","post only",
        "no mri","post only","post only","post only","post only","no mri","no mri","no mri","post only","post only",
        "post only","post only", rep("no mri",16),"post only", rep("no mri",5),"post only","post only", rep("no mri",3),
        "post only", rep("no mri",6)
        )
  )

# add IPN188 to missing patients due to missing sitmulation
elloc_miss <- elloc_miss %>% add_row( id = "IPN188", elloc_miss = "miss stim" )
elloc_pats <- elloc_pats[ elloc_pats != "IPN188" ]


# MATLAB PROCESSING ----

# read matlab files with electrode stats
d1 <-
  
  lapply(
    
    setNames(elloc_pats,elloc_pats),
    function(i)
      
      readMat( here("_data","elloc","subs",i,"ea_stats_v6.mat") )$ea.stats %>%
      `names<-`( c( "conmat", "conmat_inside_vox", "conmat_inside_hull", "patname", "atlases", "electrodes", "stimulation" ) )
    
  )


# DATA REFORMAT ----

for ( i in names(d1) ) {

  # extract IDs
  d1[[i]]$patname <- unique( unlist( d1[[i]][["patname"]] ) )
  d1[[i]]$electrodes[ , ,1]$name <- c( d1[[i]][["electrodes"]][ , ,1][["name"]] )

  # prepare atlas names
  d1[[i]]$atlases <- sub( ".nii.gz", "", unlist( d1[[i]][["atlases"]][ , ,1][["names"]] ) )
  atlnms <- d1[[i]][["atlases"]]
  
  
  ## conmats ----
  
  for ( j in names(d1[[i]])[ grepl( "conmat", names(d1[[i]]) ) ] ) {
    
    d1[[i]][[j]] <-
      
      parcol(
        d1 = d1[[i]][[j]][ ,1][[1]][[1]],
        d2 = d1[[i]][[j]][ ,2][[1]][[1]],
        colnms = atlnms
      )
    
  }
  
  
  ## electrodes ----
  
  # trajectories
  for ( j in c("coords.mm","trajectory") ) {
    
    d1[[i]]$electrodes[ , ,1][[j]] <-
      
      parcol(
        d1 = d1[[i]]$electrodes[ , ,1][[j]][ ,1][[1]][[1]],
        d2 = d1[[i]]$electrodes[ , ,1][[j]][ ,2][[1]][[1]],
        colnms = c("x","y","z"),
        contnm = ifelse( j == "coords.mm", "contact", "point" )
      )
    
  }
  
  # markers
  d1[[i]]$electrodes[ , ,1]$markers <-
    
    parcol(
      d1 = do.call( rbind.data.frame, d1[[i]][["electrodes"]][ , ,1][["markers"]][ , ,1] ),
      d2 = do.call( rbind.data.frame, d1[[i]][["electrodes"]][ , ,1][["markers"]][ , ,2] ),
      colnms = c("x","y","z")
    )
  
  # pull them together
  d1[[i]]$electrodes <- with( d1[[i]][["electrodes"]][ , ,1], list( coords.mm = coords.mm, trajectory = trajectory, name = name, markers = markers ) )
  
  
  ## stimulation ----
  
  d1[[i]]$stimulation <- d1[[i]][["stimulation"]][ , , length(d1[[i]][["stimulation"]])/3 ] # keep only the last try
  d1[[i]]$stimulation$label <- c( d1[[i]][["stimulation"]][["label"]] ) # extract session label
  
  #
  d1[[i]]$stimulation$contact <- unlist( d1[[i]][["stimulation"]][["vat"]][ , ,1]["contact", ] ) %>% `names<-`( c("right","left") )
  
  # voltage/current
  d1[[i]]$stimulation$amp <-
    
    cbind.data.frame(
      right = c( d1[[i]][["stimulation"]][["vat"]][ , ,1]["amp",1][["amp"]] ),
      left = c( d1[[i]][["stimulation"]][["vat"]][ , ,1]["amp",2][["amp"]] )
    )
  
  # volume
  d1[[i]]$stimulation$volume <-
    
    c(
      right = c( d1[[i]][["stimulation"]][["vat"]][ , ,1]["volume",1][["volume"]] ),
      left = c( d1[[i]][["stimulation"]][["vat"]][ , ,1]["volume",2][["volume"]] )
    )
  
  # Intersections
  d1[[i]]$stimulation$Interscetions <-
    
    lapply(
      
      rownames( d1[[i]]$stimulation$vat[ , ,1] )[ grepl( "Atlas", rownames( d1[[i]]$stimulation$vat[ , ,1] ) ) ],
      function(j)
        
        rbind.data.frame(
          right = d1[[i]][["stimulation"]][["vat"]][ , ,1][j, ][[1]],
          left = d1[[i]][["stimulation"]][["vat"]][ , ,1][j, ][[2]]
        ) %>%
        
        `colnames<-`(atlnms) %>%
        rownames_to_column("side") %>%
        add_column( var = j, .before = 1 ) %>%
        pivot_wider( names_from = side, values_from = all_of(atlnms) ) %>%
        column_to_rownames("var")
      
    ) %>%
    
    do.call( rbind.data.frame, . )
  
  # get rid of the rest
  d1[[i]][["stimulation"]][["vat"]] <- NULL
  d1[[i]][["stimulation"]][["efield"]] <- NULL
  
  # collapse the last level
  d1[[i]] <- d1[[i]][ , ,1]

}


# JOIN DATASETS ----

d2 <-
  
  # add info about electrode locations 
  d0 %>%
  left_join( elloc_miss, by = "id" ) %>%
  mutate( elloc = ifelse( is.na(elloc_miss), 1, 0 ) ) %>%
  
  # add time from surgery at time of VAT measurement
  left_join(
    read.csv( here("_data","20240330_redcap_data.csv"), sep = "," ) %>%
      rename( "id" = "study_id" ) %>%
      filter( redcap_event_name == "operace_arm_1" ) %>%
      filter( id %in% unique(d0$id) ) %>%
      select( id, surgery_date ) %>%
      left_join( elloc_meta %>% rename( "id" = "IPN" ) %>% select(id,date_mri) , by = "id" ) %>%
      mutate(
        surgery_date = as.Date( surgery_date, format = "%Y-%m-%d" ),
        date_mri = as.Date( date_mri, format = "%m/%d/%y" ),
        vat_time_y = time_length( difftime( date_mri, surgery_date ) , unit = "year" )
      ) %>%
      select( id, vat_time_y ),
    by = "id"
  ) %>%
  
  # add time difference from assessment to VAT measurement
  mutate( post = ifelse( ass_type == "pre", 0, 1 ), vat_after_y = time_y - vat_time_y ) %>%
  
  # add VAT/STN intersections
  left_join(

    sapply(

      names(d1),
      function(i)
        d1[[i]]$stimulation$Interscetions[1:2,1:8] %>% # only raw intersections and STN
        rownames_to_column("atlas") %>%
        pivot_wider( values_from = !atlas, names_from = atlas )

    ) %>%

      t() %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      mutate_if( is.list, as.numeric ) %>%
      mutate( across( everything(), ~ ifelse( is.nan(.x), NA, .x ) ) ), # recode NaNs
    
    by = "id"

  )

# SAVE IT ----

write.table( d2, file = here("_data","20240330_data_vat_plus.csv"), sep = ",", row.names = F, quote = F )
