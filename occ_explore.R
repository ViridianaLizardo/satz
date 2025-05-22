## Script name: Occurrence exploration of plant occurrences
## using rGBIF
##
## Author: Viridiana Lizardo
##
## Copyright (c) Viridiana Lizardo
## Email: lizardo.viridiana@gmail.com

library(taxize)
library(rgbif)
library(tidyverse)

# x<- 'Opuntia_robusta'

# 01 Function to get species ----
occ_explore<- function(x){
  # The 'x' means a species name in the format "Genus_species"
  res <- taxize::get_gbifid_(str_replace(x, '_', ' '),
                             method = 'backbone',messages = T)
  names(res) <- NULL
  
  # Filter out species that are not plants
  res_df <- as.data.frame(res) %>%
    filter(kingdom == 'Plantae', 
           rank == 'species', confidence >90)
  
  # Make sure the result is an accepted and exact species name.
  
  if('species' %in% names(res_df)){
    ids_df <-  res_df %>%
      select(usagekey, species,
             family, order, matchtype,
             status )%>%
      mutate(label = x) %>%
      filter(status == 'ACCEPTED' & 
               matchtype == 'EXACT')
  }else{cat('Not found\n')}# If it is not, it prints out a message
  
  
  # Check that it has a single result
  if(nrow(ids_df) == 1){
    res_occ <- occ_search(taxonKey = ids_df$usagekey,
                          limit = 300,
                          # polygon of interest, can be removed
                          geometry = wkt_polygon, 
                          hasCoordinate = T,
                          # just lat and lon data 
                          fields = 'minimal')$data
    # Count number of occurrences
    # Note that the maximum is 300
    if(is.null(res_occ) == T){
      n_occ <- 0} else{n_occ <- nrow(res_occ)}
    
    ids_df %>%
      mutate(occ = n_occ,
             label = x) %>% 
      # add the result to a csv with the keys and the number of 
      # found records
      write.table('gbif_keys.csv',sep = ",",
                  append = TRUE, quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    
  }else{
    # In case of multiple matches or not founding anything,
    # it saves the result to a file named 'not found.csv'
    # where all the name matches are included.
    res_df %>%
      mutate(label = x) %>% 
      select('usagekey',	'species',
             'status',	'matchtype',	'label') %>%
      write.table(., 'not found.csv',sep = ",",
                  append = TRUE, quote = FALSE,
                  col.names = F, row.names = FALSE)
    
    cat('Not found\n')
    
  }
}

## Apply the function through a list with error handling.----

gbif_keys <- function(x){try(occ_explore(x),silent = T)}


# 02 Synonims ----
## Accepted usagekey
occ_explore_02 <- function(x){
  res <- taxize::get_gbifid_(str_replace(x, '_', ' '),
                             method = 'backbone',messages = T)
  names(res) <- NULL
  
  # Filter out species that are not plants
  res_df <- as.data.frame(res) %>%
    filter(kingdom == 'Plantae', 
           rank == 'species', confidence >90) 
  
  # When specieskey is found
  if('acceptedusagekey' %in% names (res_df)){
    keys <- pull(res_df, acceptedusagekey ) %>% 
      unique()
    keys <- keys[!is.na(keys)]
    if(length(keys) ==1){ 
      res_occ <- occ_search(taxonKey = keys,
                            limit = 300,
                            # polygon of interest, can be removed
                            geometry = wkt_polygon,
                            hasCoordinate = T,
                            # just lat and lon data
                            fields = 'minimal')$data
      # Count number of occurrences
      # Note that the maximum is 300
      if(is.null(res_occ) == T){
        n_occ <- 0} else{n_occ <- nrow(res_occ)}}}else{
          write.table(x, 'No acceptedusagekey.csv',sep = ",",
                      append = TRUE, quote = FALSE,
                      col.names = F, row.names = FALSE)
          message('No acceptedusagekey')}
  
  if(is.null(res_occ)){
    write.table(x, 'forget about it.csv',sep = ",",
                append = TRUE, quote = FALSE,
                col.names = F, row.names = FALSE)
    
    cat('forget about it\n')}else{
      ids_df  <- res_df%>%
        select(acceptedusagekey, species,
               family, order, matchtype,
               status )%>% distinct() %>%
        mutate(label = x)%>%
        mutate(occ = n_occ,
               label = x) %>% 
        na.omit() %>%
        filter(matchtype == 'EXACT')%>%
        # add the result to a csv with the keys and the number of 
        # found records
        write.table('gbif_keys_06.csv',sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
    }}

gbif_keys_synonims <- function(x){try(occ_explore_02(x))}

#lapply(multiple_matches, gbif_keys_synonims)

## no DOUBTFUL
occ_explore_03 <- function(x){
  res <- taxize::get_gbifid_(str_replace(x, '_', ' '),
                             method = 'backbone',messages = T)
  names(res) <- NULL
  
  # Filter out species that are not plants no DOUBTFUL
  res_df <- as.data.frame(res) %>%
    filter(kingdom == 'Plantae', 
           rank == 'species', confidence >90,
           status != 'DOUBTFUL') 
  
  if(nrow(res_df) ==1){ 
    res_occ <- occ_search(taxonKey = res_df$usagekey,
                          limit = 300,
                          # polygon of interest, can be removed
                          geometry = wkt_polygon,
                          hasCoordinate = T,
                          # just lat and lon data
                          fields = 'minimal')$data
    # Count number of occurrences
    # Note that the maximum is 300
    if(is.null(res_occ) == T){
      n_occ <- 0} else{n_occ <- nrow(res_occ)}}else{
        write.table(x, 'No acceptedusagekey_02.csv',sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = F, row.names = FALSE)
        message('No acceptedusagekey')}
  
  if(is.null(res_occ)){
    write.table(x, 'forget about it.csv',sep = ",",
                append = TRUE, quote = FALSE,
                col.names = F, row.names = FALSE)
    
    cat('forget about it\n')}else{
      ids_df  <- res_df%>%
        select(usagekey, species,
               family, order, matchtype,
               status )%>% distinct() %>%
        mutate(label = x)%>%
        mutate(occ_mx = n_occ,
               label = x) %>% 
        # add the result to a csv with the keys and the number of 
        # found records
        write.table('gbif_keys_03.csv',sep = ",",
                    append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
    }}

gbif_keys_synonims02 <- function(x){try(occ_explore_03(x))}