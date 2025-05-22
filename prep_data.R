# SATZ - Occurrence & Phylogeny Pipeline
# --------------------------------------
# Este script descarga, limpia y organiza datos de ocurrencia y filogenia
# para una familia taxonómica definida, dentro de la South American Transition Zone (SATZ).

# -------------------------------
# Configuración inicial
# -------------------------------

# Parámetro configurable
fam <- "Poaceae"  # <- Cambia esto si quieres trabajar otra familia

# Librerías necesarias
library(tidyverse)
library(ape)
library(ggtree)
library(deeptime)
library(rgbif)
library(CoordinateCleaner)
library(sf)
library(tmap)
library(kewr)
library(pbapply)
library(V.PhyloMaker2)

sf_use_s2(FALSE)
data("World")
source('C:/Users/USER/Desktop/R projects/theme_virix.r')

# -------------------------------
# Filogenia
# -------------------------------

tree_all <- GBOTB.extended.WP
tip_species <- tips.info.TPL %>% filter(family == fam)
to_prune <- tree_all$tip.label[!tree_all$tip.label %in% tip_species$species]
tree_fam <- drop.tip(tree_all, to_prune)

# -------------------------------
# Región SATZ
# -------------------------------

satz_sf <- st_read('data/SATZ.gpkg') 
wkt_polygon <- st_as_text(satz_sf$geom) %>% wkt_parse()

# -------------------------------
# Validación taxonómica en GBIF
# -------------------------------

gbif_keys_all <- tree_fam$tip.label %>% 
  str_replace("_", " ") %>%
  name_backbone_checklist(strict = TRUE)

gbif_keys_all_clean <- gbif_keys_all %>%
  mutate(label = tree_fam$tip.label) %>%
  filter(!is.na(speciesKey))

# -------------------------------
# Descarga de ocurrencias (Exploración)
# -------------------------------

occ_explore <- function(x){
  res_occ <- occ_search(
    taxonKey = x,
    limit = 100,
    geometry = wkt_polygon, 
    hasCoordinate = TRUE,
    fields = 'minimal')$data
  
  if(is.null(res_occ)){
    output <- data.frame(speciesKey = x, occ = 0)
  } else {
    output <- res_occ %>% reframe(occ = n()) %>% mutate(speciesKey = x, .before = occ)
  }
  return(output)
}

res_gbif <- pblapply(gbif_keys_all_clean$speciesKey, occ_explore) %>%
  bind_rows() %>% left_join(gbif_keys_all) %>%
  filter(occ > 0)

# -------------------------------
# Descarga masiva desde GBIF
# -------------------------------

keys <- res_gbif %>% pull(speciesKey)

gbif_download <- occ_download(
  type="and",
  pred_in("taxonKey", keys),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_or(pred_not(pred("establishmentMeans", "MANAGED")), pred_not(pred_notnull("establishmentMeans"))),
  pred_or(pred_not(pred("establishmentMeans", "INTRODUCED")), pred_not(pred_notnull("establishmentMeans"))),
  pred_or(pred_not(pred("establishmentMeans", "INVASIVE")), pred_not(pred_notnull("establishmentMeans"))),
  pred_or(pred_not(pred("establishmentMeans", "NATURALISED")), pred_not(pred_notnull("establishmentMeans"))),
  pred_or(pred_lt("coordinateUncertaintyInMeters", 1000), pred_not(pred_notnull("coordinateUncertaintyInMeters"))),
  format = "SIMPLE_CSV"
)

occ_download_wait(gbif_download)

gbif_result <- occ_download_get(key = gbif_download) %>% occ_download_import()

# -------------------------------
# Limpieza de datos
# -------------------------------

gbif_result_clean <- gbif_result %>%
  select(occurrenceID, speciesKey, order, family, species,
         decimalLatitude, decimalLongitude, establishmentMeans) %>%
  cc_val(lat = 'decimalLatitude', lon = 'decimalLongitude') %>%
  cc_sea(lat = 'decimalLatitude', lon = 'decimalLongitude') %>%
  cc_equ(lat = 'decimalLatitude', lon = 'decimalLongitude') %>%
  cc_inst(lat = 'decimalLatitude', lon = 'decimalLongitude') %>%
  cc_zero(lat = 'decimalLatitude', lon = 'decimalLongitude') %>%
  cc_dupl(lat = 'decimalLatitude', lon = 'decimalLongitude')

# -------------------------------
# Conversión a objeto espacial
# -------------------------------

gbif_sf <- gbif_result_clean %>%
  mutate(lon = decimalLongitude,
         lat = decimalLatitude) %>%
  filter(!establishmentMeans %in% c('Introduced','Uncertain', 'Vagrant', 'NativeReintroduced')) %>% 
  st_as_sf(coords = c('decimalLongitude','decimalLatitude'), crs = st_crs(satz_sf)) %>%
  st_join(satz_sf) %>% 
  mutate(in_satz = if_else(is.na(Region), 0, 1)) %>%
  rownames_to_column('ID') %>%
  group_by(species) %>%
  mutate(n = n_distinct(ID),
         satz = sum(in_satz, na.rm = TRUE),
         perc_zt = satz / n,
         endem_class = cut(perc_zt, breaks = c(0, 0.05, 0.1, 0.5, 0.99, 1),
                           include.lowest = TRUE, right = FALSE))

# -------------------------------
# Exportar resultados
# -------------------------------

taxonomy <- st_drop_geometry(gbif_sf) %>%
  group_by(speciesKey, species, endem_class, perc_zt) %>%
  tally() %>%
  inner_join(gbif_keys_all_clean) %>%
  select(-verbatim_name, -verbatim_index, -acceptedUsageKey)

write.csv(taxonomy, paste0('data/', tolower(fam), '_taxonomy.csv'))

final_tree <- drop.tip(tree_fam, tree_fam$tip.label[!tree_fam$tip.label %in% taxonomy$label])
write.nexus(final_tree, file = paste0('data/', tolower(fam), '_final_tree.nex'))

st_drop_geometry(gbif_sf) %>%
  filter(perc_zt > 0.05) %>%
  write.csv(paste0('data/', tolower(fam), '_occurrences.csv'), row.names = FALSE)
