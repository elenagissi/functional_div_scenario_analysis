#Gissi et al. - Functional diversity at grid cell level
# version September 13, 2023

## I use 127 species and 15 functional traits that 
# I have selected through exploring the correlation among them.


rm( list=ls())

#set working directory ----
setwd("~/MEDIX/WP2_Model/DATA")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("mFD")

library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(mFD)


# species x traits data frame ----
traits_fishb_test <-
  read_delim(
    "traits.csv",
    delim = ";",
    escape_double = FALSE,
    col_types = cols(
      length_max = col_number(),
      trophic_level = col_number(),
      depth_max = col_number()
    ),
    trim_ws = TRUE
  ) %>% 
  dplyr::select(
    species,
    length_max,
    trophic_level,
    depth_max, 
    sexual_dimorphism,
    sociality,
    reproductive_strategy,
    parental_investment,
    parental_care,
    terrestriality,
    adult_mobility,
    zones,
    body_shape,
    diet_main_component,
    ext_diet_spec,
    fecundity
  )

traits <- traits_fishb_test[1:127,] %>% 
  column_to_rownames( var ="species" )  # rows must take the name of species

# Nominal traits must be factors 
# Countinuous traits are numeric
str(traits)
traits$zones <- factor(traits$zones)
traits$adult_mobility <- factor(traits$adult_mobility) 
traits$reproductive_strategy <- factor(traits$reproductive_strategy)
traits$parental_investment <- factor(traits$parental_investment)
traits$diet_main_component <- factor(traits$diet_main_component)
traits$terrestriality <- factor(traits$terrestriality)
traits$parental_care <- factor(traits$parental_care, order = TRUE, levels = c("none", "week_month", "month_year", "more_t_year"))
traits$body_shape <- factor(traits$body_shape)
traits$sexual_dimorphism <- factor(traits$sexual_dimorphism)
traits$sociality <- factor(traits$sociality, order = TRUE, levels = c("solitary", "small_groups", "medium_groups", "large_groups"))
traits$ext_diet_spec <- factor(traits$ext_diet_spec)
traits$fecundity <-
  factor(traits$fecundity,
         order = TRUE,
         levels = c(
           "<1",
           "1",
           "1-2",
           "2-5",
           "5-10",
           "10-20",
           "20-50",
           "50-100",
           "100-1000",
           "1000-10000",
           ">10000"))

summary(traits)


# Aquamaps species identification ----
#occurrence file with the species distribution 
gissi_speciesoccursum <- read_delim("gissi_speciesoccursum_FINAL.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(speciesid = SpeciesID, aqmspeccode = SpecCode, genus = Genus, species = Species, 
         fbname= FBname, occurcells = OccurCells, kingdom = Kingdom, phylum = Phylum, 
         class = Class, order = Order, family = Family) %>% 
  mutate(spp = paste(genus, species, sep=" "))

gissi_speciesoccursum2 <- gissi_speciesoccursum$spp  


# Species distribution from aquamaps ----
scen_2022_current <- read_csv("species_occur.csv") %>% 
  filter (probability >= 0.5)
View(scen_2022_current)

scen_2022_current_abund <- scen_2022_current %>%   
  group_by (sci_name) %>% 
  mutate (abund = probability / sum(probability) * 10000)

test <- scen_2022_current_abund %>%  
  group_by(sci_name) %>% 
  summarize(tot_abund=sum(abund))


is.numeric(scen_2022_current$probability) 


# Assemblages by grid cells of 1 degree ----
# with the cells (cquarecode) at 1 degree
sppass_base_cells_deg1 <- scen_2022_current_abund %>% 
  dplyr::select(s_cells_1degree, sci_name, abund) %>% 
  pivot_wider(names_from = sci_name, values_from = abund,  values_fn=sum, values_fill = NA) %>% 
  column_to_rownames( var ="s_cells_1degree" )

sppass_base_cells_deg1[is.na(sppass_base_cells_deg1)] <- 0
is.matrix(sppass_base_cells_deg1)
sppass_base_cells_deg1 <- as.matrix((sppass_base_cells_deg1))
summary(sppass_base_cells_deg1)

# see if asb have less than 6 species:
rowSums(sppass_base_cells_deg1 != 0) 
less_six <- sppass_base_cells_deg1[which(rowSums(sppass_base_cells_deg1 != 0) < 6), ]
rowSums(less_six != 0)
nrow(less_six)

#SELECT ONLY CELLS WITH MORE THAN 5 SPECIES (meaning => 6 species each)
is.data.frame(less_six)
less_six_2=as.data.frame(less_six)
is.data.frame(less_six_2)
less_six_3 <- rownames_to_column(less_six_2, var = "s_cells_1degree")  %>% as_tibble()
is.data.frame(less_six_3)

sppass_base_cells_deg1_rev1=as.data.frame(sppass_base_cells_deg1)
is.data.frame(sppass_base_cells_deg1_rev1) 
sppass_base_cells_deg1_rev2 <- rownames_to_column(sppass_base_cells_deg1_rev1, var = "s_cells_1degree")

sppass_base_cells_deg1_rev3 <-
  anti_join (sppass_base_cells_deg1_rev2, less_six_3, by = "s_cells_1degree") %>%
  column_to_rownames(var = "s_cells_1degree")

#Now I re-trasform the file in a matric to keep working on the analysis of FD
is.matrix(sppass_base_cells_deg1_rev3)
sppass_base_cells_deg1_rev3 <- as.matrix((sppass_base_cells_deg1_rev3))
summary(sppass_base_cells_deg1_rev3)




## Trait type ----
#I upload the data.frame containing the description of the traits for this scenario (INPUT 3)
traits_measure <- read_delim("traits_measure.csv", 
                             delim = ";", escape_double = FALSE, col_types = cols(`Classes code` = col_integer()), 
                             trim_ws = TRUE) 
View(traits_measure)

traits_type <- traits_measure %>%   
  dplyr::select(trait_name, trait_type)

# FUNCTIONAL DIVERSITY ANALYSIS START ----

## SUMMARIZE MY TRAITS ----

#CALCULATION:  Species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_type,   
  sp_tr      = traits, 
  stop_if_NA = FALSE)

#Traits type resulting from the species traits summary:
traits_summ$"tr_types" 
traits_summ$"mod_list"  

## SUMMARIZE ASSEMBLAGES  ----


asb_sp_species_deg1_summ <- mFD::asb.sp.summary(asb_sp_w = sppass_base_cells_deg1_rev3)

### IMPORTANT: Species occurrence per each cell:
# this is a matrix with the occurrence per each cell (0,1) of the species 
asb_sp_species_occ_deg1 <- asb_sp_species_deg1_summ$"asb_sp_occ" 

# IMPORTANT: Species total biomass in all assemblages 
#(sum of all occurrences for the entire case study area):
asb_sp_species_deg1_summ$"sp_tot_w"  

# IMPORTANT:Total biomass per assemblage (sum of occurrences of all the species by cell):
# all species occurrence is summed by each species
asb_sp_species_deg1_summ$"asb_tot_w" 

# IMPORTANT: Species richness per assemblage (no. of species per each assemblage/cell)
asb_sp_species_deg1_summ$"asb_sp_richn" 


# FUNCTIONAL DISTANCE ----
# COMPUTING DISTANCES BETWEEN SPECIES BASED ON FUNCTIONAL TRAITS

sp_dist_species <- mFD::funct.dist(
  sp_tr         = traits,
  tr_cat        = traits_type,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE) %>% 
  round(3)

sp_dist_species
image(t(as.matrix(sp_dist_species)))
View(sp_dist_species)

# FUNCTIONAL SPACES & THEIR QUALITY ----
fspaces_quality_species <- mFD::quality.fspaces(
  sp_dist             = sp_dist_species,
  maxdim_pcoa         = 10,
  deviation_weighting = c('absolute', 'squarred'),
  fdist_scaling       = FALSE,
  fdendro             = "average") 

round(fspaces_quality_species$"quality_fspaces", 3)   

# display the table gathering quality metrics:
fspaces_quality_species$"quality_fspaces"

## retrieve the functional space associated with minimal quality metric: 
apply(fspaces_quality_species$quality_fspaces, 2, which.min)

percent_tot_variance_explained <-
  (
    fspaces_quality_species[["details_fspaces"]][["pc_eigenvalues"]]$Eigenvalues /
      sum(fspaces_quality_species[["details_fspaces"]][["pc_eigenvalues"]]$Eigenvalues)
  ) * 100
percent_tot_variance_explained

cum_percent_tot_variance_explained <-  c(
  percent_tot_variance_explained[1],
  percent_tot_variance_explained[1] + percent_tot_variance_explained[2],
  percent_tot_variance_explained[1] + percent_tot_variance_explained[2] +
    percent_tot_variance_explained[3],
  percent_tot_variance_explained[1] + percent_tot_variance_explained[2] +
    percent_tot_variance_explained[3] + percent_tot_variance_explained[4],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5] +
    percent_tot_variance_explained[6],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5] +
    percent_tot_variance_explained[6] + percent_tot_variance_explained[7],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5] +
    percent_tot_variance_explained[6] + percent_tot_variance_explained[7] +
    percent_tot_variance_explained[8],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5] +
    percent_tot_variance_explained[6] + percent_tot_variance_explained[7] +
    percent_tot_variance_explained[8] + percent_tot_variance_explained[9],
  percent_tot_variance_explained[1] +
    percent_tot_variance_explained[2] + percent_tot_variance_explained[3] +
    percent_tot_variance_explained[4] + percent_tot_variance_explained[5] +
    percent_tot_variance_explained[6] + percent_tot_variance_explained[7] +
    percent_tot_variance_explained[8] + percent_tot_variance_explained[9] +
    percent_tot_variance_explained[10]
)
cum_percent_tot_variance_explained 


#I can plot the quality matrix of each space:
library("magrittr")

fspaces_quality_species$"quality_fspaces" %>%
  tibble::as_tibble(rownames = "Funct.space") %>%
  tidyr::pivot_longer(cols =! Funct.space, names_to = "quality_metric", values_to = "Quality") %>%
  ggplot2::ggplot(ggplot2::aes(x = Funct.space, y = Quality, 
                               color = quality_metric, shape = quality_metric)) +
  ggplot2::geom_point() 


# Quality metrics of spaces

# Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_species,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_3d","pcoa_4d", "pcoa_5d",
    "pcoa_6d", "pcoa_7d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance"
)


## TEST CORRELATION BETWEEN FUNCTIONAL AXES AND TRAITS ---- 

sp_faxes_coord_species <- fspaces_quality_species$"details_fspaces"$"sp_pc_coord"

species_tr_faxesA <- mFD::traits.faxes.cor(
  sp_tr          = traits[,1:10], 
  sp_faxes_coord = sp_faxes_coord_species[ , c( "PC1", "PC2", "PC3", "PC4", "PC5")], 
  plot           = TRUE)

# Print traits with significant effect:
species_tr_faxesA$"tr_faxes_stat"[which(species_tr_faxesA$"tr_faxes_stat"$"p.value" < 0.05), ]

# Return plots:
species_tr_faxesA$"tr_faxes_plot"

#I plot the other traits (in groups of 10)
species_tr_faxesB <- mFD::traits.faxes.cor(
  sp_tr          = traits[,11:15], #I plot for the first 10 traits
  sp_faxes_coord = sp_faxes_coord_species[ , c( "PC1", "PC2", "PC3", "PC4", "PC5")], 
  plot           = TRUE)

# Print traits with significant effect:
species_tr_faxesB$"tr_faxes_stat"[which(species_tr_faxesB$"tr_faxes_stat"$"p.value" < 0.05), ]

# Return plots:
species_tr_faxesB$"tr_faxes_plot"

## PLOT FUNCTIONAL SPACE   ---- 
big_plot4 <- mFD::funct.space.plot( #plot for 4 dimensions
  sp_faxes_coord  = sp_faxes_coord_species[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot4$patchwork

#Plot in 10 dimensions
big_plot_all <- mFD::funct.space.plot( #plot for 10 dimensions
  sp_faxes_coord  = sp_faxes_coord_species,
  faxes           = NULL,
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot_all$patchwork   



##  FUNCTIONAL DIVERTITY INDICES  ----
### Functional alpha diversity indices in a multidimensional space ----

#The mFD::alpha.fd.multidim() function allows computing alpha FD indices:
alpha_fd_indices_species_deg1 <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_species[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w         = sppass_base_cells_deg1_rev3,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE, 
  details_returned = TRUE)

#END

