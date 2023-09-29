# Gissi et al. 
# Functional diversity indices by ecoregions
# Version Sept. 13, 2023

# Here I calculate fucntional diversity indices by ecoregions for all the scenarios
# ...with 127 species and 15 functional traits that 
# ... were have selected through exploring the correlation among them.
# ... We used the abundance of species.

#Set working directory ----

rm( list=ls())

setwd("~/MEDIX")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("mFD")

library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(mFD)


# Species x traits dataframe ----

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
  column_to_rownames( var ="species" )  


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



# Aquamaps species identification ----
gissi_speciesoccursum <- read_delim("gissi_speciesoccursum_FINAL.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(speciesid = SpeciesID, aqmspeccode = SpecCode, genus = Genus, species = Species, 
         fbname= FBname, occurcells = OccurCells, kingdom = Kingdom, phylum = Phylum, 
         class = Class, order = Order, family = Family) %>% 
  mutate(spp = paste(genus, species, sep=" "))

View(gissi_speciesoccursum) 

gissi_speciesoccursum2 <- gissi_speciesoccursum$spp  #I will do this only to check that everything is fine
# in this version gissi_speciesoccursum2=spp


# Species distribution from Aquamaps ----
scen_2022_current <- read_csv("species_occur.csv") %>% 
  filter (probability >= 0.5)
View(scen_2022_current)

scen_2022_current_abund <- scen_2022_current %>%   
  mutate (abund = probability / sum(probability) * 10000)

test <- scen_2022_current_abund %>%  
  group_by(sci_name) %>% 
  summarize(tot_abund=sum(abund))

is.numeric(scen_2022_current$probability) 



# GROUP BY ECOREGION  
s_base_ecoreg <- scen_2022_current_abund %>% 
  dplyr::select(s_ecoreg, sci_name, abund) %>%  
  pivot_wider(names_from = sci_name, values_from = abund, values_fn=sum, values_fill = NA) %>% 
  filter(s_ecoreg != "Clipperton _s210085") %>% 
  column_to_rownames( var ="s_ecoreg" ) 

s_base_ecoreg[is.na(s_base_ecoreg)] <- 0
is.matrix(s_base_ecoreg)
s_base_ecoreg <- as.matrix((s_base_ecoreg))
summary(s_base_ecoreg)

# Check if asb have less than 6 species:
rowSums(s_base_ecoreg != 0) 
less_six <- s_base_ecoreg[which(rowSums(s_base_ecoreg != 0) < 6), ]
rowSums(less_six != 0)
nrow(less_six)



# Trait type ----
#I upload the dataframe containing the description of the traits for this scenario (INPUT 3)
traits_measure <- read_delim("traits_measure.csv", 
                             delim = ";", escape_double = FALSE, col_types = cols(`Classes code` = col_integer()), 
                             trim_ws = TRUE) 
View(traits_measure)

traits_type <- traits_measure %>%   
  dplyr::select(trait_name, trait_type)



# FUNCTIONAL DIVERSITY ANALYSIS ----

## Summarize my traits ----

#Species traits summary
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_type,   
  sp_tr      = traits, 
  stop_if_NA = FALSE)

#Traits type resulting from the species traits summary
traits_summ$"tr_types" 
traits_summ$"mod_list"  


## Summary of the assemblages * species dataframe ----
asb_sp_species_ecoreg_summ <- mFD::asb.sp.summary(asb_sp_w = s_base_ecoreg)

### IMPORTANT: Species occurrence per each cell:
asb_sp_species_ecoreg_occ <- asb_sp_species_ecoreg_summ$"asb_sp_occ" 

# IMPORTANT: Species total biomass in all assemblages 
asb_sp_species_ecoreg_summ$"sp_tot_w"  

# IMPORTANT:Total biomass per assemblage (sum of occurrences of all the species by cell):
asb_sp_species_ecoreg_summ$"asb_tot_w" 

# IMPORTANT: Species richness per assemblage (no. of species per each assemblage/cell)
asb_sp_species_ecoreg_summ$"asb_sp_richn" 



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
#This function returns a dist object with traits-based distances between all pairs of species:
image(t(as.matrix(sp_dist_species)))
View(sp_dist_species)

# FUNCTIONAL SPACES & THEIR QUALITY ----
fspaces_quality_species <- mFD::quality.fspaces(
  sp_dist             = sp_dist_species,
  maxdim_pcoa         = 10,
  deviation_weighting = c('absolute', 'squared'),
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

# cumulative percentage of total variance explained to decide what PCs consider
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
  ggplot2::geom_point() +
  ggplot2::theme (panel.background = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
         axis.line = element_line(color = "black"))

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

#TEST CORRELATION BETWEEN FUNCTIONAL AXES AND TRAITS ---- 

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



# Plot species and axis in functional space ----
gissi_speciesoccursum3 <- gissi_speciesoccursum %>% 
  dplyr::select(spp, class, order, family, fbname)

sp_faxes_coord_species2 <- sp_faxes_coord_species %>% 
  as.data.frame() %>% 
  rownames_to_column("spp") %>% 
  left_join(gissi_speciesoccursum3, by="spp")

install.packages("plotly")
install.packages("gapminder")
library(ggplot2)
library(plotly)
library(gapminder)

###PC1 - PC2 ----
sp_faxes_coord_species2_plotPC1PC2 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC1, y=PC2, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, size = 2,
    check_overlap = T
  ) +
  xlim(-0.3,0.35)+ylim(-0.3,0.25)+
    theme_bw() +             
  ylab("PC2") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC1PC2)


#the same PC1 - PC2 plot but without text labels:
sp_faxes_coord_species2_plotPC1PC2_notext <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC1, y=PC2, group=fbname)) + #I put here fbname and it works
  geom_point(aes(shape=class, color=order), size = 3) +
  xlim(-0.35,0.35)+ylim(-0.4,0.25)+
  theme_bw() +             
  ylab("PC2") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC1PC2_notext)

###PC1 - PC3 ----
sp_faxes_coord_species2_plotPC1PC3 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC1, y=PC3, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, 
    check_overlap = T
  ) +
  xlim(-0.3,0.35)+ylim(-0.3,0.3)+
  theme_bw() +             
  ylab("PC3") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

sp_faxes_coord_species2_plotPC1PC3

# with names
sp_faxes_coord_species2_plotPC1PC3_text <- sp_faxes_coord_species2 %>%
  select(PC1, PC3, order, class, fbname) %>% 
  ggplot(aes(x=PC1, y=PC3, group=fbname)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  xlim(-0.3,0.35)+ylim(-0.3,0.3)+
  theme_bw() +             
  ylab("PC3") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC1PC3_text)

### PC2 - PC3 ----
sp_faxes_coord_species2_plotPC2PC3 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC2, y=PC3, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, 
    check_overlap = T
  ) +
  xlim(-0.35,0.35)+ylim(-0.3,0.25)+
  theme_bw() +             
  ylab("PC3") +   
  xlab("PC2") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

sp_faxes_coord_species2_plotPC2PC3

#with ggplotly to see the species names:
sp_faxes_coord_species2_plotPC2PC3_text <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC2, y=PC3, group=fbname)) +
  geom_point(aes(shape=class, color=order), size = 3) +
    xlim(-0.35,0.35)+ylim(-0.3,0.25)+
  theme_bw() +             
  ylab("PC3") +   
  xlab("PC2") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC2PC3_text)



### PC1 - PC4 ----
sp_faxes_coord_species2_plotPC1PC4 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC1, y=PC4, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, 
    check_overlap = T
  ) +
  xlim(-0.3,0.35)+ylim(-0.36,0.3)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

sp_faxes_coord_species2_plotPC1PC4


#PC1 and PC4 with plotly
sp_faxes_coord_species2_plotPC1PC4_text <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC1, y=PC4, group=fbname)) +
  geom_point(aes(shape=class, color=order), size = 3) +
    xlim(-0.3,0.35)+ylim(-0.36,0.3)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC1") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC1PC4_text)




###PC2 - PC4 ----
sp_faxes_coord_species2_plotPC2PC4 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC2, y=PC4, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, 
    check_overlap = T
  ) +
  xlim(-0.35,0.25)+ylim(-0.38,0.2)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC2") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

sp_faxes_coord_species2_plotPC2PC4


#PC2 - PC4 with plotly
sp_faxes_coord_species2_plotPC2PC4_text <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC2, y=PC4, group=fbname)) +
  geom_point(aes(shape=class, color=order), size = 3) +
    xlim(-0.35,0.25)+ylim(-0.38,0.2)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC2") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC2PC4_text)

### PC3 - PC4 ----
sp_faxes_coord_species2_plotPC3PC4 <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC3, y=PC4, group=order)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  geom_text(
    label=sp_faxes_coord_species2$fbname, 
    nudge_x = 0.01, nudge_y = 0.01, 
    check_overlap = T
  ) +
  xlim(-0.25,0.28)+ylim(-0.38,0.2)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC3") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

sp_faxes_coord_species2_plotPC3PC4

#PC3 - PC4
sp_faxes_coord_species2_plotPC3PC4_text <- sp_faxes_coord_species2 %>%
  ggplot(aes(x=PC3, y=PC4, group=fbname)) +
  geom_point(aes(shape=class, color=order), size = 3) +
  xlim(-0.25,0.28)+ylim(-0.38,0.2)+
  theme_bw() +             
  ylab("PC4") +   
  xlab("PC3") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggplotly(sp_faxes_coord_species2_plotPC3PC4_text)

# PLOT FUNCTIONAL SPACE   ---- 
big_plot4 <- mFD::funct.space.plot( #plot for 4 dimensions
  sp_faxes_coord  = sp_faxes_coord_species[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "#5ec962",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "#3b528b",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "#3b528b",
  fill_vert       = "#3b528b",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "#3b528b",
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
  color_pool      = "#5ec962",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "#3b528b",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "#3b528b",
  fill_vert       = "#3b528b",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "#3b528b",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot_all$patchwork

# FUNCTIONAL DIVERTITY INDICES  ----
### Functional alpha diversity indices in a multidimensional space ----

#Compute alpha FD indexes:
alpha_fd_indices_species_ecoreg <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_species[ , c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w         = s_base_ecoreg,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE, 
  details_returned = TRUE)

##Plot functional indexes ----
plots_alpha_ecoreg <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices_species_ecoreg,
  plot_asb_nm              = c("Aleutian Islands _s202200", "Aleutian Islands _s210045"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#009392FF", asb2 = "#D0587EFF"), 
  color_vert               = c(pool = "grey50", asb1 = "#009392FF", asb2 = "#D0587EFF"),
  fill_sp                  = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
  fill_vert                = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
  color_ch                 = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
  fill_ch                  = c(pool = "white", asb1 = "#009392FF", asb2 = "#D0587EFF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

#Now I plot the indexes one by one
#FUNCTIONAL RICHNESS
plots_alpha_ecoreg$"fric"$"patchwork"

#Functional divergence
plots_alpha_ecoreg$"fdiv"$"patchwork"

#Functional specialization
plots_alpha_ecoreg$"fspe"$"patchwork"

#Functional dispersion
plots_alpha_ecoreg$"fdis"$"patchwork"

#Functional identity
plots_alpha_ecoreg$"fide"$"patchwork"

#Functional evenness
plots_alpha_ecoreg$"feve"$"patchwork"

#Functional originality
plots_alpha_ecoreg$"fori"$"patchwork"

#Functional Nearest Neighbour Distance
plots_alpha_ecoreg$"fnnd"$"patchwork"

### PLOT ECOREGIONS AXIS AND EXPLORE SPECIES ----

##### Create the function ----
library(readr)
loop_baseline_scenario <- read_delim("test_127x15_ECOREG_221103_abund/loop_baseline_scenario.csv", 
                                     delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(loop_baseline_scenario)

alpha_fd_indices_species_ecoreg[["functional_diversity_indices"]] 

plot_alpha_ecoregions <- function() {
  
  for(i in 1:nrow(loop_baseline_scenario)){
    baseline <- loop_baseline_scenario$baseline[i]
    scenario <- loop_baseline_scenario$scenario[i]
  
    plot <- mFD::alpha.multidim.plot(
         output_alpha_fd_multidim = alpha_fd_indices_species_ecoreg,
          plot_asb_nm              = c(baseline, scenario),
          ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                                       "fdiv", "fori", "fspe"),
          faxes                    = NULL,
          faxes_nm                 = NULL,
          range_faxes              = c(NA, NA),
          color_bg                 = "grey95",
          shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
          size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
          color_sp                 = c(pool = "grey50", asb1 = "#009392FF", asb2 = "#D0587EFF"), 
          color_vert               = c(pool = "grey50", asb1 = "#009392FF", asb2 = "#D0587EFF"),
          fill_sp                  = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
          fill_vert                = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
          color_ch                 = c(pool = NA, asb1 = "#009392FF", asb2 = "#D0587EFF"),
          fill_ch                  = c(pool = "white", asb1 = "#009392FF", asb2 = "#D0587EFF"),
          alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
          shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
          shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
          shape_centroid_fspe      = 23,
          color_centroid_fspe      = "black",
          size_sp_nm               = 3, 
          color_sp_nm              = "black",
          plot_sp_nm               = NULL,
          fontface_sp_nm           = "plain",
          save_file                = TRUE,
          check_input              = TRUE) 
  }
  return(plot)
}













