#SPECIES TURNOVER ----
#Gissi et al., version at Sept. 13, 2023
#to understand how species change by ecoregion and 
#species contribution to change

rm( list=ls())

#set working directory ----
setwd("~/MEDIX")

library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(mFD)
library(plotly)
library(vegan)
library(readr)
library(tibble)
library(forcats) 
library(data.table)
library(ggplot2)
library(rgdal)
library(rgeos)
library(maps)
library(maptools)
library(raster)
library(viridis)
library(colourpicker)
library("PerformanceAnalytics") #correlation plots
library("ggpubr") #to plot correlation
library("patchwork") #to combine graphs
library ("ggordiplots")
library("indicspecies")




#I start from the input files which are the following

species_full_list <-
  read_delim(
    "gissi_speciesoccursum_FINAL.csv",
    delim = ";",
    escape_double = FALSE,
    trim_ws = TRUE
  ) %>%
  rename(
    speciesid = SpeciesID,
    aqmspeccode = SpecCode,
    genus = Genus,
    species = Species,
    fbname = FBname,
    occurcells = OccurCells,
    kingdom = Kingdom,
    phylum = Phylum,
    class = Class,
    order = Order,
    family = Family
  )  %>% 
  mutate(species_name=(paste(genus,species, sep=" ")))


library(readxl)
species_order_class <- read_excel("species_order_class_ORDERING.xlsx")
View(species_order_class_ORDERING)


species_abund_by_ecoreg <- read_csv("species_abund_by_ecoreg.csv")

species_abund_by_ecoreg2 <- species_abund_by_ecoreg %>%   
  mutate(species_abund_by_ecoreg, ecoreg= as.character(purrr::map(strsplit(s_ecoreg, split = "_"), 1))) %>% 
  filter(
    ecoreg == "Aleutian Islands " |
      ecoreg == "Cortezian " |
      ecoreg == "Eastern Bering Sea " |
      ecoreg == "Gulf of Alaska " |
      ecoreg == "Hawaii " |
      ecoreg == "Magdalena Transition " |
      ecoreg == "North American Pacific Fijordland " |
      ecoreg == "Northern California " |
      ecoreg == "Oregon Washington Vancouver Coast and Shelf " |
      ecoreg == "Puget Trough Georgia Basin " |
      ecoreg == "Revillagigedos " |
      ecoreg == "Aleutian Islands " |
      ecoreg == "Southern California Bight ") %>% 
  dplyr::select(!ecoreg) %>% 
  dplyr::select(!...1) 


#I prepare the input files to compare scenarios
n_last<-7
species_abund_scenarios <- species_abund_by_ecoreg2 %>% 
  mutate(scen= substr(s_ecoreg, nchar(s_ecoreg) - n_last + 1, nchar(s_ecoreg))) %>% 
  mutate(ecoreg= as.character(purrr::map(strsplit(s_ecoreg, split = "_"), 1))) %>% 
  mutate(ecoreg = fct_relevel(
    ecoreg,
    "Aleutian Islands ",
    "Eastern Bering Sea ",
    "Gulf of Alaska ",
    "North American Pacific Fijordland ",
    "Puget Trough Georgia Basin ",
    "Oregon Washington Vancouver Coast and Shelf ",
    "Northern California ",
    "Southern California Bight ",
    "Cortezian ",
    "Magdalena Transition ",
    "Revillagigedos ",
    "Hawaii "
  ))

species_abund_s202200 <- species_abund_scenarios %>% 
  filter(scen=="s202200") %>% 
  mutate(ecoreg = fct_relevel(
    ecoreg,
    "Aleutian Islands ",
    "Eastern Bering Sea ",
    "Gulf of Alaska ",
    "North American Pacific Fijordland ",
    "Puget Trough Georgia Basin ",
    "Oregon Washington Vancouver Coast and Shelf ",
    "Northern California ",
    "Southern California Bight ",
    "Cortezian ",
    "Magdalena Transition ",
    "Revillagigedos ",
    "Hawaii "
  )) %>% 
  arrange(ecoreg) %>% 
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

species_abund_s210026 <- species_abund_scenarios %>% 
  filter(scen=="s210026") %>% 
  arrange(ecoreg) %>% 
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

species_abund_s205045 <- species_abund_scenarios %>% 
  filter(scen=="s205045") %>% 
  arrange(ecoreg) %>%
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

species_abund_s205085 <- species_abund_scenarios %>% 
  filter(scen=="s205085") %>% 
  arrange(ecoreg) %>%
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

species_abund_s210045 <- species_abund_scenarios %>% 
  filter(scen=="s210045") %>% 
  arrange(ecoreg) %>%
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

species_abund_s210085 <- species_abund_scenarios %>% 
  filter(scen=="s210085") %>% 
  arrange(ecoreg) %>%
  column_to_rownames( var ="ecoreg" ) %>% 
  dplyr::select(-c(s_ecoreg,scen))

ecoreg_column <- species_abund_scenarios %>% 
  filter(scen=="s202200") %>% 
  arrange(ecoreg) %>%
  dplyr::select(-c(s_ecoreg,scen))

#Now I apply TBI
install.packages("adespatial")
install.packages("codyn")
library("adespatial")
library("codyn")


####TBI with abundance ----

#TBI(mat1, mat2, method = "%difference", pa.tr = FALSE, nperm = 99,
# BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
#save.BC = FALSE, seed. = NULL, clock = FALSE)

tbi_abund_s210026 <- TBI(species_abund_s202200, species_abund_s210026, method = "%difference", pa.tr = FALSE, nperm = 999,
                         BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
                         save.BC = TRUE, seed. = NULL, clock = TRUE)

tbi_abund_s210026v <- tbi_abund_s210026[["TBI"]] 


tbi_abund_s210045 <- TBI(species_abund_s202200, species_abund_s210045, method = "%difference", pa.tr = FALSE, nperm = 999,
                         BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
                         save.BC = TRUE, seed. = NULL, clock = TRUE)

tbi_abund_s210045v <- tbi_abund_s210045[["TBI"]]

tbi_abund_s205045 <- TBI(species_abund_s202200, species_abund_s205045, method = "%difference", pa.tr = FALSE, nperm = 999,
                         BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
                         save.BC = TRUE, seed. = NULL, clock = TRUE)

tbi_abund_s205045v <- tbi_abund_s205045[["TBI"]]

tbi_abund_s210085 <- TBI(species_abund_s202200, species_abund_s210085, method = "%difference", pa.tr = FALSE, nperm = 999,
                         BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
                         save.BC = TRUE, seed. = NULL, clock = TRUE)

tbi_abund_s210085v <- tbi_abund_s210085[["TBI"]]

tbi_abund_s205085 <- TBI(species_abund_s202200, species_abund_s205085, method = "%difference", pa.tr = FALSE, nperm = 999,
                         BCD = TRUE, replace = FALSE, test.BC = TRUE, test.t.perm = FALSE,
                         save.BC = TRUE, seed. = NULL, clock = TRUE)

tbi_abund_s205085v <- tbi_abund_s205085[["TBI"]]



#plot b-c    
prova_plot <- plot(
  tbi_abund_s205045,
  type = "BC",
  s.names = NULL,
  pch.loss = 21,
  pch.gain = 22,
  cex.names = 1,
  col.rim = "black",
  col.bg = "gold1",
  cex.symb = 3,
  diam = TRUE,
  main = "B-C plot",
  cex.main = 1,
  cex.lab = 1,
  xlim = NULL,
  ylim = NULL,
  silent = TRUE
)

ggplotly(prova_plot)

bc_abund_s205045v <- tbi_abund_s205045[["BC"]]
bcsummary_abund_s205045v <- tbi_abund_s205045[["BCD.summary"]]

#BCD.mat gives back the percentage difference: denominator = (2A+B+C)
#B=loss of species
#C=gain of species
#D=percent difference or Bray Curtis dissimilarity

bcdmat_abund_s205045v <- tbi_abund_s205045[["BCD.mat"]]
bcdmat_abund_s205085v <- tbi_abund_s205085[["BCD.mat"]]
bcdmat_abund_s210045v <- tbi_abund_s210045[["BCD.mat"]]
bcdmat_abund_s210085v <- tbi_abund_s210085[["BCD.mat"]]
bcdmat_abund_s210026v <- tbi_abund_s210026[["BCD.mat"]]

bcdmat_abund_s205045v <- bcdmat_abund_s205045v %>% 
  rename(
    fb_s205045='B/(2A+B+C)',
    fc_s205045='C/(2A+B+C)',
    fd_s205045='D=(B+C)/(2A+B+C)',
    fchange_s205045='Change' 
  )
bcdmat_abund_s210045v <- bcdmat_abund_s210045v %>% 
  rename(
    fb_s210045='B/(2A+B+C)',
    fc_s210045='C/(2A+B+C)',
    fd_s210045='D=(B+C)/(2A+B+C)',
    fchange_s210045='Change' 
  )
bcdmat_abund_s205085v <- bcdmat_abund_s205085v %>% 
  rename(
    fb_s205085='B/(2A+B+C)',
    fc_s205085='C/(2A+B+C)',
    fd_s205085='D=(B+C)/(2A+B+C)',
    fchange_s205085='Change' 
  )
bcdmat_abund_s210085v <- bcdmat_abund_s210085v %>% 
  rename(
    fb_s210085='B/(2A+B+C)',
    fc_s210085='C/(2A+B+C)',
    fd_s210085='D=(B+C)/(2A+B+C)',
    fchange_s210085='Change' 
  )
bcdmat_abund_s210026v <- bcdmat_abund_s210026v %>% 
  rename(
    fb_s210026='B/(2A+B+C)',
    fc_s210026='C/(2A+B+C)',
    fd_s210026='D=(B+C)/(2A+B+C)',
    fchange_s210026='Change' 
  )



# Create a data frame with the vectors as columns
tbi_abund_tot <- data.frame(
  ecoreg_column$ecoreg,
  tbi_abund_s210026v,
  tbi_abund_s205045v,
  tbi_abund_s210045v,
  tbi_abund_s205085v,
  tbi_abund_s210085v
) %>% 
  rename(
    ftbi_abund_s210026=tbi_abund_s210026v,
    ftbi_abund_s205045=tbi_abund_s205045v,
    ftbi_abund_s210045=tbi_abund_s210045v,
    ftbi_abund_s205085=tbi_abund_s205085v,
    ftbi_abund_s210085=tbi_abund_s210085v
    
  )


bcdmat_abund_tot <- data.frame(
  ecoreg_column$ecoreg,
  bcdmat_abund_s210026v,
  bcdmat_abund_s205045v,
  bcdmat_abund_s210045v,
  bcdmat_abund_s205085v,
  bcdmat_abund_s210085v)


#Pivot_long to plot - I start from here without doing all the calculations

library(readr)

n_last<-7
n_last_rcp<-2

bcdmat_abund_tot_long <- bcdmat_abund_tot %>% 
  rename(ecoreg=ecoreg_column.ecoreg) %>% 
  dplyr::select(-fchange_s210026, -fchange_s205045, -fchange_s210045, -fchange_s205085, -fchange_s210085) %>% 
  pivot_longer(
    cols = starts_with("f"),
    names_to = "tbi_index_scen",
    names_prefix = "f",
    values_to = "tbi",
    values_drop_na = TRUE
  ) %>% 
  mutate(scen= substr(tbi_index_scen, nchar(tbi_index_scen) - n_last + 1, nchar(tbi_index_scen))) %>% 
  mutate(s_ecoreg = paste(ecoreg, scen, sep=" _")) %>% 
  mutate(tbi_index = str_sub(tbi_index_scen, 1, 1)) %>% 
  mutate(rcp= substr(tbi_index_scen, nchar(tbi_index_scen) - n_last_rcp + 1, nchar(tbi_index_scen))) %>% 
  mutate(year = str_sub(tbi_index_scen, 4, 7)) 

bcd_change_abund_tot_long <- bcdmat_abund_tot.csv %>% 
  rename(ecoreg=ecoreg_column.ecoreg) %>% 
  dplyr::select(ecoreg, fchange_s210026,  fchange_s205045, fchange_s210045, fchange_s205085, fchange_s210085) %>% 
  pivot_longer(
    cols = starts_with("f"),
    names_to = "tbi_change_scen",
    names_prefix = "f",
    values_to = "change_value",
    values_drop_na = TRUE
  ) %>% 
  mutate(scen= substr(tbi_change_scen, nchar(tbi_change_scen) - n_last + 1, nchar(tbi_change_scen))) %>% 
  mutate(s_ecoreg = paste(ecoreg, scen, sep=" _")) 


tbi_abund_tot_long <- tbi_abund_tot.csv %>% 
  rename(ecoreg=ecoreg_column.ecoreg) %>% 
  pivot_longer(
    cols = starts_with("f"),
    names_to = "tbi_index_scen",
    names_prefix = "f",
    values_to = "tbi",
    values_drop_na = TRUE
  ) %>% 
  mutate(scen= substr(tbi_index_scen, nchar(tbi_index_scen) - n_last + 1, nchar(tbi_index_scen))) %>% 
  mutate(s_ecoreg = paste(ecoreg, scen, sep=" _")) %>% 
  mutate(tbi_index = str_sub(tbi_index_scen, 1, 3)) %>% 
  mutate(rcp= substr(tbi_index_scen, nchar(tbi_index_scen) - n_last_rcp + 1, nchar(tbi_index_scen))) %>% 
  mutate(year = str_sub(tbi_index_scen, 12, 15))  

tbi_bcdmat_abund_tot_long <- bcdmat_abund_tot_long %>% 
  bind_rows(tbi_abund_tot_long)


ecoreg_centroids_latlong <- ecoreg_centroids %>% 
  dplyr::select(ecoreg, centroid_lat, centroid_long)


bcdmat_abund_tot_long2 <- bcdmat_abund_tot_long %>% 
  left_join(ecoreg_centroids_latlong, by="ecoreg") %>% 
  group_by(rcp, year, ecoreg, tbi_index, centroid_lat, centroid_long) 
#summarize(tbi_mean=mean(tbi))

tbi_bcdmat_abund_tot_long2 <- tbi_bcdmat_abund_tot_long %>% 
  left_join(ecoreg_centroids_latlong, by="ecoreg") %>% 
  group_by(rcp, year, ecoreg, tbi_index, centroid_lat, centroid_long)
#summarize(tbi_mean=mean(tbi))



##### plot TBI ----

#TBI WITH LOSS AND GAIN AND TURNOVER AT YEAR 2100
p_tbi_solo_byrcp_vertical_2100 <- tbi_bcdmat_abund_tot_long2 %>%
  filter (year=="2100") %>% 
  filter (tbi_index=="b" | tbi_index=="c" |tbi_index=="tbi") %>% 
  ggplot(aes(x=tbi, y=centroid_lat , group=tbi_index)) +
  geom_point(aes(color=tbi_index), shape=17, size = 4) + #color=concordance
  xlim(0,0.4) +
  scale_color_manual(values=c("#440154", "#fde725", "#808080"))+
  theme_bw() +             
  ylab("Latitude") +   
  xlab("Species turnover") +
  facet_wrap(~rcp, ncol = 3, labeller = labeller(rcp = rcp.labs)) +
  geom_vline(xintercept= 0, linetype='solid', color='black') 

p_tbi_solo_byrcp_vertical_2100





