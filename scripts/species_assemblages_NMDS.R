#SPECIES PATTERNS BY ECOREG NMDS ----
#Gissi et al., version at Sept. 13, 2023

#Analysis of what species drive the patterns across ecoregions
#I use the package NMDS 

rm( list=ls())

#set working directory ----
setwd("~/MEDIX")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("mFD")
install.packages("gg_ordiplot")
install.packages(c("gg_ordiplot"))

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
library(vegan)
library(readr)

## By order YEAR 2100 ----
species_abund_by_ecoreg_long <- read_csv("species_abund_by_ecoreg_long_ordered.csv")

species_abund_by_ecoreg <- read_delim("species_abund_by_ecoreg.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

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

View(species_abund_by_ecoreg2)  

species_abund_by_ecoreg3 <- species_abund_by_ecoreg2 %>% 
  column_to_rownames( var ="s_ecoreg" ) 


order_abund_by_ecoreg2100 <- species_abund_by_ecoreg_long %>% 
  dplyr::select(s_ecoreg, abund, order, ecoreg, scen, rcp, year) %>% 
  filter(scen=="s202200" | scen=="s210026" | scen=="s210045" | scen=="s210085") %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) %>%
  tidyr::pivot_wider(names_from = order, values_from = abund, values_fn=sum, values_fill = 0) %>% 
  dplyr::select(!ecoreg | !scen | !rcp | !year) %>% 
  subset(select = -c(ecoreg, scen, rcp, year)) %>% 
  column_to_rownames( var ="s_ecoreg" ) 


#Data preparation: df with species in columns and “entities” in lines (e.g., ecoregions, gridcells)

# NMDS
data.mds <- metaMDS(order_abund_by_ecoreg2100, distance = "bray")
data.mds 
data.mds$stress 
stressplot(data.mds) # Shepard plot 


# Significant Species
spp.fit <- envfit(data.mds, order_abund_by_ecoreg2100, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(sig.spp.scrs)

install.packages("ggordiplots")
remotes::install_github("jfq3/ggordiplots")
library (ggordiplots)

n_last = 7
species_abund_by_ecoreg2_2100 <- species_abund_by_ecoreg %>%   
  mutate(species_abund_by_ecoreg, ecoreg= as.character(purrr::map(strsplit(s_ecoreg, split = "_"), 1))) %>% 
  mutate(scen= substr(s_ecoreg, nchar(s_ecoreg) - n_last + 1, nchar(s_ecoreg))) %>%
  mutate (year = substr(scen, 2, 5)) %>%
  filter(scen =="s202200" | scen =="s210026" | scen =="s210045" | scen == "s210085") %>% 
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
  dplyr::select(!...1)

plot <- gg_ordiplot(data.mds, 
                    aes(shape = species_abund_by_ecoreg2_2100$year),
                    groups = species_abund_by_ecoreg2_2100$s_ecoreg,
                    hull = FALSE, ellipse = TRUE, spiders = FALSE,
                    kind = "sd", pt.size = 4)
plot

# refine the plot 
p_nmds_2100 <- plot$plot + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position= "right",
        legend.box="vertical") + 
  labs(colour= "Ecoregions at 2100 across RCPs") +
  theme_pubr(base_size = 16)+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3)+ 
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), 
                           cex = 4, direction = "both", segment.size = 0.25)+ 
  scale_colour_hue(h = c(180, 300))  

p_nmds_2100

ggplotly(p_nmds)






## By order YEAR 2050 ----

species_abund_by_ecoreg_long <- read_csv("species_abund_by_ecoreg_long_ordered.csv")


order_abund_by_ecoreg2050 <- species_abund_by_ecoreg_long %>% 
  dplyr::select(s_ecoreg, abund, order, ecoreg, scen, rcp, year) %>% 
  filter(scen=="s202200" | scen=="s205026" | scen=="s205045" | scen=="s205085") %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) %>%
  tidyr::pivot_wider(names_from = order, values_from = abund, values_fn=sum, values_fill = 0) %>% 
  dplyr::select(!ecoreg | !scen | !rcp | !year) %>% 
  subset(select = -c(ecoreg, scen, rcp, year)) %>% 
  column_to_rownames( var ="s_ecoreg" ) 


#Data preparation: df with species in columns and “entities” in lines 

# NMDS
data.mds <- metaMDS(order_abund_by_ecoreg2050, distance = "bray")
data.mds 
data.mds$stress 
stressplot(data.mds) # Shepard plot 


# Significant Species
spp.fit <- envfit(data.mds, order_abund_by_ecoreg2050, permutations = 999) # this fits species vectors
spp.scrs <- as.data.frame(scores(spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

head(sig.spp.scrs)

install.packages("ggordiplots")
remotes::install_github("jfq3/ggordiplots")
library (ggordiplots)

n_last = 7
species_abund_by_ecoreg2_2050 <- species_abund_by_ecoreg %>%   
  mutate(species_abund_by_ecoreg, ecoreg= as.character(purrr::map(strsplit(s_ecoreg, split = "_"), 1))) %>% 
  mutate(scen= substr(s_ecoreg, nchar(s_ecoreg) - n_last + 1, nchar(s_ecoreg))) %>%
  mutate (year = substr(scen, 2, 5)) %>%
  filter(scen =="s202200" | scen =="s205026" | scen =="s205045" | scen == "s205085") %>% 
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
  dplyr::select(!...1)

plot <- gg_ordiplot(data.mds, 
                    aes(shape = species_abund_by_ecoreg2_2050$year),
                    groups = species_abund_by_ecoreg2_2050$s_ecoreg,
                    hull = FALSE, ellipse = TRUE, spiders = FALSE,
                    kind = "sd", pt.size = 4)
plot



# refine the plot 
p_nmds_2050 <- plot$plot + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position= "right",
        legend.box="vertical") + 
  labs(colour= "Ecoregions at 2050 across RCPs") +
  theme_pubr(base_size = 16)+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3)+ 
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), 
                           cex = 4, direction = "both", segment.size = 0.25)+ 
  scale_colour_hue(h = c(180, 300))  



p_nmds_2050

#END

