#Gissi et al. 
#version September 13, 2023
#this file is to plot the functional indices by ecoregion 
# from the results of the analysis of functional indices

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




# Species distribution from Aquamaps for all the scenarios ----

scen_2022_current <- read_csv("species_occur.csv") %>% 
  filter (probability >= 0.5)
View(scen_2022_current)

scen_2022_current_abund <- scen_2022_current %>%   
  group_by (sci_name) %>% 
  mutate (abund = probability / sum(probability) * 10000)

# GROUP BY ECOREGION  
s_base_ecoreg <- scen_2022_current_abund %>% 
  dplyr::select(s_ecoreg, sci_name, abund) %>% 
  tidyr::pivot_wider(names_from = sci_name, values_from = abund, values_fn=sum, values_fill = NA) %>%   
  filter(s_ecoreg != "Clipperton _s210085") 
s_base_ecoreg[is.na(s_base_ecoreg)] <- 0

rowSums(s_base_ecoreg != 0) 
less_six <- s_base_ecoreg[which(rowSums(s_base_ecoreg != 0) < 6), ]
rowSums(less_six != 0)
nrow(less_six)


# FD PLOT BY ECOREGION ----

#Simpson's diversity index ----
s_base_spdiv <- scen_2022_current_abund %>% 
  dplyr::select(s_ecoreg, sci_name, abund) %>% 
  group_by(s_ecoreg, sci_name) %>% 
  summarise (abund_ecoreg=sum(abund))

ecoreg_sp_div <- s_base_spdiv %>% 
  group_by(s_ecoreg) %>%
  filter(abund_ecoreg>0) %>%
  summarise(N=sum(abund_ecoreg),
            shannon.di=-sum((abund_ecoreg/sum(abund_ecoreg))*log(abund_ecoreg/sum(abund_ecoreg))),
            simpson.di=1-sum((abund_ecoreg/sum(abund_ecoreg))^2),
            inv.simpson.di=1/sum((abund_ecoreg/sum(abund_ecoreg))^2)) %>%
  arrange(-shannon.di)
ecoreg_sp_div

#dataframe with the functional diversity indices from previous script:
fd_ind_ecoreg <- read_delim("FDindexes_ecoreg.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 
colnames(fd_ind_ecoreg)[1] <- "s_ecoreg"

fd_ind_ecoreg_tot <- fd_ind_ecoreg %>% 
  left_join(ecoreg_sp_div, "s_ecoreg")

library(tibble)
n_last<-7
fd_ind_ecoreg_tot$sp_richn <- as.numeric(fd_ind_ecoreg_tot$sp_richn)
str(fd_ind_ecoreg_tot)

fd_ind_ecoreg_tot2 <- fd_ind_ecoreg_tot %>%  #dataframe for calculation
  mutate(scen= substr(s_ecoreg, nchar(s_ecoreg) - n_last + 1, nchar(s_ecoreg))) %>% 
  mutate(fd_ind_ecoreg_tot, ecoreg= as.character(purrr::map(strsplit(s_ecoreg, split = "_"), 1))) %>% 
  mutate(fsp_richn = sp_richn/127) %>% 
  mutate(fredund = (simpson.di - fmpd)) #redundancy 

colnames(fd_ind_ecoreg_tot2)[which(names(fd_ind_ecoreg_tot2) == "shannon.di")] <- "fshandi" 
colnames(fd_ind_ecoreg_tot2)[which(names(fd_ind_ecoreg_tot2) == "simpson.di")] <- "fsimpdi" 



fd_ind_ecoreg_tot2_selection <- fd_ind_ecoreg_tot2 %>% 
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
      ecoreg == "Southern California Bight ")


#prepare dataframe to plot
fd_ind_ecoreg_long <- fd_ind_ecoreg_tot2 %>% 
  pivot_longer(
    cols = starts_with("f"),
    names_to = "indices",
    names_prefix = "f",
    values_to = "fdi",
    values_drop_na = TRUE
  )
fd_ind_ecoreg_long$fdi <- as.numeric(fd_ind_ecoreg_long$fdi)

#filter on the 12 ecoregions of the study
library(forcats)

fd_ind_ecoreg_long3 <- fd_ind_ecoreg_long %>%
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
      ecoreg == "Southern California Bight "
  ) %>%
  mutate(
    ecoreg = fct_relevel(
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
    )
  )  



### FRIC ABSOLUTE DIFFERENCE ----

n_last_rcp = 2
fd_ind_ecoreg_diff_fric <- fd_ind_ecoreg_tot2 %>% 
  mutate (rcp = substr(scen, nchar(scen) - n_last_rcp + 1, nchar(scen))) %>% 
  mutate (year = substr(scen, 1, 5)) %>%
    dplyr::select(ecoreg, rcp, year, fric) %>%
  pivot_wider(names_from = year, values_from = fric)

write_csv(fd_ind_ecoreg_diff_fric, file = " fd_ind_ecoreg_diff_fric2.csv")

#upload
fd_ind_ecoreg_diff_fric2 <- read_delim(" fd_ind_ecoreg_diff_fric_mod.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_fric3 <-  fd_ind_ecoreg_diff_fric2 %>% 
  mutate(diff_2050=s2050-s2022, diff_2100=s2100-s2022) %>% 
  rename(diff_2022 = s2022)

#no change at the baseline
fd_ind_ecoreg_diff_fric3$diff_2022 <- 0

n_last_year <- 4
fd_ind_ecoreg_diff_fric4 <-  fd_ind_ecoreg_diff_fric3 %>%
  dplyr::select(ecoreg, rcp, diff_2022, diff_2050, diff_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fric_diff",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot
fd_ind_ecoreg_diff_fric5 <-  fd_ind_ecoreg_diff_fric4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) 




### SHANNON DIVERSITY ABSOLUTE DIFFERENCE ----
n_last_rcp = 2
fd_ind_ecoreg_diff_fshandi <- fd_ind_ecoreg_tot2 %>% 
  mutate (rcp = substr(scen, nchar(scen) - n_last_rcp + 1, nchar(scen))) %>% 
  mutate (year = substr(scen, 1, 5)) %>%
  dplyr::select(ecoreg, rcp, year, fshandi) %>%
  pivot_wider(names_from = year, values_from = fshandi)

#save
write_csv(fd_ind_ecoreg_diff_fshandi, file = " fd_ind_ecoreg_diff_fshandi.csv")

#upload
fd_ind_ecoreg_diff_fshandi2 <- read_delim(" fd_ind_ecoreg_diff_fshandi_mod.csv", 
                                            delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_fshandi3 <-  fd_ind_ecoreg_diff_fshandi2 %>% 
  mutate(diff_2050=s2050-s2022, diff_2100=s2100-s2022) %>% 
  rename(diff_2022 = s2022)

#no change with the baseline
fd_ind_ecoreg_diff_fshandi3$diff_2022 <- 0

n_last_year <- 4

fd_ind_ecoreg_diff_fshandi4 <-  fd_ind_ecoreg_diff_fshandi3 %>%
  dplyr::select(ecoreg, rcp, diff_2022, diff_2050, diff_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fshandi_diff",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot
fd_ind_ecoreg_diff_fshandi5 <-  fd_ind_ecoreg_diff_fshandi4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) 


### FREDUND ABSOLUTE DIFFERENCE ----

n_last_rcp = 2
fd_ind_ecoreg_diff_fredund <- fd_ind_ecoreg_tot2 %>% 
  mutate (rcp = substr(scen, nchar(scen) - n_last_rcp + 1, nchar(scen))) %>% 
  mutate (year = substr(scen, 1, 5)) %>%
  dplyr::select(ecoreg, rcp, year, fredund) %>%
  pivot_wider(names_from = year, values_from = fredund)

#save and add baseline values to calculate change
write_csv(fd_ind_ecoreg_diff_fredund, file = " fd_ind_ecoreg_diff_fredund2.csv")

#upload
fd_ind_ecoreg_diff_fredund2 <- read_delim(" fd_ind_ecoreg_diff_fredund_mod.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_fredund3 <-  fd_ind_ecoreg_diff_fredund2 %>% 
  mutate(diff_2050=s2050-s2022, diff_2100=s2100-s2022) %>% 
  rename(diff_2022 = s2022)

#no change with baseline
fd_ind_ecoreg_diff_fredund3$diff_2022 <- 0

n_last_year <- 4
fd_ind_ecoreg_diff_fredund4 <-  fd_ind_ecoreg_diff_fredund3 %>%
  dplyr::select(ecoreg, rcp, diff_2022, diff_2050, diff_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fredund_diff",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot 
fd_ind_ecoreg_diff_fredund5 <-  fd_ind_ecoreg_diff_fredund4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) 

 
### FRIC RELATIVE DIFFERENCE ----

#upload file
fd_ind_ecoreg_diff_fric2 <- read_delim(" fd_ind_ecoreg_diff_fric_mod.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_perc_fric3 <-  fd_ind_ecoreg_diff_fric2 %>%
  mutate(                         #absolute difference
    diff_2022 = s2022 - s2022,
    diff_2050 = s2050 - s2022,
    diff_2100 = s2100 - s2022
  ) %>%
  mutate(                         #percent change
    diff_perc_2022 = diff_2022 / s2022 * 100,
    diff_perc_2050 = diff_2050 / s2022 * 100,
    diff_perc_2100 = diff_2100 / s2022 * 100
  )

#dataframe with percent change
n_last_year <- 4
fd_ind_ecoreg_diff_perc_fric4 <-  fd_ind_ecoreg_diff_perc_fric3 %>%
  dplyr::select(ecoreg, rcp, diff_perc_2022, diff_perc_2050, diff_perc_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fric_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot
fd_ind_ecoreg_diff_perc_fric5 <-  fd_ind_ecoreg_diff_perc_fric4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) 

# New facet label names for RCP variable
rcp.labs <- c("RCP 2.6", "RCP 4.5", "RCP 8.5")
names(rcp.labs) <- c("26", "45", "85")



#####FRIC RELATIVE DIFFERENCE BY LATITUDE ----

ecoreg_centroids <- read_csv(" ecoreg_centroids.csv")

#NO OUTLIERS
fd_ind_ecoreg_diff_perc_fric6_lat <- fd_ind_ecoreg_diff_perc_fric5 %>% 
  left_join(ecoreg_centroids, by='ecoreg')

p_fric_diff_perc_bylat <- fd_ind_ecoreg_diff_perc_fric6_lat %>%
  ggplot(aes(x=fric_diff_perc, y=centroid_lat, group=year)) +
  geom_point(aes(shape=year, color=year), size = 3) +
  scale_shape_manual(values=c(17, 17, 17))+
  scale_color_manual(values=c('#C9B0CB', '#A37EA9', '#440154FF'))+
  scale_size_manual(values=c(20,20,20))+
  xlim(-50,75) +
  theme_bw() +             # Changing the theme to get rid of the grey background
  ylab("Latitude") +   
  xlab("Functional richness percent change (∆%)") +
  #scale_y_discrete(labels=c('AI', 'EBS', 'GA', 'NAPF', 'PTGB', 'OWVCS', 'NC', 'SCB', 'C', 'MT', 'R', 'H'))+
  facet_wrap(~rcp, ncol = 3, labeller = labeller(rcp = rcp.labs)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Ecoregion",  #funziona! per includere 
                                         breaks = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$centroid_lat),
                                         labels = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$ecoreg)))

p_fric_diff_perc_bylat


#COMPLETE FULL WITH OUTLIERS
p_fric_diff_perc_bylat_full <- fd_ind_ecoreg_diff_perc_fric6_lat %>%
  ggplot(aes(x=fric_diff_perc, y=centroid_lat, group=year)) +
  geom_point(aes(shape=year, color=year), size = 3) +
  scale_shape_manual(values=c(17, 17, 17))+
  scale_color_manual(values=c('#C9B0CB', '#A37EA9', '#440154FF'))+
  scale_size_manual(values=c(20,20,20))+
 # xlim(-50,75) +
  theme_bw() +             # Changing the theme to get rid of the grey background
  ylab("Latitude") +   
  xlab("Functional richness percent change (∆%)") +
  #scale_y_discrete(labels=c('AI', 'EBS', 'GA', 'NAPF', 'PTGB', 'OWVCS', 'NC', 'SCB', 'C', 'MT', 'R', 'H'))+
  facet_wrap(~rcp, ncol = 3, labeller = labeller(rcp = rcp.labs)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Ecoregion",  #funziona! per includere 
                                         breaks = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$centroid_lat),
                                         labels = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$ecoreg)))

p_fric_diff_perc_bylat_full 




### SHANNON INDEX RELATIVE DIFFERENCE ----

#dataframe
fd_ind_ecoreg_diff_fshandi2 <- read_delim(" fd_ind_ecoreg_diff_fshandi_mod.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_perc_fshandi3 <-  fd_ind_ecoreg_diff_fshandi2 %>%
  mutate(                        
    diff_2022 = s2022 - s2022,
    diff_2050 = s2050 - s2022,
    diff_2100 = s2100 - s2022
  ) %>%
  mutate(                        
    diff_perc_2022 = diff_2022 / s2022 * 100,
    diff_perc_2050 = diff_2050 / s2022 * 100,
    diff_perc_2100 = diff_2100 / s2022 * 100
  )

#dataframe with changes
n_last_year <- 4
fd_ind_ecoreg_diff_perc_fshandi4 <-  fd_ind_ecoreg_diff_perc_fshandi3 %>%
  dplyr::select(ecoreg, rcp, diff_perc_2022, diff_perc_2050, diff_perc_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fshandi_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot
fd_ind_ecoreg_diff_perc_fshandi5 <-  fd_ind_ecoreg_diff_perc_fshandi4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) #questa parte mi serve per dare l'ordine ai grafici


#####SHANNON INDEX RELATIVE DIFFERENCE BY LATITUDE ----

#complete
fd_ind_ecoreg_diff_perc_fshandi6_lat <- fd_ind_ecoreg_diff_perc_fshandi5 %>% 
  left_join(ecoreg_centroids, by='ecoreg')

p_fshandi_diff_perc_bylat <- fd_ind_ecoreg_diff_perc_fshandi6_lat %>%
  ggplot(aes(x=fshandi_diff_perc, y=centroid_lat, group=year)) +
  geom_point(aes(shape=year, color=year), size = 3) +
  scale_shape_manual(values=c(17, 17, 17))+
  scale_color_manual(values=c( '#9C9FFE', '#4E53F7', '#1115A1'))+
  scale_size_manual(values=c(20,20,20))+
  #xlim(-20,15) +
  theme_bw() +             # Changing the theme to get rid of the grey background
  ylab("Latitude") +   
  xlab("Species diversity percent change (∆%)") +
  #scale_y_discrete(labels=c('AI', 'EBS', 'GA', 'NAPF', 'PTGB', 'OWVCS', 'NC', 'SCB', 'C', 'MT', 'R', 'H'))+
  facet_wrap(~rcp, ncol = 3, labeller = labeller(rcp = rcp.labs)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Ecoregion",  #funziona! per includere 
                                         breaks = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$centroid_lat),
                                         labels = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$ecoreg)))

p_fshandi_diff_perc_bylat





### FREDUND RELATIVE DIFFERENCE ----

#input
fd_ind_ecoreg_diff_fredund2 <- read_delim(" fd_ind_ecoreg_diff_fredund_mod.csv", 
                                       delim = ";", escape_double = FALSE, trim_ws = TRUE) 

fd_ind_ecoreg_diff_perc_fredund3 <-  fd_ind_ecoreg_diff_fredund2 %>%
  mutate(                         
    diff_2022 = s2022 - s2022,
    diff_2050 = s2050 - s2022,
    diff_2100 = s2100 - s2022
  ) %>%
  mutate(                         
    diff_perc_2022 = diff_2022 / s2022 * 100,
    diff_perc_2050 = diff_2050 / s2022 * 100,
    diff_perc_2100 = diff_2100 / s2022 * 100
  )

#percent change
n_last_year <- 4
fd_ind_ecoreg_diff_perc_fredund4 <-  fd_ind_ecoreg_diff_perc_fredund3 %>%
  dplyr::select(ecoreg, rcp, diff_perc_2022, diff_perc_2050, diff_perc_2100) %>% 
  pivot_longer(
    cols = starts_with("diff"),     
    names_to = "yearbad",
    names_prefix = "diff",
    values_to = "fredund_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(year= substr(yearbad, nchar(yearbad) - n_last_year + 1, nchar(yearbad))) %>% 
  dplyr::select(!yearbad) 

#order to plot
fd_ind_ecoreg_diff_perc_fredund5 <-  fd_ind_ecoreg_diff_perc_fredund4 %>% 
  filter( ecoreg == "Aleutian Islands" |
            ecoreg == "Cortezian" |
            ecoreg == "Eastern Bering Sea" |
            ecoreg == "Gulf of Alaska" |
            ecoreg == "Hawaii" |
            ecoreg == "Magdalena Transition" |
            ecoreg == "North American Pacific Fijordland" |
            ecoreg == "Northern California" |
            ecoreg == "Oregon Washington Vancouver Coast and Shelf" |
            ecoreg == "Puget Trough Georgia Basin" |
            ecoreg == "Revillagigedos" |
            ecoreg == "Aleutian Islands" |
            ecoreg == "Southern California Bight"
  ) %>% 
  mutate(ecoreg = fct_relevel(ecoreg, "Aleutian Islands", "Eastern Bering Sea", 
                              "Gulf of Alaska", "North American Pacific Fijordland",
                              "Puget Trough Georgia Basin", "Oregon Washington Vancouver Coast and Shelf",
                              "Northern California", "Southern California Bight",
                              "Cortezian", "Magdalena Transition", "Revillagigedos",
                              "Hawaii")) 

# New facet label names for RCP variable
rcp.labs <- c("RCP 2.6", "RCP 4.5", "RCP 8.5")
names(rcp.labs) <- c("26", "45", "85")



#####FREDUNDANCY RELATIVE DIFFERENCE BY LATITUDE ----

#complete
fd_ind_ecoreg_diff_perc_fredund6_lat <- fd_ind_ecoreg_diff_perc_fredund5 %>% 
  left_join(ecoreg_centroids, by='ecoreg')

p_fredund_diff_perc_bylat <- fd_ind_ecoreg_diff_perc_fredund6_lat %>%
  ggplot(aes(x=fredund_diff_perc, y=centroid_lat, group=year)) +
  geom_point(aes(shape=year, color=year), size = 3) +
  scale_shape_manual(values=c(17, 17, 17))+
  scale_color_manual(values=c('#E478FC', '#C034DF', '#730D89'))+
  scale_size_manual(values=c(20,20,20))+
  theme_bw() +             
  ylab("Latitude") +   
  xlab("Functional redundancy percent change (∆%)") +
  facet_wrap(~rcp, ncol = 3, labeller = labeller(rcp = rcp.labs)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_y_continuous(sec.axis = sec_axis(~ ., name = "Ecoregion",  
                                         breaks = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$centroid_lat),
                                         labels = unique(fd_ind_ecoreg_diff_perc_fshandi6_lat$ecoreg)))

p_fredund_diff_perc_bylat
ggsave("16_fredund_diff_perc_bylat_nooutliers.tiff", units="in", width=8, height=3, dpi=300, compression = 'lzw')


fredund_diff_perc <- fd_ind_ecoreg_diff_perc_fredund5 %>% 
  filter(year!="2022") %>% 
  mutate(scen= paste("s",year,rcp, sep="")) %>%  
  mutate (s_ecoreg = paste(ecoreg,scen,sep= " _")) %>% 
  dplyr::select(s_ecoreg, fredund_diff_perc)



### Join dataframes ----

fd_ind_ecoreg_diff_perc_fsp_richn6 <- fd_ind_ecoreg_diff_perc_fsp_richn5 %>%
  mutate(ecoreg_rcp_year = paste(ecoreg, rcp, year, sep="_"))

fd_ind_ecoreg_diff_perc_feve6 <- fd_ind_ecoreg_diff_perc_feve5 %>%
  mutate(ecoreg_rcp_year = paste(ecoreg, rcp, year, sep="_")) %>%
  dplyr::select(ecoreg_rcp_year, feve_diff_perc)

fd_ind_ecoreg_diff_perc_fric6 <- fd_ind_ecoreg_diff_perc_fric5 %>% 
  mutate(ecoreg_rcp_year = paste(ecoreg, rcp, year, sep="_")) %>%
  dplyr::select(ecoreg_rcp_year, fric_diff_perc)

fd_ind_ecoreg_diff_perc_fshandi6 <- fd_ind_ecoreg_diff_perc_fshandi5 %>% 
  mutate(ecoreg_rcp_year = paste(ecoreg, rcp, year, sep="_")) %>%
  dplyr::select(ecoreg_rcp_year, fshandi_diff_perc)

fd_ing_ecoreg_diff_perc_scatterplot <- fd_ind_ecoreg_diff_perc_fsp_richn6 %>%
  left_join(fd_ind_ecoreg_diff_perc_feve6, by = "ecoreg_rcp_year") %>% 
  left_join(fd_ind_ecoreg_diff_perc_fric6, by= "ecoreg_rcp_year") %>% 
  left_join(fd_ind_ecoreg_diff_perc_fshandi6, by= "ecoreg_rcp_year")



#END ----


