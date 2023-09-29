#Gissi et al. - Functional diversity at grid cell level
# version September 13, 2023

## I use 127 species and 15 functional traits that 
# I have selected through exploring the correlation among them.


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


#SPECIES ABUNDANCE BY GRID CELL ----
## Species distribution from aquamaps
scen_2022_current <- read_csv("species_occur.csv") %>% 
  filter (probability >= 0.5)
View(scen_2022_current)

scen_2022_current_abund <- scen_2022_current %>%   
  group_by (sci_name) %>% 
  mutate (abund = probability / sum(probability) * 10000)

# GROUP BY CELLS   
sppass_base_cells_deg1 <- scen_2022_current_abund %>% 
  dplyr::select(s_cells_1degree, sci_name, abund) %>% 
  pivot_wider(names_from = sci_name, values_from = abund,  values_fn=sum, values_fill = NA) %>% 
  column_to_rownames( var ="s_cells_1degree" )

sppass_base_cells_deg1[is.na(sppass_base_cells_deg1)] <- 0
is.matrix(sppass_base_cells_deg1)
sppass_base_cells_deg1 <- as.matrix((sppass_base_cells_deg1))
summary(sppass_base_cells_deg1)



#DIVERSITY INDICES ----
#I upload the file with the indices calculated before
fdi_1deg <- read_csv("fdi_1deg.csv")
View(fdi_1deg) 
colnames(fdi_1deg)[1] <- "c1deg"

n_last<-7
fdi_deg1_1 <- fdi_1deg %>% 
  mutate(scen= substr(c1deg, nchar(c1deg) - n_last + 1, nchar(c1deg))) %>% 
  mutate(fdi_1deg, csquare_1deg= as.character(purrr::map(strsplit(c1deg, split = " "), 1)))

#file with coordinates
csquare_1degree <- read_delim("csquare_1degree.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  rename(csquare_1deg=c1deg)

n_last_rcp = 2
fdi_deg1_2 <- fdi_deg1_1 %>% 
  left_join(csquare_1degree, by="csquare_1deg") %>% 
  mutate(fsp_richn = sp_richn / 127) %>% 
  mutate (rcp = substr(scen, nchar(scen) - n_last_rcp + 1, nchar(scen))) %>% 
  mutate (year = substr(scen, 1, 5)) 


### F richness ----
#I prepare the file to plot the difference
n_last_rcp = 2
fdi_d1_diff_fric <- fdi_deg1_2 %>% 
  mutate (rcp = substr(scen, nchar(scen) - n_last_rcp + 1, nchar(scen))) %>% 
  mutate (year = substr(scen, 1, 5)) 


#I extract the files to make the difference betwen scenarios of the fric
fric_202200 <- fdi_d1_diff_fric %>% 
  filter(year== "s2022") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_2022=fric) 

fric_205045 <- fdi_d1_diff_fric %>% 
  filter(scen== "s205045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_205045=fric)%>% 
  dplyr::select(csquare_1deg, fric_205045)

fric_205085 <- fdi_d1_diff_fric %>% 
  filter(scen== "s205085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_205085=fric)%>% 
  dplyr::select(csquare_1deg, fric_205085)  

fric_210026 <- fdi_d1_diff_fric %>% 
  filter(scen== "s210026") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_210026=fric)%>% 
  dplyr::select(csquare_1deg, fric_210026)

fric_210045 <- fdi_d1_diff_fric %>% 
  filter(scen== "s210045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_210045=fric)%>% 
  dplyr::select(csquare_1deg, fric_210045)

fric_210085 <- fdi_d1_diff_fric %>% 
  filter(scen== "s210085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fric) %>%
  rename (fric_210085=fric)%>% 
  dplyr::select(csquare_1deg, fric_210085)

#I unify the tables to create one unique dataframe from which I make the differences
fric_1deg_diff_tot <- fric_202200 %>% 
  full_join(fric_205045, by="csquare_1deg") %>% 
  full_join(fric_205085, by="csquare_1deg") %>% 
  full_join(fric_210026, by="csquare_1deg") %>% 
  full_join(fric_210045, by="csquare_1deg") %>% 
  full_join(fric_210085, by="csquare_1deg") %>% 
  mutate (fric_diff_45_2050 = fric_205045 - fric_2022) %>% 
  mutate (fric_diff_45_2100 = fric_210045 - fric_2022) %>% 
  mutate (fric_diff_85_2050 = fric_205085 - fric_2022) %>% 
  mutate (fric_diff_85_2100 = fric_210085 - fric_2022) %>%
  mutate (fric_diff_26_2100 = fric_210026 - fric_2022)


#prepare the dataframe to plot scenarios together
n_last_year <- 4
fric_1deg_diff_tot_long <-  fric_1deg_diff_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fric_diff_45_2050, fric_diff_45_2100, 
                fric_diff_85_2050, fric_diff_85_2100,
                fric_diff_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fric_diff_"),     
    names_to = "rcp_year",
    names_prefix = "fric_diff_",
    values_to = "fric_diff",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 1, 2)) %>% 
  mutate(year = substr(rcp_year, 4, 7))

#diff per rcp 45
plot_fric_1deg_diff_45 <- fric_1deg_diff_tot_long %>% 
  filter (rcp == "45") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff)) +
  geom_raster() +
  facet_wrap(~year, ncol = 1 ) +
  scale_fill_gradientn(colours = colorspace::diverge_hcl(7), 
                       limits = c(-0.55, 0.55),
                       breaks = c(- 0.55, - 0.22, 0, 0.22, 0.55)) +
  theme_minimal()


plot_fric_1deg_diff_45

#diff per rcp 85
plot_fric_1deg_diff_85 <- fric_1deg_diff_tot_long %>% 
  filter (rcp == "85") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff)) +
  geom_raster() +
  facet_wrap(~year, ncol = 1 ) +
  scale_fill_gradientn(colours = colorspace::diverge_hcl(7), 
                       limits = c(-0.55, 0.55),
                       breaks = c(- 0.55, - 0.22, 0, 0.22, 0.55)) +
  theme_minimal()

plot_fric_1deg_diff_85




# Percent change from baseline scenario ----

## Fric % difference - 5 classes ><5% ----

n_last_year <- 4
fric_1deg_diff_perc_tot_long_5classes <-  fric_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fric_diff_perc_45_2050, fric_diff_perc_45_2100, 
                fric_diff_perc_85_2050, fric_diff_perc_85_2100,
                fric_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fric_diff_perc"),     
    names_to = "rcp_year",
    names_prefix = "fric_diff_perc",
    values_to = "fric_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fric_diff_interval = cut(fric_diff_perc,
                                  breaks = c(-Inf, -5, -0.1, 0.1, 5, Inf),
                                  labels = c("< -5", "-5 to 0", "0", "0 to 5","> 5"),
                                  include.lowest = TRUE)) 

   
##plots

plot_fric_1deg_diff_perc_205045_5classes <- fric_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fric_1deg_diff_perc_205045_5classes


plot_fric_1deg_diff_perc_205085_5classes <- fric_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fric_1deg_diff_perc_205085_5classes


plot_fric_1deg_diff_perc_210045_5classes <- fric_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fric_1deg_diff_perc_210045_5classes


plot_fric_1deg_diff_perc_210085_5classes <- fric_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fric_1deg_diff_perc_210085_5classes


plot_fric_1deg_diff_perc_210026_5classes <- fric_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_26_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fric_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fric_1deg_diff_perc_210026_5classes


#species diversity + functional richness + functional redundancy

p_1deg_fric_5classes <- ggarrange ( plot_fric_1deg_diff_perc_210026_5classes,
                                    plot_fric_1deg_diff_perc_210045_5classes,
                                    plot_fric_1deg_diff_perc_210085_5classes,
                                    plot_fric_1deg_diff_perc_205045_5classes,
                                    plot_fric_1deg_diff_perc_205085_5classes,
                                         ncol = 3,
                                         nrow = 2,
                               common.legend = TRUE,
                                         align = c("v")
)

# Display the composite plot
p_1deg_fric_5classes



## Boxplot fric by latitude ----

#First I need to connect cells with the ecoregions to define coastal by ecoregion
scen_2022_current <- read_csv("species_occur.csv") %>% 
  filter (probability >= 0.5) %>% 
  dplyr::select(cells_1degree, ecoreg, s_ecoreg)

colnames(scen_2022_current)[1] <- "csquare_1deg" 

# Remove duplicates by single column
scen_2022_current2 <- scen_2022_current[!duplicated(scen_2022_current$csquare_1deg), ]


base_gam <- read_csv("base_corr.csv")  

n_last_yearrcp=6

base_gam2 <- base_gam %>% 
  mutate(yearrcp= substr(c1deg, nchar(c1deg) - n_last_yearrcp + 1, nchar(c1deg)))

View(base_gam)

base_gam3 <- base_gam2 %>% 
  mutate(centre_longitude2=centre_longitude %% 360) %>%  
  mutate(cotang_abs = fric/fshandi) %>% 
  mutate(fredund = (fsimp - fmpd)) %>% 
  left_join(scen_2022_current2, by="csquare_1deg") %>% 
  mutate(coastal_ecoreg = if_else(ecoreg == "High seas", 
                                  "highseas", "coastal_ecoreg")) %>% 
  dplyr::select(csquare_1deg, coastal_depth, coastal_distance, ecoreg,        
                coastal_ecoreg)

# Remove duplicates by single column
base_gam4  <- base_gam3 [!duplicated(base_gam3), ]



#Now I unify the definition of coastal grid cells with the data from fric 

fric_1deg_diff_perc_tot_long_5classes_2 <- fric_1deg_diff_perc_tot_long_5classes %>% 
  left_join(base_gam4, by="csquare_1deg") %>% 
  mutate (coasttype_general = "general") %>% 
  mutate (coasttype_ecoreg = coastal_ecoreg)


#Now I plot the boxplot
plot_bp_fric_1deg_diff_perc_210026_5classes <- fric_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (rcp_year == "_26_2100") %>%
  ggplot(aes(x=fric_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()
  

plot_bp_fric_1deg_diff_perc_210026_5classes 

#plot only gridcells within ecoregions
plot_bp_coastalecoreg_fric_1deg_diff_perc_210026_5classes <- fric_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (rcp_year == "_26_2100") %>%
  filter(coastal_ecoreg == "coastal_ecoreg") %>% 
  ggplot(aes(x=fric_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()

plot_bp_coastalecoreg_fric_1deg_diff_perc_210026_5classes


# Filter the data for the year and rcp of interest
fric_1deg_diff_perc_tot_long_5classes_3 <- fric_1deg_diff_perc_tot_long_5classes_2 %>%
    pivot_longer(
    cols = starts_with("coasttype"),     
    names_to = "bar_type",
    names_prefix = "coasttype",
    values_to = "subdivision",
    values_drop_na = TRUE
  ) 

# Create the combined boxplot
combined_boxplot_fric <- fric_1deg_diff_perc_tot_long_5classes_3 %>%
  filter (year == "2100") %>% 
  filter(subdivision=="coastal_ecoreg" | subdivision=="general") %>% 
  ggplot(aes(x = fric_diff_interval, y = centre_latitude)) +
  geom_boxplot(aes(fill = subdivision), position = position_dodge(width = 0.75)) +
  theme_minimal() +
  facet_wrap(~rcp, ncol = 3) +
  ylab("Latitude") +
  xlab("F. richness change (∆%)") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("general" = "grey42", "coastal_ecoreg" = "grey80"))  

# Plot the combined boxplot
print(combined_boxplot_fric)



#Shannon's diversity index ----
install.packages("vegan")
library("vegan")

scen_2022_current_abund

sp_div <- scen_2022_current_abund %>%
  dplyr::select (sci_name, s_cells_1degree, abund) %>% 
  group_by(s_cells_1degree) %>%
  filter(abund>0) %>%
  summarise(N=sum(abund),
            shannon.di=-sum((abund/sum(abund))*log(abund/sum(abund))),
            simpson.di=1-sum((abund/sum(abund))^2),
            inv.simpson.di=1/sum((abund/sum(abund))^2)) %>%
  arrange(-shannon.di)
sp_div

#PREPARE THE DB TO MAKE THE DIFFERENCE AND ABSOLUTE DIFF OF SHANDI
#change variables
sp_div_2 <- sp_div %>% 
  rename(c1deg=s_cells_1degree, fshandi=shannon.di, fsimp=simpson.di, fsimpinv=inv.simpson.di)

sp_div_3 <- fdi_deg1_2 %>% 
  dplyr::select(c1deg, csquare_1deg, centre_latitude, centre_longitude, scen, rcp, year, fmpd) %>% 
  right_join(sp_div_2, by="c1deg") %>% 
  mutate(fredund = (fsimp - fmpd))

#### Shandi abs difference ----
#I extract the files to make the difference betwen scenarios of the fshandi
fshandi_202200 <- sp_div_3 %>% 
  filter(year== "s2022") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_2022=fshandi) 

fshandi_205045 <- sp_div_3 %>% 
  filter(scen== "s205045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_205045=fshandi)%>% 
  dplyr::select(csquare_1deg, fshandi_205045)

fshandi_205085 <- sp_div_3 %>% 
  filter(scen== "s205085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_205085=fshandi)%>% 
  dplyr::select(csquare_1deg, fshandi_205085)  

fshandi_210026 <- sp_div_3 %>% 
  filter(scen== "s210026") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_210026=fshandi)%>% 
  dplyr::select(csquare_1deg, fshandi_210026)

fshandi_210045 <- sp_div_3 %>% 
  filter(scen== "s210045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_210045=fshandi)%>% 
  dplyr::select(csquare_1deg, fshandi_210045)

fshandi_210085 <- sp_div_3 %>% 
  filter(scen== "s210085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fshandi) %>%
  rename (fshandi_210085=fshandi)%>% 
  dplyr::select(csquare_1deg, fshandi_210085)

#I unify the tables to create one unique dataframe from which I make the differences
fshandi_1deg_diff_tot <- fshandi_202200 %>% 
  full_join(fshandi_205045, by="csquare_1deg") %>% 
  full_join(fshandi_205085, by="csquare_1deg") %>% 
  full_join(fshandi_210026, by="csquare_1deg") %>% 
  full_join(fshandi_210045, by="csquare_1deg") %>% 
  full_join(fshandi_210085, by="csquare_1deg") %>% 
  mutate (fshandi_diff_45_2050 = fshandi_205045 - fshandi_2022) %>% 
  mutate (fshandi_diff_45_2100 = fshandi_210045 - fshandi_2022) %>% 
  mutate (fshandi_diff_85_2050 = fshandi_205085 - fshandi_2022) %>% 
  mutate (fshandi_diff_85_2100 = fshandi_210085 - fshandi_2022) %>%
  mutate (fshandi_diff_26_2100 = fshandi_210026 - fshandi_2022)


## Shandi % difference ----
fshandi_1deg_diff_perc_tot <- fshandi_1deg_diff_tot %>% 
  mutate (fshandi_diff_perc_45_2050 = fshandi_diff_45_2050/fshandi_2022*100) %>% 
  mutate (fshandi_diff_perc_45_2100 = fshandi_diff_45_2100/fshandi_2022*100) %>% 
  mutate (fshandi_diff_perc_85_2050 = fshandi_diff_85_2050/fshandi_2022*100) %>% 
  mutate (fshandi_diff_perc_85_2100 = fshandi_diff_85_2100/fshandi_2022*100) %>%
  mutate (fshandi_diff_perc_26_2100 = fshandi_diff_26_2100/fshandi_2022*100)

#prepare the file 
n_last_year <- 4
fshandi_1deg_diff_perc_tot_long <-  fshandi_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fshandi_diff_perc_45_2050, fshandi_diff_perc_45_2100, 
                fshandi_diff_perc_85_2050, fshandi_diff_perc_85_2100,
                fshandi_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fshandi_diff_perc"),     
    names_to = "rcp_year",
    names_prefix = "fshandi_diff_perc",
    values_to = "fshandi_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fshandi_diff_perc_range = 
           case_when (fshandi_diff_perc < -30 ~ "s.div<-30%",
                      fshandi_diff_perc >= -30 & fshandi_diff_perc < -5 ~ "c-30%<s.div<-5%",
                      fshandi_diff_perc >= -5 & fshandi_diff_perc < 0 ~ "c-5%<s.div<0%",
                      fshandi_diff_perc == 0 ~ "s.div=0%",
                      fshandi_diff_perc > 0 & fshandi_diff_perc < 5 ~ "c0%<s.div<5%",
                      fshandi_diff_perc >= 5 & fshandi_diff_perc < 30 ~ "c5%<s.div<30%",
                      fshandi_diff_perc > 30 ~ "s.div>30%")
  ) %>% 
  mutate(fshandi_diff_perc_range = fct_relevel(fshandi_diff_perc_range, "s.div<-30%", "c-30%<s.div<-5%","c-5%<s.div<0%",
                                               "c0%<s.div<5%", "c5%<s.div<30%", 
                                            "s.div>30%")) 


## Shandi % difference - 5 classes ><5% ----
n_last_year <- 4
fshandi_1deg_diff_perc_tot_long_5classes <-  fshandi_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fshandi_diff_perc_45_2050, fshandi_diff_perc_45_2100, 
                fshandi_diff_perc_85_2050, fshandi_diff_perc_85_2100,
                fshandi_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fshandi_diff_perc"),     
    names_to = "rcp_year",
    names_prefix = "fshandi_diff_perc",
    values_to = "fshandi_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fshandi_diff_interval = cut(fshandi_diff_perc,
                                  breaks = c(-Inf, -5, -0.1, 0.1, 5, Inf),
                                  labels = c("< -5", "-5 to 0", "0", "0 to 5","> 5"),
                                  include.lowest = TRUE)) 


#plot

plot_fshandi_1deg_diff_perc_205045_5classes <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fshandi_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fshandi_1deg_diff_perc_205045_5classes


plot_fshandi_1deg_diff_perc_205085_5classes <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fshandi_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fshandi_1deg_diff_perc_205085_5classes


plot_fshandi_1deg_diff_perc_210045_5classes <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fshandi_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fshandi_1deg_diff_perc_210045_5classes


plot_fshandi_1deg_diff_perc_210085_5classes <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fshandi_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fshandi_1deg_diff_perc_210085_5classes


plot_fshandi_1deg_diff_perc_210026_5classes <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_26_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fshandi_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#31688E","#21918C","#35B779", "#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fshandi_1deg_diff_perc_210026_5classes


#species diversity + functional richness + functional redundancy

p_1deg_fshandi_5classes <- ggarrange ( plot_fshandi_1deg_diff_perc_210026_5classes,
                                    plot_fshandi_1deg_diff_perc_210045_5classes,
                                    plot_fshandi_1deg_diff_perc_210085_5classes,
                                    plot_fshandi_1deg_diff_perc_205045_5classes,
                                    plot_fshandi_1deg_diff_perc_205085_5classes,
                                    ncol = 3,
                                    nrow = 2,
                                    common.legend = TRUE,
                                    align = c("v")
)

# Display the composite plot
p_1deg_fshandi_5classes


## Boxplot shandi by latitude ----

#Now I unify the definition of coastal grid cells with the data from fshandi - I need just to copy this passage in the other files

fshandi_1deg_diff_perc_tot_long_5classes_2 <- fshandi_1deg_diff_perc_tot_long_5classes %>% 
  left_join(base_gam4, by="csquare_1deg") %>% 
  mutate (coasttype_general = "general") %>% 
  mutate (coasttype_ecoreg = coastal_ecoreg)


#HERE I TRY TO PLOT ALL THE BOXPLOT TOGETHER BY RCP

#Now I plot the boxplot
plot_bp_fshandi_1deg_diff_perc_2100_5classes <- fshandi_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (year == "2100") %>%
  ggplot(aes(x=fshandi_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~rcp, ncol = 3)+
  ylab("Latitude") +   
  xlab("Species diversity change (∆%)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1))

plot_bp_fshandi_1deg_diff_perc_2100_5classes


plot_bp_fshandi_1deg_diff_perc_2100_5classes_coastalecoreg <- fshandi_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (year == "2100") %>%
  filter (coastal_ecoreg =="coastal_ecoreg") %>% 
  ggplot(aes(x=fshandi_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~rcp, ncol = 3)+
  ylab("Latitude") +   
  xlab("Species diversity change (∆%)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1))


plot_bp_fshandi_1deg_diff_perc_2100_5classes_coastalecoreg

p_bp_1deg_fshandi_5classes_2100_tot_ecoreg <- ggarrange (plot_bp_fshandi_1deg_diff_perc_2100_5classes, 
                                                      plot_bp_fshandi_1deg_diff_perc_2100_5classes_coastalecoreg,
                                                      ncol = 1,
                                                      nrow = 2,
                                                      common.legend = TRUE,
                                                      align = c("v")
)

p_bp_1deg_fshandi_5classes_2100_tot_ecoreg


# Filter the data for the year and rcp of interest
fshandi_1deg_diff_perc_tot_long_5classes_3 <- fshandi_1deg_diff_perc_tot_long_5classes_2 %>%
  pivot_longer(
    cols = starts_with("coasttype"),     
    names_to = "bar_type",
    names_prefix = "coasttype",
    values_to = "subdivision",
    values_drop_na = TRUE
  ) 

# Create the combined boxplot
combined_boxplot_fshandi <- fshandi_1deg_diff_perc_tot_long_5classes_3 %>%
  filter (year == "2100") %>% 
  filter(subdivision=="coastal_ecoreg" | subdivision=="general") %>% 
  ggplot(aes(x = fshandi_diff_interval, y = centre_latitude)) +
  geom_boxplot(aes(fill = subdivision), position = position_dodge(width = 0.75)) +
  theme_minimal() +
  facet_wrap(~rcp, ncol = 3) +
  ylab("Latitude") +
  xlab("Species diversity change (∆%)") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("general" = "grey42", "coastal_ecoreg" = "grey80"))  

# Plot the combined boxplot
print(combined_boxplot_fshandi)


#### FRedundancy abs difference ----
fredund_202200 <- sp_div_3 %>% 
  filter(year== "s2022") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_2022=fredund) 

fredund_205045 <- sp_div_3 %>% 
  filter(scen== "s205045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_205045=fredund)%>% 
  dplyr::select(csquare_1deg, fredund_205045)

fredund_205085 <- sp_div_3 %>% 
  filter(scen== "s205085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_205085=fredund)%>% 
  dplyr::select(csquare_1deg, fredund_205085)  

fredund_210026 <- sp_div_3 %>% 
  filter(scen== "s210026") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_210026=fredund)%>% 
  dplyr::select(csquare_1deg, fredund_210026)

fredund_210045 <- sp_div_3 %>% 
  filter(scen== "s210045") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_210045=fredund)%>% 
  dplyr::select(csquare_1deg, fredund_210045)

fredund_210085 <- sp_div_3 %>% 
  filter(scen== "s210085") %>% 
  dplyr::select(csquare_1deg, centre_latitude, centre_longitude, rcp, year, fredund) %>%
  rename (fredund_210085=fredund)%>% 
  dplyr::select(csquare_1deg, fredund_210085)

#I unify the tables to create one unique dataframe from which I make the differences
fredund_1deg_diff_tot <- fredund_202200 %>% 
  full_join(fredund_205045, by="csquare_1deg") %>% 
  full_join(fredund_205085, by="csquare_1deg") %>% 
  full_join(fredund_210026, by="csquare_1deg") %>% 
  full_join(fredund_210045, by="csquare_1deg") %>% 
  full_join(fredund_210085, by="csquare_1deg") %>% 
  mutate (fredund_diff_45_2050 = fredund_205045 - fredund_2022) %>% 
  mutate (fredund_diff_45_2100 = fredund_210045 - fredund_2022) %>% 
  mutate (fredund_diff_85_2050 = fredund_205085 - fredund_2022) %>% 
  mutate (fredund_diff_85_2100 = fredund_210085 - fredund_2022) %>%
  mutate (fredund_diff_26_2100 = fredund_210026 - fredund_2022)


## FRedund % difference ----

fredund_1deg_diff_perc_tot <- fredund_1deg_diff_tot %>% 
  mutate (fredund_diff_perc_45_2050 = fredund_diff_45_2050/fredund_2022*100) %>% 
  mutate (fredund_diff_perc_45_2100 = fredund_diff_45_2100/fredund_2022*100) %>% 
  mutate (fredund_diff_perc_85_2050 = fredund_diff_85_2050/fredund_2022*100) %>% 
  mutate (fredund_diff_perc_85_2100 = fredund_diff_85_2100/fredund_2022*100) %>%
  mutate (fredund_diff_perc_26_2100 = fredund_diff_26_2100/fredund_2022*100)


#prepare file
n_last_year <- 4
fredund_1deg_diff_perc_tot_long <-  fredund_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fredund_diff_perc_45_2050, fredund_diff_perc_45_2100, 
                fredund_diff_perc_85_2050, fredund_diff_perc_85_2100,
                fredund_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fredund_diff_perc"),   
    names_to = "rcp_year",
    names_prefix = "fredund_diff_perc",
    values_to = "fredund_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fredund_diff_perc_range = 
           case_when (fredund_diff_perc < -30 ~ "f.redund<-30%",
                      fredund_diff_perc >= -30 & fredund_diff_perc < -5 ~ "c-30%<f.redund<-5%",
                      fredund_diff_perc >= -5 & fredund_diff_perc < 0 ~ "c-5%<f.redund<0%",
                      fredund_diff_perc == 0 ~ "f.redund=0%",
                      fredund_diff_perc > 0 & fredund_diff_perc < 5 ~ "c0%<f.redund<5%",
                      fredund_diff_perc >= 5 & fredund_diff_perc < 30 ~ "c5%<f.redund<30%",
                      fredund_diff_perc > 30 ~ "f.redund>30%")
  ) %>% 
  mutate(fredund_diff_perc_range = fct_relevel(fredund_diff_perc_range, "f.redund<-30%", "c-30%<f.redund<-5%","c-5%<f.redund<0%",
                                               "c0%<f.redund<5%", "c5%<f.redund<30%", 
                                               "f.redund>30%")) 



##FRedund % difference - 5 classes ><5% ----

n_last_year <- 4
fredund_1deg_diff_perc_tot_long_5classes <-  fredund_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fredund_diff_perc_45_2050, fredund_diff_perc_45_2100, 
                fredund_diff_perc_85_2050, fredund_diff_perc_85_2100,
                fredund_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fredund_diff_perc"),     
    names_to = "rcp_year",
    names_prefix = "fredund_diff_perc",
    values_to = "fredund_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fredund_diff_interval = cut(fredund_diff_perc,
                                     breaks = c(-Inf, -5, -0.1, 0.1, 5, Inf),
                                     labels = c("< -5", "-5 to 0", "0", "0 to 5","> 5"),
                                     include.lowest = TRUE)) 


#plot
plot_fredund_1deg_diff_perc_205045_5classes <- fredund_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fredund_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c( "#FDE725","#35B779", "#21918C","#31688E","#440154")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fredund_1deg_diff_perc_205045_5classes


plot_fredund_1deg_diff_perc_205085_5classes <- fredund_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2050") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fredund_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c( "#FDE725","#35B779", "#21918C","#31688E","#440154")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fredund_1deg_diff_perc_205085_5classes


plot_fredund_1deg_diff_perc_210045_5classes <- fredund_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_45_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fredund_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c( "#FDE725","#35B779", "#21918C","#31688E","#440154")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fredund_1deg_diff_perc_210045_5classes


plot_fredund_1deg_diff_perc_210085_5classes <- fredund_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_85_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fredund_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values =c( "#FDE725","#35B779", "#21918C","#31688E","#440154")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fredund_1deg_diff_perc_210085_5classes


plot_fredund_1deg_diff_perc_210026_5classes <- fredund_1deg_diff_perc_tot_long_5classes %>% 
  filter (rcp_year == "_26_2100") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = fredund_diff_interval)) +
  geom_raster() +
  scale_fill_manual(values = c( "#FDE725","#35B779", "#21918C","#31688E","#440154")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_fredund_1deg_diff_perc_210026_5classes


#species diversity + functional richness + functional redundancy

p_1deg_fredund_5classes <- ggarrange ( plot_fredund_1deg_diff_perc_210026_5classes,
                                       plot_fredund_1deg_diff_perc_210045_5classes,
                                       plot_fredund_1deg_diff_perc_210085_5classes,
                                       plot_fredund_1deg_diff_perc_205045_5classes,
                                       plot_fredund_1deg_diff_perc_205085_5classes,
                                       ncol = 3,
                                       nrow = 2,
                                       common.legend = TRUE,
                                       align = c("v")
)

# Display the composite plot
p_1deg_fredund_5classes





## Boxplot fredund by latitude ----

n_last_year <- 4
fredund_1deg_diff_perc_tot_long_5classes_invert <-  fredund_1deg_diff_perc_tot %>%
  dplyr::select(csquare_1deg, centre_longitude, centre_latitude,
                fredund_diff_perc_45_2050, fredund_diff_perc_45_2100, 
                fredund_diff_perc_85_2050, fredund_diff_perc_85_2100,
                fredund_diff_perc_26_2100) %>% 
  pivot_longer(
    cols = starts_with("fredund_diff_perc"),     
    names_to = "rcp_year",
    names_prefix = "fredund_diff_perc",
    values_to = "fredund_diff_perc",
    values_drop_na = TRUE
  ) %>% 
  mutate(fredund_diff_perc_invert = fredund_diff_perc * (-1)) %>% 
  mutate(rcp = substr(rcp_year, 2, 3)) %>% 
  mutate(year = substr(rcp_year, 5, 8)) %>% 
  mutate(fredund_diff_interval = cut(fredund_diff_perc_invert,  
                                     breaks = c(-Inf, -5, -0.1, 0.1, 5, Inf),
                                     labels = c("< -5", "-5 to 0", "0", "0 to 5","> 5"),
                                     include.lowest = TRUE)) 


#Now I unify the definition of coastal grid cells with the data from fredund 

fredund_1deg_diff_perc_tot_long_5classes_2 <- fredund_1deg_diff_perc_tot_long_5classes_invert %>% 
  left_join(base_gam4, by="csquare_1deg") %>% 
  mutate (coasttype_general = "general") %>% 
  mutate (coasttype_ecoreg = coastal_ecoreg)


#Now I plot the boxplot
plot_bp_fredund_1deg_diff_perc_2100_5classes <- fredund_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (year == "2100") %>%
  ggplot(aes(x=fredund_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~rcp, ncol = 3)+
  ylab("Latitude") +   
  xlab("F. redundancy change (∆%)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1))

plot_bp_fredund_1deg_diff_perc_2100_5classes


plot_bp_fredund_1deg_diff_perc_2100_5classes_coastalecoreg <- fredund_1deg_diff_perc_tot_long_5classes_2 %>% 
  filter (year == "2100") %>%
  filter (coastal_ecoreg =="coastal_ecoreg") %>% 
  ggplot(aes(x=fredund_diff_interval , y=centre_latitude)) +
  geom_boxplot()+
  theme_minimal()+
  facet_wrap(~rcp, ncol = 3)+
  ylab("Latitude") +   
  xlab("F. redundancy change (∆%)") +
  theme(text = element_text(size=10), axis.text.x = element_text(angle = 45, hjust=1))


plot_bp_fredund_1deg_diff_perc_2100_5classes_coastalecoreg

p_bp_1deg_fredund_5classes_2100_tot_ecoreg <- ggarrange (plot_bp_fredund_1deg_diff_perc_2100_5classes, 
                                                         plot_bp_fredund_1deg_diff_perc_2100_5classes_coastalecoreg,
                                                         ncol = 1,
                                                         nrow = 2,
                                                         common.legend = TRUE,
                                                         align = c("v")
)

p_bp_1deg_fredund_5classes_2100_tot_ecoreg


# Filter the data for the year and rcp of interest
fredund_1deg_diff_perc_tot_long_5classes_3 <- fredund_1deg_diff_perc_tot_long_5classes_2 %>%
  pivot_longer(
    cols = starts_with("coasttype"),     
    names_to = "bar_type",
    names_prefix = "coasttype",
    values_to = "subdivision",
    values_drop_na = TRUE
  ) 

# Create the combined boxplot
combined_boxplot_fredund <- fredund_1deg_diff_perc_tot_long_5classes_3 %>%
  filter (year == "2100") %>% 
  filter(subdivision=="coastal_ecoreg" | subdivision=="general") %>% 
  ggplot(aes(x = fredund_diff_interval, y = centre_latitude)) +
  geom_boxplot(aes(fill = subdivision), position = position_dodge(width = 0.75)) +
  theme_minimal() +
  facet_wrap(~rcp, ncol = 3) +
  ylab("Latitude") +
  xlab("F. redudancy change (∆%)") +
  theme(text = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("general" = "grey42", "coastal_ecoreg" = "grey80")) 

# Plot the combined boxplot
print(combined_boxplot_fredund)



p_bp_fdi_tot <- ggarrange (combined_boxplot_fshandi, 
                           combined_boxplot_fric,
                           combined_boxplot_fredund,
                           ncol = 1,
                           nrow = 3,
                           common.legend = TRUE,
                           align = c("v")
)

p_bp_fdi_tot






#COTANG fric /shandi ----

fric_1deg_diff_perc_tot_long_cotang <- fric_1deg_diff_perc_tot_long %>%
  na.omit() %>%
  mutate(c1deg = paste(csquare_1deg, rcp_year, sep = ""))

fshandi_1deg_diff_perc_tot_long_cotang <- fshandi_1deg_diff_perc_tot_long %>%
  na.omit() %>%
  mutate(c1deg = paste(csquare_1deg, rcp_year, sep = ""))

fredund_1deg_diff_perc_tot_long_cotang <- fredund_1deg_diff_perc_tot_long %>%
  na.omit() %>%
  mutate(c1deg = paste(csquare_1deg, rcp_year, sep = ""))

cotang_1deg <- fshandi_1deg_diff_perc_tot_long_cotang %>%
  dplyr::select(c1deg, fshandi_diff_perc, fshandi_diff_perc_range) %>% 
 full_join(fric_1deg_diff_perc_tot_long_cotang, by = "c1deg") %>% 
  full_join(fredund_1deg_diff_perc_tot_long_cotang, by = "c1deg")




#Fredund percentages by rcp and year ----

percent_fredund_byyearrcp <- cotang_1deg %>% 
  dplyr::select(rcp_year.y, fredund_diff_perc, fredund_diff_perc_range, c1deg, centre_longitude.x, centre_latitude.x)




#PLOT COTANG ----

###### Cotang with categories ----
#I start from the file for the cotang plot:
library(readr)
cotang_1deg_csv <- read_csv("test_127x15_GRID_221103_abund/base_corr.csv") %>% 
  mutate(cotang = fric_diff_perc/fshandi_diff_perc) %>% 
  mutate(cotang_range = 
           case_when (cotang < -30 ~ "f.ric/s.div.<-30",
                      cotang >= -30 & cotang < -5 ~ "-30<f.ric/s.div.<-5",
                      cotang >= -5 & cotang < -1 ~ "-5<f.ric/s.div.<-1",
                      cotang == -1 ~ "f.ric/s.div.= -1",
                      cotang > -1 & cotang < 0 ~ "-1<f.ric/s.div.<0",
                      cotang >= 0 & cotang < 1 ~ "0<f.ric/s.div.<1",
                      cotang == 1 ~ "f.ric/s.div.= 1",
                      cotang > -1 & cotang < 5 ~ "1<f.ric/s.div.<5",
                      cotang >= 5 & cotang < 30 ~ "5<f.ric/s.div.<30",
                      cotang > 30 ~ "f.ric/s.div.>30")
  ) %>% 
  mutate(cotang_range = fct_relevel(cotang_range, "f.ric/s.div.<-30", "-30<f.ric/s.div.<-5",
                                     "-5<f.ric/s.div.<-1", "-1<f.ric/s.div.<0", "0<f.ric/s.div.<1",
                                    "1<f.ric/s.div.<5", "5<f.ric/s.div.<30", "f.ric/s.div.>30" )) %>% 
  mutate (concordance = 
            case_when (fric_diff_perc > 0 & fshandi_diff_perc > 0 ~ "f.ric & sp.div increase",
                       fric_diff_perc < 0 & fshandi_diff_perc < 0 ~ "f.ric & sp.div decrease",
                       fric_diff_perc > 0 & fshandi_diff_perc < 0 ~ "f.ric increases & sp.div decreases",
                       fric_diff_perc < 0 & fshandi_diff_perc > 0 ~ "f.ric decreases & sp.div increases",
                        fric_diff_perc = 0 & fshandi_diff_perc < 0 ~"f.ric doesn't change & sp.ric decreases",
                        fric_diff_perc = 0 & fshandi_diff_perc > 0 ~ "f.ric doesn't change & sp.ric increases"
                        ))  %>% 
  mutate(concordance = if_else(is.na(concordance) & fshandi_diff_perc < 0, 
                               "f.ric doesn't change & sp.ric decreases",
                               if_else(is.na(concordance) & fshandi_diff_perc > 0,
                                       "f.ric doesn't change & sp.ric increases",
                                       if_else(is.na(concordance) & fshandi_diff_perc == 0,
                                               "f.ric & sp.ric don't change",
                                               concordance)))) %>% 
  #this part is for the combination of cotang and concordance:
  mutate (cotang_concord = 
            case_when (#both increase
                      fric_diff_perc > 0 & fshandi_diff_perc > 0 & cotang > 1 ~ "f.ric & s.div increase & |f.ric| > |s.div|",
                       fric_diff_perc > 0 & fshandi_diff_perc > 0 & cotang < 1 & cotang > 0 ~ "f.ric & s.div increase & |f.ric| < |s.div|",
                       #both decrease
                       fric_diff_perc < 0 & fshandi_diff_perc < 0 & cotang > 1 ~ "f.ric & s.div decrease & |f.ric| > |s.div|",
                       fric_diff_perc < 0 & fshandi_diff_perc < 0 & cotang < 1 & cotang > 0 ~ "f.ric & s.div decrease & |f.ric| < |s.div|",
                       #fric increases and shandi decreases
                       fric_diff_perc > 0 & fshandi_diff_perc < 0 & cotang < -1 ~ "f.ric increases & s.div decreases & |f.ric| > |s.div|",
                       fric_diff_perc > 0 & fshandi_diff_perc < 0 & cotang < 0 & cotang > -1 ~ "f.ric increases & s.div decreases & |f.ric| < |s.div|",
                       #fric decreases and shandi increases
                       fric_diff_perc < 0 & fshandi_diff_perc > 0 & cotang < -1 ~ "f.ric decreases & s.div increases & |f.ric| > |s.div|",
                       fric_diff_perc < 0 & fshandi_diff_perc > 0 & cotang < 0 & cotang > -1 ~ "f.ric decreases & s.div increases & |f.ric| < |s.div|", 
                       )) %>% 
  mutate(cotang_concord = if_else(is.na(cotang_concord) & fshandi_diff_perc < 0, 
                               "f.ric doesn't change & sp.ric decreases",
                               if_else(is.na(cotang_concord) & fshandi_diff_perc > 0,
                                       "f.ric doesn't change & sp.ric increases",
                                       if_else(is.na(cotang_concord) & fshandi_diff_perc == 0,
                                               "f.ric & sp.ric don't change",
                                               cotang_concord)))) %>% 
  mutate(cotang_concord = fct_relevel(cotang_concord, 
                                    "f.ric & s.div increase & |f.ric| > |s.div|", 
                                    "f.ric & s.div increase & |f.ric| < |s.div|", 
                                    "f.ric & s.div decrease & |f.ric| > |s.div|",
                                    "f.ric & s.div decrease & |f.ric| < |s.div|",
                                    "f.ric increases & s.div decreases & |f.ric| > |s.div|",
                                    "f.ric increases & s.div decreases & |f.ric| < |s.div|",
                                    "f.ric decreases & s.div increases & |f.ric| > |s.div|",
                                    "f.ric decreases & s.div increases & |f.ric| < |s.div|",
                                    "f.ric doesn't change & sp.ric decreases",
                                    "f.ric doesn't change & sp.ric increases"
                                      ))
  
test_contang <- cotang_1deg_csv %>% 
  filter(is.na(cotang_concord))



#####P cotang percent value ----
#percent change
plot_cotang_1deg_diff_perc_85 <- cotang_1deg_csv %>% 
  filter (rcp == "85") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_range)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39","#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_cotang_1deg_diff_perc_85
ggsave("16_cotang_1deg_diff_perc_85.tiff", plot_cotang_1deg_diff_perc_85 , width = 5, height = 5, dpi = 300)


#percent change
plot_cotang_1deg_diff_perc_45 <- cotang_1deg_csv %>% 
  filter (rcp == "45") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_range)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39","#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_cotang_1deg_diff_perc_45
ggsave("16_cotang_1deg_diff_perc_45.tiff", plot_cotang_1deg_diff_perc_45 , width = 5, height = 5, dpi = 300)

#percent change
plot_cotang_1deg_diff_perc_26 <- cotang_1deg_csv %>% 
  filter (rcp == "26") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_range)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39","#FDE725")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))


plot_cotang_1deg_diff_perc_26
ggsave("16_cotang_1deg_diff_perc_26.tiff", plot_cotang_1deg_diff_perc_26 , width = 5, height = 5, dpi = 300)

#####Plot concordance ----
#Concordance  rcp 85
plot_concordance_1deg_diff_perc_85 <- cotang_1deg_csv %>% 
  filter (rcp == "85") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = concordance)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#FDE725", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))

plot_concordance_1deg_diff_perc_85

#Concordance rcp 45
plot_concordance_1deg_diff_perc_45 <- cotang_1deg_csv %>% 
  filter (rcp == "45") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = concordance)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#FDE725", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))

plot_concordance_1deg_diff_perc_45

#Concordance rcp 26
plot_concordance_1deg_diff_perc_26 <- cotang_1deg_csv %>% 
  filter (rcp == "26") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = concordance)) +
  geom_raster() +
  scale_fill_manual(values = c("#440154", "#FDE725", "#46327E", "#365C8D","#277F8E", "#1FA187","#4AC16D", "#A0DA39", "#408AB3"
)) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))

plot_concordance_1deg_diff_perc_26


#####P cotang & concordance categories ----

library(devtools)
install_github("ericpante/marmap")
install.packages(mapdata)
install.packages("mapdata")
library(marmap) 
library(mapdata)

#Concordance nel caso di rcp 85
plot_cotang_concord_1deg_diff_perc_85 <- cotang_1deg_csv %>% 
  filter (rcp == "85") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_concord)) +
  geom_raster() +
  scale_fill_manual(values = c("#FDE725", "#6DCD59", "#440154", "#3B528B","#D55E00", "#F07F27","#408AB3", "#56B4E9", "#F025FD", "#F025FD")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75)) 
 
plot_cotang_concord_1deg_diff_perc_85

#Concordance rcp 45
plot_cotang_concord_1deg_diff_perc_45 <- cotang_1deg_csv %>% 
  filter (rcp == "45") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_concord)) +
  geom_raster() +
  scale_fill_manual(values = c("#FDE725", "#6DCD59", "#440154", "#3B528B","#D55E00", "#F07F27","#408AB3", "#56B4E9", "#F025FD", "#F025FD")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))

plot_cotang_concord_1deg_diff_perc_45


#Concordance rcp 26
plot_cotang_concord_1deg_diff_perc_26 <- cotang_1deg_csv %>% 
  filter (rcp == "26") %>% 
  ggplot(aes(centre_longitude, centre_latitude, fill = cotang_concord)) +
  geom_raster() +
  scale_fill_manual(values = c("#FDE725", "#6DCD59", "#440154", "#3B528B","#D55E00", "#F07F27","#408AB3", "#56B4E9", "#F025FD", "#F025FD")) +
  facet_wrap(~year, ncol = 1 ) +
  theme_minimal()+
  borders() + 
  coord_cartesian(xlim = c(-169.25, -70), ylim= c(10, 67.75))

plot_cotang_concord_1deg_diff_perc_26




#GAM ----
library(mgcv)
library(readr)
library(tidyverse)
library(dismo)

base_gam2 <- read_csv("test_127x15_GRID_221103_abund/base_corr.csv")  

n_last_yearrcp=6

base_gam <- base_gam2 %>% 
    mutate(yearrcp= substr(c1deg, nchar(c1deg) - n_last_yearrcp + 1, nchar(c1deg)))

 View(base_gam)
 
 
  base_gam2 <- base_gam %>% 
    mutate(centre_longitude2=centre_longitude %% 360) %>%  
    mutate(cotang_abs = fric/fshandi) %>% 
    mutate(fredund = (fsimp - fmpd))
  
  


#### Fric - GAM models ----
#0 - only temperature
fric_sst_mod <- mgcv::gam(fric ~ s(sst_1deg),
            data=base_gam2
                )

#1 - spatial model
fric_latlong_mod <- mgcv::gam(fric ~ s(centre_longitude2, centre_latitude, 
                                       k=100, bs="gp", m=c(3,1)),
                                  data=base_gam2
                                    )

#2 - spatiotemporal model   
fric_latlong_temp_mod <-
  mgcv::gam(fric ~ s(centre_longitude2,centre_latitude,
                        k = 100, bs = "gp", m = c(3, 1)
                          ) +
                yearrcp, 
            data = base_gam2)   

summary(fric_latlong_temp_mod )

#3 - temperature and spatial model
fric_sst_latlong_mod <- mgcv::gam(fric ~ s(sst_1deg) +
                                    s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)),
                                  data=base_gam2
                                    )

#4 - temperature and spatiotemporal model 
fric_sst_latlong_temp_mod <- mgcv::gam(fric ~ s(sst_1deg) +
                                    s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                   yearrcp,
                                  data=base_gam2
)

#5 - only temperature and depth
fric_sst_depth_mod <- mgcv::gam(fric ~ s(sst_1deg) + s(depth_mean_1d),
                          data=base_gam2
)

#6 - spatial model and depth
fric_latlong_depth_mod <- mgcv::gam(fric ~ s(centre_longitude2, centre_latitude, 
                                       k=100, bs="gp", m=c(3,1)) +
                                      s(depth_mean_1d),
                              data=base_gam2
)


#7 - temperature and spatial model and depth
fric_sst_latlong_depth_mod <- mgcv::gam(fric ~ s(sst_1deg) +
                                               s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                               s(depth_mean_1d),
                                  data=base_gam2
)

#8 - temperature and spatial model and depth
fric_sst_latlong_temp_depth_mod <- mgcv::gam(fric ~ s(sst_1deg) +
                                          s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                          s(depth_mean_1d) + yearrcp,
                                        data=base_gam2
)
summary(fric_sst_latlong_temp_depth_mod) 

#MODEL SELECTION
fric_sst_rsq <- summary(fric_sst_mod)$r.sq ## R-squared
fric_sst_AIC <- AIC(fric_sst_mod) ## AIC

fric_latlong_rsq <- summary(fric_latlong_mod)$r.sq ## R-squared
fric_latlong_AIC <- AIC(fric_latlong_mod) ## AIC

fric_sst_latlong_rsq <- summary(fric_sst_latlong_mod)$r.sq ## R-squared
fric_sst_latlong_AIC <- AIC(fric_sst_latlong_mod) ## AIC

fric_sst_depth_rsq <- summary(fric_sst_depth_mod)$r.sq ## R-squared
fric_sst_depth_AIC <- AIC(fric_sst_depth_mod) ## AIC

fric_latlong_depth_rsq <- summary(fric_latlong_depth_mod)$r.sq ## R-squared
fric_latlong_depth_AIC <- AIC(fric_latlong_depth_mod) ## AIC

fric_sst_latlong_depth_rsq <- summary(fric_sst_latlong_depth_mod)$r.sq ## R-squared
fric_sst_latlong_depth_AIC <- AIC(fric_sst_latlong_depth_mod) ## AIC

fric_latlong_temp_mod_rsq <- summary(fric_latlong_temp_mod)$r.sq ## R-squared
fric_latlong_temp_mod_AIC <- AIC(fric_latlong_temp_mod) ## AIC

fric_sst_latlong_temp_mod_rsq <- summary(fric_sst_latlong_temp_mod)$r.sq ## R-squared
fric_sst_latlong_temp_mod_AIC <- AIC(fric_sst_latlong_temp_mod) ## AIC

fric_sst_latlong_temp_depth_mod_rsq <- summary(fric_sst_latlong_temp_depth_mod)$r.sq ## R-squared
fric_sst_latlong_temp_depth_mod_AIC <- AIC(fric_sst_latlong_temp_depth_mod) ## AIC


model <- c("fric_sst", "fric_latlong", "fric_sst_latlong", "fric_sst_depth", "fric_latlong_depth",
      "fric_sst_latlong_depth", "fric_latlong_temp_mod", "fric_sst_latlong_temp_mod", "fric_sst_latlong_temp_depth_mod")
rsq <- c(fric_sst_rsq, fric_latlong_rsq, fric_sst_latlong_rsq, 
         fric_sst_depth_rsq, fric_latlong_depth_rsq, fric_sst_latlong_depth_rsq,
         fric_latlong_temp_mod_rsq, fric_sst_latlong_temp_mod_rsq, fric_sst_latlong_temp_depth_mod_rsq)
AIC <- c(fric_sst_AIC, fric_latlong_AIC, fric_sst_latlong_AIC, 
         fric_sst_depth_AIC, fric_latlong_depth_AIC, fric_sst_latlong_depth_AIC,
         fric_latlong_temp_mod_AIC, fric_sst_latlong_temp_mod_AIC, fric_sst_latlong_temp_depth_mod_AIC)

mod_fric_selection <- data.frame(model, rsq, AIC)
mod_fric_selection


#Plot smoothers


plot(fric_sst_latlong_temp_depth_mod, scale = 0, scheme = 2, rug=TRUE)

#RESIDUALS  
gam.check(fric_sst_latlong_temp_depth_mod)
summary(fric_sst_latlong_temp_depth_mod$residuals)
shapiro.test(fric_sst_latlong_temp_depth_mod$residuals)

#significance of model terms
anova(fric_sst_latlong_depth_mod)
anova(fric_sst_latlong_temp_depth_mod)

#Characterizing model accuracy

#Root mean squared error (RMSE)
# create a new data frame with the same predictor variables as base_gam
new_data_fric_f <- data.frame(sst_1deg = base_gam2$sst_1deg,
                       centre_longitude2 = base_gam2$centre_longitude2,
                       centre_latitude = base_gam2$centre_latitude,
                       depth_mean_1d = base_gam2$depth_mean_1d,
                       yearrcp = base_gam2$yearrcp )

# generate predictions using the new data frame
base_gam2$pred_fric_f <- predict(fric_sst_latlong_temp_depth_mod, newdata = new_data_fric_f, type = "response")

# Compute RMSE
summary(base_gam2)

#remove missing values
na_rows_pred_fric_f <- which(is.na(base_gam2$pred_fric_f))
na_rows_pred_fric_f

#I run the rmse 
base_gam_naomit_fric_f <- base_gam[!is.na(base_gam$pred_fric_f), ]
rmse_fric_f <- sqrt(mean((base_gam_naomit_fric_f$fric - base_gam_naomit_fric_f$pred_fric_f)^2))

rmse_fric_f 



#### Shandi - GAM models ----
#0 - only temperature
fshandi_sst_mod <- mgcv::gam(fshandi ~ s(sst_1deg),
                          data=base_gam2
)

#1 - spatial model
fshandi_latlong_mod <- mgcv::gam(fshandi ~ s(centre_longitude2, centre_latitude, 
                                       k=100, bs="gp", m=c(3,1)),
                              data=base_gam2
)

#2 - spatiotemporal model   
fshandi_latlong_temp_mod <-
  mgcv::gam(fshandi ~ s(centre_longitude2,centre_latitude,
                     k = 100, bs = "gp", m = c(3, 1)
  ) +
    yearrcp, 
  data = base_gam2)   

summary(fshandi_latlong_temp_mod )

#3 - temperature and spatial model
fshandi_sst_latlong_mod <- mgcv::gam(fshandi ~ s(sst_1deg) +
                                    s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)),
                                  data=base_gam2
)

#4 - temperature and spatiotemporal model 
fshandi_sst_latlong_temp_mod <- mgcv::gam(fshandi ~ s(sst_1deg) +
                                         s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                         yearrcp,
                                       data=base_gam2
)

#5 - only temperature and depth
fshandi_sst_depth_mod <- mgcv::gam(fshandi ~ s(sst_1deg) + s(depth_mean_1d),
                                data=base_gam2
)

#6 - spatial model and depth
fshandi_latlong_depth_mod <- mgcv::gam(fshandi ~ s(centre_longitude2, centre_latitude, 
                                             k=100, bs="gp", m=c(3,1)) +
                                      s(depth_mean_1d),
                                    data=base_gam2
)


#7 - temperature and spatial model and depth
fshandi_sst_latlong_depth_mod <- mgcv::gam(fshandi ~ s(sst_1deg) +
                                          s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                          s(depth_mean_1d),
                                        data=base_gam2
)

#8 - temperature and spatial model and depth
fshandi_sst_latlong_temp_depth_mod <- mgcv::gam(fshandi ~ s(sst_1deg) +
                                               s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                               s(depth_mean_1d) + yearrcp,
                                             data=base_gam2
)



#MODEL SELECTION
fshandi_sst_rsq <- summary(fshandi_sst_mod)$r.sq ## R-squared
fshandi_sst_AIC <- AIC(fshandi_sst_mod) ## AIC

fshandi_latlong_rsq <- summary(fshandi_latlong_mod)$r.sq ## R-squared
fshandi_latlong_AIC <- AIC(fshandi_latlong_mod) ## AIC

fshandi_sst_latlong_rsq <- summary(fshandi_sst_latlong_mod)$r.sq ## R-squared
fshandi_sst_latlong_AIC <- AIC(fshandi_sst_latlong_mod) ## AIC

fshandi_sst_depth_rsq <- summary(fshandi_sst_depth_mod)$r.sq ## R-squared
fshandi_sst_depth_AIC <- AIC(fshandi_sst_depth_mod) ## AIC

fshandi_latlong_depth_rsq <- summary(fshandi_latlong_depth_mod)$r.sq ## R-squared
fshandi_latlong_depth_AIC <- AIC(fshandi_latlong_depth_mod) ## AIC

fshandi_sst_latlong_depth_rsq <- summary(fshandi_sst_latlong_depth_mod)$r.sq ## R-squared
fshandi_sst_latlong_depth_AIC <- AIC(fshandi_sst_latlong_depth_mod) ## AIC

fshandi_latlong_temp_mod_rsq <- summary(fshandi_latlong_temp_mod)$r.sq ## R-squared
fshandi_latlong_temp_mod_AIC <- AIC(fshandi_latlong_temp_mod) ## AIC

fshandi_sst_latlong_temp_mod_rsq <- summary(fshandi_sst_latlong_temp_mod)$r.sq ## R-squared
fshandi_sst_latlong_temp_mod_AIC <- AIC(fshandi_sst_latlong_temp_mod) ## AIC

fshandi_sst_latlong_temp_depth_mod_rsq <- summary(fshandi_sst_latlong_temp_depth_mod)$r.sq ## R-squared
fshandi_sst_latlong_temp_depth_mod_AIC <- AIC(fshandi_sst_latlong_temp_depth_mod) ## AIC


model <- c("fshandi_sst", "fshandi_latlong", "fshandi_sst_latlong", "fshandi_sst_depth", "fshandi_latlong_depth",
           "fshandi_sst_latlong_depth", "fshandi_latlong_temp_mod", "fshandi_sst_latlong_temp_mod", "fshandi_sst_latlong_temp_depth_mod")
rsq <- c(fshandi_sst_rsq, fshandi_latlong_rsq, fshandi_sst_latlong_rsq, 
         fshandi_sst_depth_rsq, fshandi_latlong_depth_rsq, fshandi_sst_latlong_depth_rsq,
         fshandi_latlong_temp_mod_rsq, fshandi_sst_latlong_temp_mod_rsq, fshandi_sst_latlong_temp_depth_mod_rsq)
AIC <- c(fshandi_sst_AIC, fshandi_latlong_AIC, fshandi_sst_latlong_AIC, 
         fshandi_sst_depth_AIC, fshandi_latlong_depth_AIC, fshandi_sst_latlong_depth_AIC,
         fshandi_latlong_temp_mod_AIC, fshandi_sst_latlong_temp_mod_AIC, fshandi_sst_latlong_temp_depth_mod_AIC)

mod_fshandi_selection <- data.frame(model, rsq, AIC)
mod_fshandi_selection


#Plot smoothers
plot(fshandi_sst_latlong_temp_depth_mod, scale = 0, scheme = 2, rug=TRUE)

summary(fshandi_sst_latlong_temp_depth_mod)


#RESIDUALS 
gam.check(fshandi_sst_latlong_temp_depth_mod)
summary(fshandi_sst_latlong_temp_depth_mod$residuals)
shapiro.test(fshandi_sst_latlong_temp_depth_mod$residuals)


#significance of model terms
anova(fshandi_sst_latlong_temp_depth_mod)


#COMPUTE RMSE
# create a new data frame 
new_data_fshandi_f <- data.frame(sst_1deg = base_gam2$sst_1deg,
                              centre_longitude2 = base_gam2$centre_longitude2,
                              centre_latitude = base_gam2$centre_latitude,
                              depth_mean_1d = base_gam2$depth_mean_1d,
                              yearrcp = base_gam2$yearrcp )

# generate predictions using the new data frame
base_gam2$pred_fshandi_f <- predict(fshandi_sst_latlong_temp_depth_mod, newdata = new_data_fshandi_f, type = "response")

# Compute RMSE

summary(base_gam2)

na_rows_pred_fshandi_f <- which(is.na(base_gam2$pred_fshandi_f))
na_rows_pred_fshandi_f

base_gam_naomit_fshandi_f <- base_gam2[!is.na(base_gam2$pred_fshandi_f), ]
rmse_fshandi_f <- sqrt(mean((base_gam_naomit_fshandi_f$fshandi - base_gam_naomit_fshandi_f$pred_fshandi_f)^2))

rmse_fshandi_f 








#### F.redund - GAM models ----
#0 - only temperature
fredund_sst_mod <- mgcv::gam(fredund ~ s(sst_1deg),
                          data=base_gam2
)

#1 - spatial model
fredund_latlong_mod <- mgcv::gam(fredund ~ s(centre_longitude2, centre_latitude, 
                                       k=100, bs="gp", m=c(3,1)),
                              data=base_gam2
)

#2 - spatiotemporal model   
fredund_latlong_temp_mod <-
  mgcv::gam(fredund ~ s(centre_longitude2,centre_latitude,
                     k = 100, bs = "gp", m = c(3, 1)
  ) +
    yearrcp, 
  data = base_gam2)   

summary(fredund_latlong_temp_mod )

#3 - temperature and spatial model
fredund_sst_latlong_mod <- mgcv::gam(fredund ~ s(sst_1deg) +
                                    s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)),
                                  data=base_gam2
)

#4 - temperature and spatiotemporal model 
fredund_sst_latlong_temp_mod <- mgcv::gam(fredund ~ s(sst_1deg) +
                                         s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                         yearrcp,
                                       data=base_gam2
)

#5 - only temperature and depth
fredund_sst_depth_mod <- mgcv::gam(fredund ~ s(sst_1deg) + s(depth_mean_1d),
                                data=base_gam2
)

#6 - spatial model and depth
fredund_latlong_depth_mod <- mgcv::gam(fredund ~ s(centre_longitude2, centre_latitude, 
                                             k=100, bs="gp", m=c(3,1)) +
                                      s(depth_mean_1d),
                                    data=base_gam2
)


#7 - temperature and spatial model and depth
fredund_sst_latlong_depth_mod <- mgcv::gam(fredund ~ s(sst_1deg) +
                                          s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                          s(depth_mean_1d),
                                        data=base_gam2
)

#8 - temperature and spatial model and depth
fredund_sst_latlong_temp_depth_mod <- mgcv::gam(fredund ~ s(sst_1deg) +
                                               s(centre_longitude2, centre_latitude, k=100, bs="gp", m=c(3,1)) +
                                               s(depth_mean_1d) + yearrcp,
                                             data=base_gam2
)
summary(fredund_sst_latlong_temp_depth_mod) 


#MODEL SELECTION
fredund_sst_rsq <- summary(fredund_sst_mod)$r.sq ## R-squared
fredund_sst_AIC <- AIC(fredund_sst_mod) ## AIC

fredund_latlong_rsq <- summary(fredund_latlong_mod)$r.sq ## R-squared
fredund_latlong_AIC <- AIC(fredund_latlong_mod) ## AIC

fredund_sst_latlong_rsq <- summary(fredund_sst_latlong_mod)$r.sq ## R-squared
fredund_sst_latlong_AIC <- AIC(fredund_sst_latlong_mod) ## AIC

fredund_sst_depth_rsq <- summary(fredund_sst_depth_mod)$r.sq ## R-squared
fredund_sst_depth_AIC <- AIC(fredund_sst_depth_mod) ## AIC

fredund_latlong_depth_rsq <- summary(fredund_latlong_depth_mod)$r.sq ## R-squared
fredund_latlong_depth_AIC <- AIC(fredund_latlong_depth_mod) ## AIC

fredund_sst_latlong_depth_rsq <- summary(fredund_sst_latlong_depth_mod)$r.sq ## R-squared
fredund_sst_latlong_depth_AIC <- AIC(fredund_sst_latlong_depth_mod) ## AIC

fredund_latlong_temp_mod_rsq <- summary(fredund_latlong_temp_mod)$r.sq ## R-squared
fredund_latlong_temp_mod_AIC <- AIC(fredund_latlong_temp_mod) ## AIC

fredund_sst_latlong_temp_mod_rsq <- summary(fredund_sst_latlong_temp_mod)$r.sq ## R-squared
fredund_sst_latlong_temp_mod_AIC <- AIC(fredund_sst_latlong_temp_mod) ## AIC

fredund_sst_latlong_temp_depth_mod_rsq <- summary(fredund_sst_latlong_temp_depth_mod)$r.sq ## R-squared
fredund_sst_latlong_temp_depth_mod_AIC <- AIC(fredund_sst_latlong_temp_depth_mod) ## AIC


model <- c("fredund_sst", "fredund_latlong", "fredund_sst_latlong", "fredund_sst_depth", "fredund_latlong_depth",
           "fredund_sst_latlong_depth", "fredund_latlong_temp_mod", "fredund_sst_latlong_temp_mod", "fredund_sst_latlong_temp_depth_mod")
rsq <- c(fredund_sst_rsq, fredund_latlong_rsq, fredund_sst_latlong_rsq, 
         fredund_sst_depth_rsq, fredund_latlong_depth_rsq, fredund_sst_latlong_depth_rsq,
         fredund_latlong_temp_mod_rsq, fredund_sst_latlong_temp_mod_rsq, fredund_sst_latlong_temp_depth_mod_rsq)
AIC <- c(fredund_sst_AIC, fredund_latlong_AIC, fredund_sst_latlong_AIC, 
         fredund_sst_depth_AIC, fredund_latlong_depth_AIC, fredund_sst_latlong_depth_AIC,
         fredund_latlong_temp_mod_AIC, fredund_sst_latlong_temp_mod_AIC, fredund_sst_latlong_temp_depth_mod_AIC)

mod_fredund_selection <- data.frame(model, rsq, AIC)
mod_fredund_selection


#Plot smoothers
plot(fredund_sst_latlong_temp_depth_mod, scale = 0, scheme = 2, rug=TRUE)

#RESIDUALS 
gam.check(fredund_sst_latlong_temp_depth_mod)
summary(fredund_sst_latlong_temp_depth_mod$residuals)
shapiro.test(fredund_sst_latlong_temp_depth_mod$residuals)

#significance of model terms
anova(fredund_sst_latlong_depth_mod)
anova(fredund_sst_latlong_temp_depth_mod)

#Characterizing model accuracy

# create a new data frame with the same predictor variables as base_gam
new_data <- data.frame(sst_1deg = base_gam2$sst_1deg,
                       centre_longitude2 = base_gam2$centre_longitude2,
                       centre_latitude = base_gam2$centre_latitude,
                       depth_mean_1d = base_gam2$depth_mean_1d,
                       yearrcp = base_gam2$yearrcp)

# generate predictions using the new data frame
base_gam2$pred <- predict(fredund_sst_latlong_temp_depth_mod, newdata = new_data, type = "response")

# Compute RMSE

summary(base_gam2)

na_rows_pred <- which(is.na(base_gam2$pred))
na_rows_pred


base_gam_naomit <- base_gam2[!is.na(base_gam2$pred), ]
rmse <- sqrt(mean((base_gam_naomit$fredund - base_gam_naomit$pred)^2))

rmse 

#I compare with the other models:

base_gam$pred2 <- predict(fredund_sst_latlong_mod, newdata = new_data, type = "response")

# Compute RMSE
summary(base_gam)
na_rows_pred <- which(is.na(base_gam$pred))
na_rows_pred

base_gam <- na.omit(base_gam)
rmse <- sqrt(mean((base_gam$fredund - base_gam$pred)^2))

rmse 

#END




