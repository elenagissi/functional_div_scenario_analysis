### Gissi et al. - Traits compilation
#version Sept. 13, 2023

# Workspace ----
rm( list=ls())
setwd("~/MEDIX")

install.packages("tidyverse")
library(tidyverse)
install.packages("dplyr")
library(dplyr)
install.packages("tidyr")
library(tidyr)


#Traits from databases in literature ----

library(readr)
trinidade_full_select <- read_delim("trinidade_full_select.csv", 
                                                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
albouy_select <- read_delim("albouy_select.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE)

waechter_select <- read_csv("waechter_select.csv")


butt_select_transpose2 <- read_csv("butt_select.csv")

colnames(butt_select_transpose2)
names(butt_select_transpose2)[names(butt_select_transpose2) == "spp_gp"] <- "spp_list_gissi"

colnames(waechter_select)
names(waechter_select) [names(waechter_select) == "SPP"] <- "spp_M_Trat_list_gissi"

#Join tables by gissi_spp
traits_join1 <- full_join(gissi_spp, butt_select_transpose2, by = "spp_list_gissi")
traits_join2 <- full_join(traits_join1, trinidade_full_select, by = "spp_list_gissi")
colnames(traits_join2)
traits_join3 <-
  full_join(traits_join2, waechter_select, by = "spp_M_Trat_list_gissi")
traits_full <- full_join(traits_join3, albouy_select, by = "spp_M_Trat_list_gissi")
colnames(traits_full) 

#clean
traits_full_colnames<- colnames(traits_full)
traits_full_colnames

traits_select <- traits_full %>%
  select(
    COD_ID,
    spp_list_gissi,
    spp_M_list_gissi,
    spp_M_Trat_list_gissi,
    adult_body_mass_body_size,
    al_adult_weight_max,
    Body_mass,
    adult_mobility,
    age_to_1st_reproduction_generation_time,
    tri_age_at_maturity,
    are_there_sub_populations_reproduct,
    can_the_sex_ratio_be_altered_by_temperature,
    depth_min_max,
    tri_depth,
    al_foraging_water_depth,
    Depth_category,
    extreme_diet_specialization,
    fecundity,
    lifetime_reproductive_opportunities,
    al_inter_litter,
    max_age,
    parental_investment,
    tri_parental_investment,
    ph,
    post_birth_hatching_parental_dependence,
    al_weaning,
    reproductive_strategy,
    sub_pop_dependence_breeding,
    sub_pop_dependence_foraging,
    terrestrial_and_marine_life_stages,
    al_terrestriality,
    thermal_sensitivity_to_heat_spikes_heat_waves,
    thermal_tolerance_range,
    tri_thermal_tolerance,
    tri_trophic_level,
    tri_body_shape,
    Body_shape,
    tri_water_column,
    al_dimorphism,
    al_social_group_mean,
    al_main_diet,
    Maximum_size,
    Trophic_group,
    Caudal_fin
  )




# Traits from Fishbase ----

install.packages("devtools")
library(devtools)
remotes::install_github("ropensci/rfishbase")
library(rfishbase)
library(dplyr)

gissi_spp <- read_delim("gissi_spp.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)

#running load_taxa()
fishbase <- load_taxa(server="fishbase")

#create a vector with species list
spp <- gissi_spp[,'spp_M_list_gissi']
#change name to the vector in spp
names(spp)[names(spp) == "spp_M_list_gissi"] <- "spp"

#validate names of my species list in fishbase 
sppv <- validate_names(spp$spp, server="fishbase")
sppv 



## Species table ----
sppv_species_fishb <- species(sppv, server="fishbase") %>% 
  filter(!is.na(Species))
sppv_species_fishb

sppv_species_fishb_select <- sppv_species_fishb %>% 
    select(SpecCode, Species, Genus, SpeciesRefNo, FBname, BodyShapeI, 
           DemersPelag, AnaCat, DepthRangeShallow, DepthRangeDeep, 
           DepthRangeComShallow, DepthRangeComDeep, LongevityWild, Vulnerability, 
           Length, LTypeMaxM, LengthFemale, LTypeMaxF, CommonLength, LTypeComM, 
           CommonLengthF, LTypeComF, Weight, WeightFemale, Importance, 
           PriceCateg, PriceReliability, PD50)
  
  

 #I rename the traits with my code
  traits_g1_species_fishb <- sppv_species_fishb_select %>% 
    rename( body_shape = BodyShapeI, water_column_position = DemersPelag, adult_mobility = AnaCat, depth_min = DepthRangeShallow, depth_max = DepthRangeDeep, 
         depth_com_min = DepthRangeComShallow, depth_com_max = DepthRangeComDeep, longevity = LongevityWild, vulnerability = Vulnerability, 
         length_max = Length, length_com = CommonLength, weight_max = Weight, uniqueness_pd50 = PD50)

## Ecology ----
ecology_fishb <- ecology(sppv, server="fishbase") %>% 
  filter(!is.na(Species))

troph_fishb <- ecology_fishb %>% 
  select(SpecCode, Species, FoodTroph, FoodSeTroph, DietTroph, DietSeTroph, FeedingType)

 

watercolumn_fishb <- ecology_fishb %>% 
  select(SpecCode, Species, Neritic, SupraLittoralZone, Saltmarshes, 
         LittoralZone, TidePools, Intertidal, SubLittoral, Caves, Oceanic, 
         Epipelagic, Mesopelagic, Bathypelagic, Abyssopelagic, Hadopelagic,
         Estuaries, Mangroves, MarshesSwamps, CaveAnchialine, Stream, Lakes, 
         Cave, Cave2)



traits_g2_troph_fishb <- troph_fishb %>% 
  rename(trophic_level_food = FoodTroph, trophic_level_diet = DietTroph)
  

## Reproduction ----
reproduct_fishb <- reproduction(sppv, server="fishbase") %>% 

reproduct_fishb_select <- reproduct_fishb %>% 
  select(SpecCode, Species, ReproducRefNo, ReproMode, Fertilization, RepGuild1,
         RepGuild2, ParentalCare, ParentalCareQ, AddInfos) %>% 
  mutate(parental_investment = paste(RepGuild1, RepGuild2, sep = "_")) %>% 
  rename(reproductive_strategy = ReproMode, parental_care_sex = ParentalCare) 


## Maturity ----
maturit_fishb <- maturity(sppv, server="fishbase") %>% 

maturit_fishb_select <- maturit_fishb %>% 
  select (SpecCode, Species, Sex, Lm, LengthMatMin, LengthMatMin2, tm, 
          AgeMatMin, AgeMatMin2, Comment)

popgrowth_fishb <- popgrowth(sppv, server="fishbase") %>% 
  filter(!is.na(Species)) 


## Length at maturity ----
library(dplyr)
maturit_fishb_lm_unsexed <- maturit_fishb %>%
  group_by(Species) %>%
  summarise(LengthMatMin_ave = min(LengthMatMin, na.rm = TRUE),
            LengthMatMin2_ave = max(LengthMatMin2, na.rm = TRUE),
            Lm_ave = mean(Lm, na.rm = TRUE),
            AgeMatMin_ave = min(AgeMatMin, na.rm = TRUE),
            AgeMatMin2_ave = max(AgeMatMin2, na.rm = TRUE),
            tm_ave = mean(tm, na.rm = TRUE)) %>% 
  rename(lm_ave_unsex = Lm_ave, length_mat_min_unsex = LengthMatMin_ave, 
         length_mat_max_unsex = LengthMatMin2_ave, age_mat_min_unsex = AgeMatMin_ave,
         age_mat_max_unsex = AgeMatMin2_ave, tm_ave_unsex = tm_ave)



## Swimming mode ----
swimming_fishbase <- swimming(sppv, server="fishbase") %>% 
  filter(!is.na(Species))   
str(swimming_fishbase)

swim_fishb_select <- swimming_fishbase %>% 
  select (SpecCode, Species, AdultType, AdultMode) %>% 
  rename (swimming_mode = AdultMode)


 
## Fecundity ----
fecund_fishb <- fecundity(sppv, server="fishbase") %>% 
  filter(!is.na(Species))   
str(fecund_fishb)

fecund_fishb_select <- fecund_fishb %>% 
select(SpecCode, Species, FecundityMin, FecundityMax, FecundityMean, SpawningCycles, AddInfos)

fecund_fishb_select_spp <- fecund_fishb_select %>%
  group_by(Species) %>%
  dplyr::summarize(fecundity_min = min(FecundityMin, na.rm = TRUE),
                   fecundity_max = max(FecundityMax, na.rm = TRUE),
                   fecundity_mean = mean(FecundityMean, na.rm = TRUE),
                   spawning_cycles_mean = mean(SpawningCycles, na.rm = TRUE))



## Diet traits ---- 
diet_fishb <- diet(sppv, server="fishbase") %>% 
  filter(!is.na(Species))

diet_items_fishb <- diet_items(server="fishbase") 

diet_items_fishb_max <- diet_items_fishb %>% 
  group_by(DietCode) %>% 
  filter(DietPercent== max(DietPercent)) %>% 
  select(DietCode, FoodI, FoodII, FoodIII, DietPercent, ItemName)

diet_fishb_s <- diet_fishb %>% 
  select(SpecCode, Species, DietCode, StockCode, SampleStage, Troph, seTroph, 
         SizeMin, SizeMax, FishLength) 

 diet_fishb_spp_max <- left_join(diet_fishb_s, diet_items_fishb_max, by="DietCode")
 
 diet_fishb_spp_max_synthesis <-  diet_fishb_spp_max %>%
   group_by(Species, FoodI) %>%
   summarise(FoodII = paste(FoodII, collapse = " - ")) %>% 
   group_by(Species) %>% 
   summarise(FoodI=)
 
 
 diet_fishb_spp_max_synthesis2 <-  diet_fishb_spp_max %>%
   group_by(Species) %>%
   summarise(FoodI = paste (FoodI, collapse=" - "), FoodII = paste (FoodII, collapse=" - "))
 
library(stringr)
             
 diet_fishb_spp_max_synthesis3 <-  diet_fishb_spp_max_synthesis2 %>% 
   mutate(FoodII_new = str_extract(FoodII, "fish"),
          FoodI_new = str_extract(FoodI, "zooplankton"))
             
 
 ##Join from Fishbase ----
  traits_fishb <- full_join(traits_g1_species_fishb, traits_g2_troph_fishb, by = "Species") %>% 
   full_join(reproduct_fishb_select, by = "Species") %>% 
   full_join(maturit_fishb_lm_unsexed, by = "Species") %>% 
   full_join(maturit_fishb_lm_sexed_transpose, by = "Species") %>% 
   full_join(swim_fishb_select, by = "Species") %>% 
   full_join(fecund_fishb_select_spp, by = "Species") %>% 
   full_join(diet_fishbase_cleaned, by = "Species")  


# Traits from Sealifebase ----

install.packages("devtools")
library(devtools)
remotes::install_github("ropensci/rfishbase")
library(rfishbase)
library(dplyr)
library(tidyr)
options(FISHBASE_API = "http://fishbase.ropensci.org/sealifebase")

setwd("~/MEDIX/WP2_Model/DATA")
library(readr)
gissi_spp <- read_delim("gissi_spp.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

#running load_taxa() 
sealifebase <- load_taxa(server="sealifebase")
str(sealifebase)

sealifebase_phylum <- sealifebase %>% 
  group_by(Phylum) %>% 
  count()

sealifebase_class <- sealifebase %>% 
  group_by(Class) %>% 
  count()

#species list
spp <- gissi_spp[,'spp_M_list_gissi']
names(spp)[names(spp) == "spp_M_list_gissi"] <- "spp"
sppv <- validate_names(spp$spp, server="sealifebase")
sppv 

## Species table
sppv_species_sealb <- species(sppv, server="sealifebase") %>% 
  filter(!is.na(Species))
sppv_species_sealb

sppv_species_sealb_select <- sppv_species_sealb %>%  TRAITS
  select(SpecCode, Species, Genus, SpeciesRefNo, FBname, BodyShapeI, 
         DemersPelag, AnaCat, DepthRangeShallow, DepthRangeDeep, 
         DepthRangeComShallow, DepthRangeComDeep, LongevityWild, Vulnerability, 
         Length, LTypeMaxM, LengthFemale, LTypeMaxF, CommonLength, LTypeComM, 
         CommonLengthF, LTypeComF, Weight, WeightFemale, Importance, 
         PriceCateg, PriceReliability, Comments)

#I rename the traits with my code
traits_g1_species_sealb <- sppv_species_sealb_select %>% 
  rename( body_shape = BodyShapeI, water_column_position = DemersPelag, adult_mobility = AnaCat, depth_min = DepthRangeShallow, depth_max = DepthRangeDeep, 
          depth_com_min = DepthRangeComShallow, depth_com_max = DepthRangeComDeep, longevity = LongevityWild, vulnerability = Vulnerability, 
          length_max = Length, length_com = CommonLength, weight_max = Weight)

## Ecology ----
ecology_sealb <- ecology(sppv, server="sealifebase") %>% 
  filter(!is.na(Species))

troph_sealb <- ecology_sealb %>% 
  select(SpecCode, Species, FoodTroph, FoodSeTroph, DietTroph, DietSeTroph, FeedingType)

 

watercolumn_sealb <- ecology_sealb %>% 
  select(SpecCode, Species, Neritic, SupraLittoralZone, saltmarshes, 
         LittoralZone, TidePools, Intertidal, SubLittoral, Oceanic, 
         Epipelagic, mesopelagic, bathypelagic, abyssopelagic, hadopelagic, 
         Estuaries, Mangroves, MarshesSwamps, CaveAnchialine, Stream, Lakes, 
         Cave, Cave2)

#I rename the traits with my code
traits_g2_troph_sealb <- troph_sealb %>% 
  rename(trophic_level_food = FoodTroph, trophic_level_diet = DietTroph)


## Reproduction ----
reproduct_sealb <- reproduction(sppv, server="sealifebase") %>% 
  filter(!is.na(Species))

reproduct_sealb_select <- reproduct_sealb %>% 
  select(SpecCode, Species, ReproducRefNo, ReproMode, Fertilization, RepGuild1,
         RepGuild2, AddInfos) %>% 
  mutate(parental_investment = paste(RepGuild1, RepGuild2, sep = "_")) %>% 
  rename(reproductive_strategy = ReproMode) 



## Maturity ----
maturit_sealb <- maturity(sppv, server="sealifebase") %>% 
  filter(!is.na(Species)) 

maturit_sealb_select <- maturit_sealb %>% 
  select (SpecCode, Species, Sex, Lm, LengthMatMin, LengthMatMin2, tm, 
          AgeMatMin, AgeMatMin2, Comment)

popgrowth_sealb <- popgrowth(sppv, server="sealifebase") %>% 
  filter(!is.na(Species)) 


## Length at maturity ----

maturit_sealb_lm_unsexed <- maturit_sealb %>%
  group_by(Species) %>%
  summarize(LengthMatMin_ave = min(LengthMatMin, na.rm = TRUE),
            LengthMatMin2_ave = max(LengthMatMin2, na.rm = TRUE),
            Lm_ave = mean(Lm, na.rm = TRUE),
            AgeMatMin_ave = min(AgeMatMin, na.rm = TRUE),
            AgeMatMin2_ave = max(AgeMatMin2, na.rm = TRUE),
            tm_ave = mean(tm, na.rm = TRUE)) %>% 
  rename(lm_ave_unsex = Lm_ave, length_mat_min_unsex = LengthMatMin_ave, 
         length_mat_max_unsex = LengthMatMin2_ave, age_mat_min_unsex = AgeMatMin_ave,
         age_mat_max_unsex = AgeMatMin2_ave, tm_ave_unsex = tm_ave) 



## Fecundity ----
fecund_sealb <- fecundity(sppv, server="sealifebase") %>% 
  filter(!is.na(Species))   
str(fecund_sealb)

fecund_sealb_select <- fecund_sealb %>% 
  select(SpecCode, Species, FecundityMin, FecundityMax, FecundityMean, SpawningCycles, AddInfos)

fecund_sealb_select_spp <- fecund_sealb_select %>%
  group_by(Species) %>%
  summarize(fecundity_min = min(FecundityMin, na.rm = TRUE),
            fecundity_max = max(FecundityMax, na.rm = TRUE),
            fecundity_mean = mean(FecundityMean, na.rm = TRUE),
            spawning_cycles_mean = mean(SpawningCycles, na.rm = TRUE))



## Diet traits ---- 
diet_sealb <- diet(sppv, server="sealifebase") %>% 
  filter(!is.na(Species))

diet_items_sealb <- diet_items(server="sealifebase") 

diet_items_sealb_max <- diet_items_sealb %>% 
  group_by(DietCode) %>% 
  filter(DietPercent== max(DietPercent)) %>% 
  select(DietCode, FoodI, FoodII, FoodIII, DietPercent, ItemName)

diet_sealb_s <- diet_sealb %>% 
  select(SpecCode, Species, DietCode, StockCode, SampleStage, Troph, seTroph, 
         SizeMin, SizeMax, FishLength) 

diet_sealb_spp_max <- left_join(diet_sealb_s, diet_items_sealb_max, by="DietCode")

diet_sealb_spp_max_synthesis <-  diet_sealb_spp_max %>%
  group_by(Species, FoodI) %>%
  summarise(FoodII = paste(FoodII, collapse = " - ")) %>% 
  group_by(Species) %>% 
  summarise(FoodI)


diet_sealb_spp_max_synthesis2 <-  diet_sealb_spp_max %>%
  group_by(Species) %>%
  summarise(FoodI = paste (FoodI, collapse=" - "), FoodII = paste (FoodII, collapse=" - "))

library(stringr)

diet_sealb_spp_max_synthesis3 <-  diet_sealb_spp_max_synthesis2 %>% 
  mutate(FoodII_new = str_extract(FoodII, "fish"),
         FoodI_new = str_extract(FoodI, "zooplankton"))



## Join from Sealifebase ----
#Partial combination of all the traits from sealifebase 
traits_sealb_clean2 <- full_join(traits_g1_species_sealb, traits_g2_troph_sealb, by = "Species") %>% 
  full_join(reproduct_sealb_select, by = "Species") %>% 
  full_join(maturit_sealb_lm_unsexed, by = "Species") %>% 
  full_join(fecund_sealb_select_spp, by = "Species")  %>% 
  full_join(diet_fishbase_cleaned, by = "Species")  

#END
  