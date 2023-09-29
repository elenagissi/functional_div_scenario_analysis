# CORRELATION OF SPECIES TRAITS
# version manuscript Sept 13, 2023

# Workspace ----
rm( list=ls())

# set working directory 
setwd("~/MEDIX/WP2_Model/DATA")
library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Species by traits dataframe ----
traits_fishb_test <-
  read_delim(
    "int_fd/traits.csv",
    delim = ";",
    escape_double = FALSE,
    col_types = cols(
      length_max = col_number(),
      weight_max = col_number(),
      trophic_level = col_number(),
      depth_max = col_number()
    ),
    trim_ws = TRUE
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
traits$parental_care <- factor(traits$parental_care, order = TRUE, levels = c("none", "week_month", "month_year", ">year"))
traits$fecundity <- factor(traits$fecundity, order = TRUE, levels = c("<1", "1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100-1000", "1000-10000", ">10000" ))
traits$reproductive_mode <- factor(traits$reproductive_mode)
traits$body_shape <- factor(traits$body_shape)
traits$sexual_dimorphism <- factor(traits$sexual_dimorphism)
traits$sociality <- factor(traits$sociality, order = TRUE, levels = c("solitary", "small_groups", "medium_groups", "large_groups"))
traits$length_at_mat <- factor(traits$length_at_mat)
traits$ext_diet_spec <- factor(traits$ext_diet_spec)
traits$length_at_mat <-as.numeric(traits$length_at_mat)

# Trait type ----
traits_measure <- read_delim("int_fd/traits_measure.csv", 
                             delim = ";", escape_double = FALSE, col_types = cols(`Classes code` = col_integer()), 
                             trim_ws = TRUE) 


#Continuous variables - correlation ----
traits_continuous <- traits %>% 
  dplyr::select(length_max, weight_max, length_at_mat, trophic_level, depth_max) 

#correlation of continous variables towards Spearman coefficient:
install.packages("Hmisc")
library("Hmisc")

corr_continuos_base <- cor(traits_continuous)
corr_continuos <- rcorr(as.matrix(traits_continuous), type = "spearman")

# Extract the correlation coefficients
corr_continuos$r
# Extract p-values
corr_continuos$P

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
corr_continuos_matrix <- flattenCorrMatrix(corr_continuos$r, corr_continuos$P)

### Draw scatterplot ----
install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(traits_continuous, histogram=TRUE, method="spearman", pch=19)


#Categorical variables vs continuous variables - correlation ----

# Correlation through the Kruskal-wallis-test in r
install.packages("rstatix")
library(rstatix)
library(ggpubr)  
library(broom)

categorical_vars <- c('sexual_dimorphism',
                      'sociality',
                      'reproductive_strategy',
                      'reproductive_mode',
                      'parental_investment',
                      'parental_care',
                      'terrestriality',
                      'adult_mobility',
                      'zones',
                     'body_shape',
                      'diet_main_component',
                      'ext_diet_spec',
                     'fecundity')
numerical_vars <- c('length_max',
                    'weight_max',
                    'length_at_mat',
                    'trophic_level',
                    'depth_max')

#Kruskal test between categorical and numerical variables
kret = list() 
i <- 1 
for (c in categorical_vars) {
  for (n in numerical_vars) { 
    f <- as.formula(paste(n, '~', c)) 
    kret[[deparse(f)]] <- kruskal.test(f, data=traits) 
    i <- i + 1 
    } 
}


# Effect size 
keff = list() 
i <- 1 
for (c in categorical_vars) {
  for (n in numerical_vars) { 
    f <- as.formula(paste(n, '~', c))
    keff[[deparse(f)]] <- kruskal_effsize(f, data=traits) 
    i <- i + 1 
  } 
}


corr_KW <- as.data.frame(do.call (rbind, kret)) 
corr_KW$data.name <- as.character(corr_KW$data.name) 

corr_KW_2 <- corr_KW %>% 
  mutate(trait1= as.character(purrr::map(strsplit(data.name, split = " by "), 1)))%>% 
  mutate(trait2= as.character(purrr::map(strsplit(data.name, split = " by "), 2)))
corr_KW_2 <- as.matrix(corr_KW_2) 


corr_KWeff <- as.data.frame(do.call (rbind, keff)) 
corr_KWeff$data.name <- as.character(corr_KWeff$data.name)

corr_KWeff_2 <- as.matrix(corr_KWeff) 



#Categorical variables among them ----
traits_categ <- traits %>%
  dplyr::select(
    sexual_dimorphism,
    sociality,
    reproductive_strategy,
    reproductive_mode,
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


#The following package computes correlation together
install.packages('GoodmanKruskal')
library('GoodmanKruskal')

GKmatrix1 <- GKtauDataframe(traits_categ)  
plot(GKmatrix1)

#Contingency tables are reported in the supplementary methods of the manuscript.
# END