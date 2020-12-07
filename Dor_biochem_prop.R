#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(tidyr); packageVersion("tidyr")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Extract metadata from Phyloseq
#####################################
#look only into samples of the fish culturing period (i.e., parameters from more than 1 pool)
n.samples<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  group_by(Year,Month) %>% 
  summarise(n = n())

Dor_biochem <- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  left_join(n.samples) %>% 
  filter(n >1) %>% 
  mutate(Redfield = as.numeric(as.numeric(NO3_NO2_N_L)/as.numeric(TP_ug_L)))

Dor_biochem.summary<-Dor_biochem %>% 
  melt(id = "location", measure.vars =c("Chl_a_mg_L","NO3_NO2_N_L", "TP_ug_L", "Redfield")) %>% 
  group_by(location, variable) %>% 
  summarise(mean = mean(as.numeric(value), na.rm = T), se = se(as.numeric(value), na.rm = T),
            n= length(value[!is.na(value)]))


Dor_fish_sum<- Dor_biochem %>%   melt(id = c("location","Year"), measure.vars =c("Fish_biomass_g_pond","Food_Kg_pond")) %>% 
  group_by(location, Year, variable) %>% 
  summarise(sum = sum(as.numeric(value), na.rm = T), n= length(value[!is.na(value)]))
