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


Dor_fish_sum<- Dor_biochem %>%   #melt(id = c("location","Year"), measure.vars =c("Fish_biomass_g_pond","Food_Kg_pond")) %>% 
  group_by(location, Year) %>% 
  summarise(Food.sum = sum(as.numeric(Food_Kg_pond), na.rm = T), Biomass = max(as.numeric(Fish_biomass_g_pond), na.rm = T), n= length(Food_Kg_pond[!is.na(Food_Kg_pond)]))


#####################################
#model food and fish biomass
#####################################
Fish_metadata<- read.csv2("tables/Fish_biomass_data.csv")


#build linear model
fit1 <- lm(Total.biomass..kg. ~ Day, data = subset(Fish_metadata, location == "D1."))


ggplot(Fish_metadata, aes(x = Day, y = Total.biomass..kg., group = location))+
  geom_point()+
  facet_wrap(location~., scales = "free_y")+
  stat_smooth(method = "lm", col = "blue") +
  labs(title = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 2),
                     "Intercept =",signif(fit1$coef[[1]],2),
                     " Slope =",signif(fit1$coef[[2]], 2),
                     " P =",signif(summary(fit1)$coef[2,4], 2)))







