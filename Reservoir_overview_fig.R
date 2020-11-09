#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")
library(ggpubr); packageVersion("ggpubr")

#load colors and functions
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Seasonal dynamic in the reservoir
#####################################
#generate plot for physicochemical parameters
Res_parameters_long<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Season, Year, Month, Temp_degC,
         #Food_Kg_pond, Fish_biomass_g_pond, O2_mg_L, pH, Ammonia_ug_L, # at the moment excluded from the plot
         NO3_NO2_N_L, TP_ug_L, MC_ug_L) %>% 
  filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")))

# plot
Res_parameters.p<- ggplot(data = Res_parameters_long, aes(x= Month, y = value, group = variable))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  scale_fill_manual(values = tol24rainbow)+
  scale_colour_manual(values = tol24rainbow)+
  facet_grid(variable~Year, scales = "free_y")+
  theme_bw()+
  theme(legend.position="bottom")


#generate plot for pigments
Res_pigmets_long<- data.frame(sample_data(Dor_ps.prev)) %>% 
  select(location, Season, Year, Month,
         Diatoxanthin_mg_L,Dinoxanthin_mg_L,Fucoxanthin_mg_L,
         b_caroten_mg_L,Lutein_mg_L,Zeaxanthin_mg_L) %>% 
  filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Type ="Pigments",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")))

Res_chl_long <- data.frame(sample_data(Dor_ps.prev)) %>% 
  select(location, Season, Year, Month,Chl_a_mg_L, Chl_b_mg_L) %>% 
  filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Type ="Chlorophyll",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")))

Res_pig_chl <- bind_rows(Res_pigmets_long,Res_chl_long)

# plot
Res_pig_chl.p<- ggplot(data = Res_pig_chl, aes(x= Month, y = value, colour = variable, group = variable))+
  geom_line(size = 1)+
  geom_point(size = 2)+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = tol24rainbow)+
  facet_grid(Type~Year, scales = "free_y")+
  theme_bw()+
  theme(legend.position="none")


#Plot barplots of communities
#calculate proportions
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)
Dor_ps.ra.long$Abundance <- Dor_ps.ra.long$Abundance*100

#fix unclassified lineages 
Dor_ps.ra.long$Class <- as.character(Dor_ps.ra.long$Class)
Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Class)] <- paste(Dor_ps.ra.long$Phylum[is.na(Dor_ps.ra.long$Class)],"uc", sep = "_")

Dor_ps.ra.long$Order <- as.character(Dor_ps.ra.long$Order)
Dor_ps.ra.long$Order[is.na(Dor_ps.ra.long$Order)] <- paste(Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Order)],"uc", sep = "_")

Dor_ps.ra.long$Species <- as.character(Dor_ps.ra.long$Species)
Dor_ps.ra.long$Species[is.na(Dor_ps.ra.long$Species)] <- paste(Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Species)],"uc", sep = "_")

#calculate abundance for each Class
Dor_ps.ra.long.agg <- Dor_ps.ra.long %>% 
  select(location,Month,Year,OTU,Class,Abundance)%>%
  group_by(location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

#remove below 1% ra
taxa_classes <- unique(Dor_ps.ra.long.agg$Class)
Dor_ps.ra.long.agg$Class[Dor_ps.ra.long.agg$Abund.total<2] <- "Other taxa"
Dor_ps.ra.long.agg$Class <- factor(Dor_ps.ra.long.agg$Class,
                                   levels=c(taxa_classes,"Other taxa"))
Dor_ps.ra.long.agg$Class<- droplevels(Dor_ps.ra.long.agg$Class)

#reorder location
Dor_ps.ra.long.agg$location <- factor(Dor_ps.ra.long.agg$location,
                                      levels=c("Res.","V2.","D1."))

#extract reservoir data
Res_class.long <- Dor_ps.ra.long.agg %>%  filter(location == "Res.")


#Plot 
Res_class_barplot <- ggplot(Res_class.long, aes(x = Month, y = Abund.total, fill = Class)) + 
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  facet_grid(location~Year, space= "fixed") +
  theme_bw()+
  theme(legend.position="bottom")


#####################################
#Generate merged plot for the reservoir
#####################################
ggarrange(Res_parameters.p, Res_pig_chl.p, Res_class_barplot, heights = c(2,1,2),
          ncol = 1, nrow = 3, align = "v")







Dor_ps.class.agg <- Dor_ps.ra.long.agg %>% 
  select(location, Year, Month,Class, Abund.total) %>% 
  melt(id=c("location", "Year", "Month","Class"), measure.vars = "Abund.total") %>% 
  mutate(Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         Family = factor(Class),
         variable= "Bacteria")


#merge all data together
merged_df<- bind_rows(Res_metadata_long,Res_pigmets, Dor_ps.class.agg, Res_chl) %>% 
  mutate(variable = factor(variable, levels = c("Bacteria","Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                                "Pigments","Chl","Food_Kg_pond", "Fish_biomass_g_pond", "O2_mg_L", "pH"),
                           labels = c("Bacteria","Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                      "Pigments","Chl","Food_Kg_pond", "Fish_biomass_g_pond", "O2_mg_L", "pH")),
         Family = factor(Family, levels = c(levels(Dor_ps.class.agg$Family),levels(Res_pigmets$Family),levels(Res_chl$Family),"none"),
                         labels = c(levels(Dor_ps.class.agg$Family),levels(Res_pigmets$Family),levels(Res_chl$Family),"none"))) %>% 
  filter(location %in% c("Res.","V2.","D1."))


merged_df.p<- ggplot()+
  geom_col(data = filter(merged_df, variable %in% c("Bacteria")), aes(x= Month, y = value, fill = Family, group = variable))+
  geom_line(data = filter(merged_df, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(merged_df, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(merged_df, variable %in% c("Chl")), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(merged_df, variable %in% c("Chl")), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(merged_df, Family =="none"), aes(x= Month, y = value, group = variable),size = 1)+
  geom_point(data = filter(merged_df, Family =="none"), aes(x= Month, y = value, group = variable), size = 2)+
  #geom_col(data = filter(test, variable =="Pigments"), aes(x= Month, y = value,  fill = Family, group = variable))+
  scale_fill_manual(values = tol24rainbow)+
  scale_colour_manual(values = tol24rainbow)+
  facet_grid(variable~location+Year, scales = "free_y")+
  theme_bw()







Res_bac.p<- ggplot()+
  geom_col(data = filter(merged_df, variable %in% c("Bacteria")), aes(x= Month, y = value, fill = Family, group = variable))+
  scale_fill_manual(values = class_col)+
  facet_grid(variable~location+Year, scales = "free_y")+
  theme_bw()


  


ggarrange(Res_par.p, Res_bac.p, heights = c(2,1),
            ncol = 1, nrow = 2, align = "v")


