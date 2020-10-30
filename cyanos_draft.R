
#transform data
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)


#calculate abundance for each Class
Dor_ps.ra.long_Cyanos <- Dor_ps.ra.long %>% 
  select(location, Season, Month,Year,OTU,Class, Family, Genus, Abundance)%>%
  filter(Class == "Cyanobacteriia") %>% 
  group_by(location,Year,Season,Month,Family, Genus) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  filter(Abund.total>0)

#reorder location
Dor_ps.ra.long_Cyanos$location <- factor(Dor_ps.ra.long_Cyanos$location,
                                         levels=c("Res.","V2.","D1."))


barplots_cyanos<- ggplot(Dor_ps.ra.long_Cyanos, aes(x = Month, y = Abund.total, fill = Family))+
  geom_col()+
  #ylim(0,100)+
  facet_grid(location~Year, space= "fixed")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom")

ggsave("./figures/Run1_barplots_cyanos.pdf", 
       plot = barplots_cyanos,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)







Dor_ps.prev.no.na <- Dor_ps.ra %>% 
  subset_samples(
    !is.na(MC_ug_L) & 
      !is.na(Fucoxanthin_mg_L) &
      !is.na(Dinoxanthin_mg_L) &
      !is.na(Lutein_mg_L) &
      !is.na(Zeaxanthin_mg_L) &
      !is.na(Chl_b_mg_L) &
      !is.na(NO3_NO2_N_L) &
      !is.na(TP_ug_L) &
      !is.na(Diatoxanthin_mg_L) &
      !is.na(b_caroten_mg_L) &
      !is.na(Ammonia_ug_L) &
      !is.na(Temp_degC)&
      !is.na(Chl_a_mg_L))

Dor_ps.prev.no.na <- prune_taxa(taxa_sums(Dor_ps.prev.no.na)>0,Dor_ps.prev.no.na)




#extrct metadata
env.par<- c("Temp_degC", "Chl_a_mg_L", "Chl_b_mg_L",
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L",
            "b_caroten_mg_L","MC_ug_L","Lutein_mg_L","Zeaxanthin_mg_L")


Res_metadata_long<- as(sample_data(Dor_ps.prev.no.na),"data.frame") %>% 
  select(location, Season, Year, Month, Temp_degC,Chl_a_mg_L,
         Ammonia_ug_L, NO3_NO2_N_L, TP_ug_L,MC_ug_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = "none",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014")))


Res_pigmets<- data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
  select(location, Season, Year, Month,#Chl_a_mg_L, Chl_b_mg_L,
         Diatoxanthin_mg_L,Dinoxanthin_mg_L,Fucoxanthin_mg_L,
         b_caroten_mg_L,Lutein_mg_L,Zeaxanthin_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = factor(variable),
         variable ="Pigments",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014")))


#calculate abundance for each Class
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev.no.na, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)

Cyanos_melt <- Dor_ps.ra.long %>% 
  select(location, Season, Month,Year,OTU,Class, Family, Abundance)%>%
  filter(Class == "Cyanobacteriia") %>% 
  group_by(location,Year,Season,Month,Family) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  select(location, Season, Year, Month,Family, Abund.total) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month","Family"), measure.vars = "Abund.total") %>% 
  mutate(Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014")),
         Family = factor(Family),
         variable= "Cyanobacteria")




test <- bind_rows(Res_metadata_long,Res_pigmets, Cyanos_melt) %>% 
  mutate(variable = factor(variable, levels = c("Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                                "Cyanobacteria","Pigments"),labels = c("Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                                                                       "Cyanobacteria","Pigments")),
         Family = factor(Family, levels = c(levels(Cyanos_melt$Family),levels(Res_pigmets$Family),"none"),
                         labels = c(levels(Cyanos_melt$Family),levels(Res_pigmets$Family),"none")))


ggplot()+
  geom_col(data = filter(test, variable %in% c("Cyanobacteria")), aes(x= Month, y = value, fill = Family, group = variable))+
  geom_line(data = filter(test, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(test, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(test, Family =="none"), aes(x= Month, y = value, group = variable),size = 1)+
  geom_point(data = filter(test, Family =="none"), aes(x= Month, y = value, group = variable), size = 2)+
  #geom_col(data = filter(test, variable =="Pigments"), aes(x= Month, y = value,  fill = Family, group = variable))+
  scale_fill_manual(values = tol24rainbow)+
  scale_colour_manual(values = tol24rainbow)+
  facet_grid(variable~location+Year, scales = "free_y")+
  theme_bw()



test<- as(sample_data(Dor_ps.prev.no.na),"data.frame")

