#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")

#load colors and functions
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Subset samples that have metadata
####################################
#remove NAs in metadata
sample_data(Dor_ps.prev)[sample_data(Dor_ps.prev) == "#N/A"]<- NA
#identify ASVs that are present only in the substet of samples
Dor_ps.prev.sub <- Dor_ps.prev %>% 
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
#extract the list of ASVs that are present in these samples
ASVs_list<- taxa_names(prune_taxa(taxa_sums(Dor_ps.prev.sub)>0,Dor_ps.prev.sub))

#transform to geometric mean and remove 
Dor_ps.prev_gm<- phyloseq_gm_mean_trans(Dor_ps.prev)
Dor_ps.prev.no.na<-Dor_ps.prev_gm %>% 
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

#remove ASVs that were not observed in these samples
Dor_ps.prev.no.na<-prune_taxa(ASVs_list,Dor_ps.prev.no.na)

#####################################
#Check correlations between environmental parameters
####################################
env.par<- c("Temp_degC", "Chl_a_mg_L", "Chl_b_mg_L",
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L",
            "b_caroten_mg_L","MC_ug_L","Lutein_mg_L","Zeaxanthin_mg_L","Food_Kg_pond")

#scale parameters
metadata.scaled <- data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate_at(all_of(env.par),as.numeric) %>%
  mutate_if(is.numeric, scale_par)


#correlation between the env parameters
envpar_corr <- metadata.scaled %>% 
  select(all_of(env.par))%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr.pvalues<- envpar_corr %>% cor_get_pval()

#plot and save
png(file="figures/par_corr.png",
    width=1024, height=1024)
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)
dev.off()


#subset phys. parameters
phys_par.scaled<-metadata.scaled %>% 
  select(Temp_degC,  Food_Kg_pond, 
         Ammonia_ug_L, NO3_NO2_N_L,TP_ug_L)

#subset pigments and mycrocysteins
pig_par.scaled<-metadata.scaled %>% 
  select(Dinoxanthin_mg_L,Chl_a_mg_L, Chl_b_mg_L,Fucoxanthin_mg_L,
         b_caroten_mg_L,MC_ug_L,Lutein_mg_L,
         Diatoxanthin_mg_L,Zeaxanthin_mg_L)

#####################################
#RDA analysis of physical parameters
#####################################
#RDA analysis
Dor_ps.prev.rda.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., data = phys_par.scaled) # model including all variables 

Dor_ps.prev.rda.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ 1,
                         data = phys_par.scaled) # model containing only species matrix and intercept
Dor_ps.prev.rda.sel.os <- ordistep(Dor_ps.prev.rda.0, scope = formula(Dor_ps.prev.rda.all), direction = 'both',
                                   permutations = how(nperm = 999), steps = 100) #stepwise selection

#generate an RDA plot 
Dor_ps.rda.scores <- vegan::scores(Dor_ps.prev.rda.sel.os,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.sites <- data.frame(Dor_ps.rda.scores$sites)
Dor_ps.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Season <- as.character(sample_data(Dor_ps.prev.no.na)$Season )
Dor_ps.rda.sites <- Dor_ps.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate(Season = factor(Season, levels= c("Winter","Spring",
                                           "Summer","Autumn")))

#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot*5
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.sel.os$CCA$eig))

#Plot 
Dor_ps.rda.plot <- ggplot() +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, shape = Year), 
             fill = "black", size = 5) +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Year), 
             size = 4) +
  geom_text(data = Dor_ps.rda.sites,aes(x = RDA1, y = RDA2,label = location), 
            nudge_y= -0.3,size=3)+
  scale_colour_manual(values = c("Winter"="darkblue",
                                 "Spring"="lightblue",
                                 "Summer"="orange",
                                 "Autumn"="darkred")) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1.1),
            aes(x, y, label = rownames(Dor_ps.rda.arrows)),color="black",alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")


ggsave("./figures/RDA_step_phys.png", 
       plot = Dor_ps.rda.plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#RDA analysis of pigments
#####################################
#RDA analysis
Dor_ps.prev.rda.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., data = pig_par.scaled) # model including all variables 

Dor_ps.prev.rda.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ 1,
                         data = pig_par.scaled) # model containing only species matrix and intercept
Dor_ps.prev.rda.sel.os <- ordistep(Dor_ps.prev.rda.0, scope = formula(Dor_ps.prev.rda.all), direction = 'both',
                                   permutations = how(nperm = 999), steps = 100) #stepwise selection

#generate an RDA plot 
Dor_ps.rda.scores <- vegan::scores(Dor_ps.prev.rda.sel.os,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.sites <- data.frame(Dor_ps.rda.scores$sites)
Dor_ps.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Season <- as.character(sample_data(Dor_ps.prev.no.na)$Season )
Dor_ps.rda.sites <- Dor_ps.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate(Season = factor(Season, levels= c("Winter","Spring",
                                           "Summer","Autumn")))

#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot*5
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.sel.os$CCA$eig))

#Plot 
Dor_ps.rda.plot <- ggplot() +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, shape = Year), 
             fill = "black", size = 5) +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, colour = Season, shape = Year), 
             size = 4) +
  geom_text(data = Dor_ps.rda.sites,aes(x = RDA1, y = RDA2,label = location), 
            nudge_y= -0.3,size=3)+
  scale_colour_manual(values = c("Winter"="darkblue",
                                 "Spring"="lightblue",
                                 "Summer"="orange",
                                 "Autumn"="darkred")) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1.1),
            aes(x, y, label = rownames(Dor_ps.rda.arrows)),color="black",alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")


ggsave("./figures/RDA_step_pigments.png", 
       plot = Dor_ps.rda.plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



#####################################
#RDA analysis of Cyanobacteria
#####################################
#extract the data from RDA and consolidated with taxonomy
Dor_ps.rda.species <- data.frame(Dor_ps.rda.scores$species)
Dor_ps.rda.species$OTU <- as.character(rownames(Dor_ps.rda.species))

#Get the taxonomy of ASVs with at least 1% of seqeunces 
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev.sub, function(x) x / sum(x))
Dor_ps.ra.long <- psmelt(Dor_ps.ra)
Abund_ASVs <- Dor_ps.ra.long %>% filter(Abundance >0.03) %>% 
  select(OTU,Class, Order, Family, Genus) %>% 
  unique()

#merge taxonomy with the RDA data and filter only the cynos
Dor_ps.rda.species<- left_join(Dor_ps.rda.species[,c("RDA1","RDA2","OTU")],Abund_ASVs, by = "OTU") %>% # %>% na.omit()
                        filter(Class== "Cyanobacteriia")

#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.sel.os$CCA$eig))

#plot
Dor_ps.rda.species.plot <- ggplot() +
  geom_point(data = Dor_ps.rda.species, aes(x = RDA1, y = RDA2, fill = Family), 
             shape =21, size = 4) +
  geom_text(data = Dor_ps.rda.species,aes(x = RDA1, y = RDA2,label = str_sub(Genus, 0,4) ),
            nudge_y= -0.03,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  scale_fill_manual(values = tol21rainbow) +
  #scale_x_reverse()+ 
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1),
            aes(x, y, label = rownames(Dor_ps.rda.arrows)),color="black",alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("./figures/RDA_cyanos.png", 
       plot = Dor_ps.rda.species.plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)






#####################################
#Seasonal dynamic
#####################################
Res_metadata_long<- as(sample_data(Dor_ps.prev.no.na),"data.frame") %>% 
  select(location, Season, Year, Month, Temp_degC,Food_Kg_pond, Fish_biomass_g_pond, O2_mg_L, pH,
         Ammonia_ug_L, NO3_NO2_N_L, TP_ug_L, MC_ug_L) %>% 
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

Res_chl <- data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
  select(location, Season, Year, Month,Chl_a_mg_L, Chl_b_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = factor(variable),
         variable ="Chl",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014")))

#calculate abundance of cyanos
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

#merge all data together
merged_df<- bind_rows(Res_metadata_long,Res_pigmets, Cyanos_melt, Res_chl) %>% 
  mutate(variable = factor(variable, levels = c("Cyanobacteria","Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                                "Pigments","Chl","Food_Kg_pond", "Fish_biomass_g_pond", "O2_mg_L", "pH"),
                           labels = c("Cyanobacteria","Temp_degC","NO3_NO2_N_L","TP_ug_L","MC_ug_L","Ammonia_ug_L",
                                      "Pigments","Chl","Food_Kg_pond", "Fish_biomass_g_pond", "O2_mg_L", "pH")),
         Family = factor(Family, levels = c(levels(Cyanos_melt$Family),levels(Res_pigmets$Family),levels(Res_chl$Family),"none"),
                         labels = c(levels(Cyanos_melt$Family),levels(Res_pigmets$Family),levels(Res_chl$Family),"none"))) %>% 
  filter(location %in% c("Res.","V2."))


merged_df.p<- ggplot()+
  geom_col(data = filter(merged_df, variable %in% c("Cyanobacteria")), aes(x= Month, y = value, fill = Family, group = variable))+
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

ggsave("./figures/cyanos_overview.png", 
       plot = merged_df.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)
