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
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#transform to geometric mean and remove 
Dor_ps.prev_gm<- phyloseq_gm_mean_trans(Dor_ps.prev_run1)
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

#####################################
#Check correlations between environmental parameters
####################################
env.par<- c("Temp_degC", "Chl_a_mg_L", "Chl_b_mg_L", "N:P",
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L",
            "b_caroten_mg_L","MC_ug_L","Lutein_mg_L","Zeaxanthin_mg_L")

#scale parameters
metadata.scaled <- data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate("N:P"=(as.numeric(NO3_NO2_N_L)/as.numeric(TP_ug_L))) %>% 
  mutate_at(all_of(env.par),as.numeric) %>%
  mutate_if(is.numeric, scale_par)


#correlation between the env parameters
envpar_corr <- metadata.scaled %>% 
  select(all_of(env.par))%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr.pvalues<- envpar_corr %>% cor_get_pval()

#plot and save
pdf(file="figures/par_corr.pdf")
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)
dev.off()


#subset phys. parameters and mycrocysteins
phys_par.scaled<-metadata.scaled %>% 
  select(Temp_degC, "N:P",
         Ammonia_ug_L,MC_ug_L, 
         NO3_NO2_N_L,
         TP_ug_L)

#subset pigments 
pig_par.scaled<-metadata.scaled %>% 
  select(Chl_a_mg_L, Chl_b_mg_L,
         Fucoxanthin_mg_L,
         b_caroten_mg_L,#Lutein_mg_L,Dinoxanthin_mg_L,
         Diatoxanthin_mg_L,Zeaxanthin_mg_L)

#####################################
#RDA analysis of communities
#####################################
#stepwise RDA analysis
Dor_ps.prev.rda.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., data = phys_par.scaled) # model including all variables 

#Dor_ps.prev.rda.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ 1,
#                         data = phys_par.scaled) # model containing only species matrix and intercept
#Dor_ps.prev.rda.sel.os <- ordistep(Dor_ps.prev.rda.0, scope = formula(Dor_ps.prev.rda.all), direction = 'both',
#                                   permutations = how(nperm = 999), steps = 100) #stepwise selection

R2.phys<-RsquareAdj(Dor_ps.prev.rda.all)
P2.phys<-anova(Dor_ps.prev.rda.all)

#generate an RDA plot 
Dor_ps.rda.scores <- vegan::scores(Dor_ps.prev.rda.all,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.sites <- data.frame(Dor_ps.rda.scores$sites)
Dor_ps.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Mic.Season <- as.character(sample_data(Dor_ps.prev.no.na)$Mic.Season )
Dor_ps.rda.sites <- Dor_ps.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate(Season = factor(Mic.Season, levels= c("Wet","Dry")))

#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot*5
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * summary(Dor_ps.prev.rda.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.plot <- ggplot() +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor_ps.rda.sites, aes(x = RDA1, y = RDA2, colour = Mic.Season, shape = location), 
             size = 3) +
  geom_text(data = Dor_ps.rda.sites,aes(x = RDA1, y = RDA2,label = paste(Month,gsub("20","",Year))), 
          nudge_y= -0.8,size=5)+
  scale_colour_manual(values = c("Wet"="darkblue",
                                 "Dry"="orange")) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50")+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1.1),
            aes(x, y, label = substr(rownames(Dor_ps.rda.arrows), 1, 1)),color="gray50", size=4)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


#RDA analysis of pigments
Dor_ps.prev.rda.pig.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., data = pig_par.scaled) # model including all variables 

#Dor_ps.prev.rda.pig.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ 1,
#                             data = pig_par.scaled) # model containing only species matrix and intercept
#Dor_ps.prev.rda.pig.sel.os <- ordistep(Dor_ps.prev.rda.pig.0, scope = formula(Dor_ps.prev.rda.pig.all), direction = 'both',
#                                       permutations = how(nperm = 999), steps = 100) #stepwise selection

#summary of RDA
summary(Dor_ps.prev.rda.pig.all, display = NULL)

R2.pig<-RsquareAdj(Dor_ps.prev.rda.pig.all)
P2.pig<-anova(Dor_ps.prev.rda.pig.all)

#generate an RDA plot 
Dor_ps.rda.pig.scores <- vegan::scores(Dor_ps.prev.rda.pig.all,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.pig.sites <- data.frame(Dor_ps.rda.pig.scores$sites)
Dor_ps.rda.pig.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.pig.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Mic.Season <- as.character(sample_data(Dor_ps.prev.no.na)$Mic.Season)
Dor_ps.rda.pig.sites <- Dor_ps.rda.pig.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na), by = "Sample.ID") %>% 
  mutate(Season = factor(Mic.Season, levels= c("Wet","Dry")))

#Draw biplots
Dor_ps.rda.pig.arrows<- Dor_ps.rda.pig.scores$biplot*5
colnames(Dor_ps.rda.pig.arrows)<-c("x","y")
Dor_ps.rda.pig.arrows <- as.data.frame(Dor_ps.rda.pig.arrows)
Dor_ps.rda.pig.evals <- 100 * summary(Dor_ps.prev.rda.pig.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.pig.plot <- ggplot() +
  geom_point(data = Dor_ps.rda.pig.sites, aes(x = RDA1, y = RDA2, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor_ps.rda.pig.sites, aes(x = RDA1, y = RDA2, colour = Mic.Season, shape = location), 
             size = 3) +
  geom_text(data = Dor_ps.rda.pig.sites,aes(x = RDA1, y = RDA2,label = paste(Month,gsub("20","",Year))), 
            nudge_y= -0.8,size=5)+
  scale_colour_manual(values = c("Wet"="darkblue",
                                 "Dry"="orange")) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.pig.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.pig.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.pig.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50")+
  geom_text(data=as.data.frame(Dor_ps.rda.pig.arrows*1.1),
            aes(x, y, label = substr(rownames(Dor_ps.rda.pig.arrows), 1, 1)),color="gray50", size=4)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")



ggarrange(Dor_ps.rda.plot, Dor_ps.rda.pig.plot, widths = c(1,1),
          ncol = 2, nrow = 1, align = "hv", legend = "bottom", common.legend = TRUE)

ggsave("./figures/Env_par.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#RDA analysis of ASVs
#####################################
#import list of enriched ASVs
enriched_ASV<- read.csv("tables/enriched_ASVs.csv")

#generate an RDA plot of species
Dor_ps.rda.species <- data.frame(Dor_ps.rda.scores$species)
Dor_ps.rda.species$ASV <- as.character(rownames(Dor_ps.rda.species))
Dor_tax <- data.frame(as(tax_table(Dor_ps.prev.no.na), "matrix"), ASV = rownames(tax_table(Dor_ps.prev.no.na)))
Dor_ps.rda.species <- Dor_ps.rda.species %>%
  left_join(Dor_tax, by= "ASV") 

#calculate sequence proportions of each class and reduce the colouring of the plot only to the abundant classes
Dor_ps.prev.ra <- transform_sample_counts(Dor_ps.prev, function(x) x / sum(x))
Dor_ps.prev.long <- psmelt(Dor_ps.prev.ra)
Dor_ps.ra.long.agg <- Dor_ps.prev.long %>% 
  select(location,Month,Year,OTU,Class,Abundance)%>%
  group_by(location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

taxa_classes <- unique(Dor_ps.ra.long.agg$Class)
Dor_ps.ra.long.agg$Class[Dor_ps.ra.long.agg$Abund.total<0.02] <- "Other taxa"
Dor_ps.ra.long.agg$Class <- factor(Dor_ps.ra.long.agg$Class,
                                   levels=c(taxa_classes,"Other taxa"))
Dor_ps.ra.long.agg$Class<- droplevels(Dor_ps.ra.long.agg$Class)
Dor_ps.rda.species.sub <- Dor_ps.rda.species %>% mutate(Class =  ifelse(Class %in% Dor_ps.ra.long.agg$Class, Class, "Other taxa")) %>% 
                          mutate(Class = factor(Class,levels=c(taxa_classes,"Other taxa")),
                                 Alpha = case_when(ASV %in% enriched_ASV$ASV ~ 1, 
                                                   TRUE ~ 0.5))
                                                                           
#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * summary(Dor_ps.prev.rda.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.species.plot <- ggplot() +
  geom_point(data = subset(Dor_ps.rda.species.sub,Alpha == 0.5), aes(x = RDA1, y = RDA2), alpha = 0.2,
             fill = "gray", size = 1) +
  geom_point(data = subset(Dor_ps.rda.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2), alpha = 1,
             fill = "black", size = 5) +
  geom_point(data = subset(Dor_ps.rda.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2, colour = Class), alpha = 1,
             size = 4) +
  #geom_text(data = subset(Dor_ps.rda.species.sub, Alpha == 1),aes(x = RDA1, y = RDA2,label = Genus), 
  #          nudge_y= -0.01,size=3)+
  scale_colour_manual(values = class_col) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1.1),
            aes(x, y, label = substr(rownames(Dor_ps.rda.arrows), 1, 1)),color="gray50",alpha=0.5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#generate an RDA plot of species
Dor_ps.rda.pig.species <- data.frame(Dor_ps.rda.pig.scores$species)
Dor_ps.rda.pig.species$ASV <- as.character(rownames(Dor_ps.rda.pig.species))
Dor_tax <- data.frame(as(tax_table(Dor_ps.prev.no.na), "matrix"), ASV = rownames(tax_table(Dor_ps.prev.no.na)))
Dor_ps.rda.pig.species <- Dor_ps.rda.pig.species %>%
  left_join(Dor_tax, by= "ASV") 

#calculate sequence proportions of each class and reduce the colouring of the plot only to the abundant classes
Dor_ps.rda.pig.species.sub <- Dor_ps.rda.pig.species %>% mutate(Class =  ifelse(Class %in% Dor_ps.ra.long.agg$Class, Class, "Other taxa")) %>% 
  mutate(Class = factor(Class,levels=c(taxa_classes,"Other taxa")),
         Alpha = case_when(ASV %in% enriched_ASV$ASV ~ 1, 
                           TRUE ~ 0.5))

#Draw biplots
Dor_ps.rda.pig.arrows<- Dor_ps.rda.pig.scores$biplot
colnames(Dor_ps.rda.pig.arrows)<-c("x","y")
Dor_ps.rda.pig.arrows <- as.data.frame(Dor_ps.rda.pig.arrows)
Dor_ps.rda.pig.evals <- 100 * summary(Dor_ps.prev.rda.pig.all)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.pig.plot <- ggplot() +
  geom_point(data = subset(Dor_ps.rda.pig.species.sub,Alpha == 0.5), aes(x = RDA1, y = RDA2), alpha = 0.2,
             fill = "gray", size = 1) +
  geom_point(data = subset(Dor_ps.rda.pig.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2), alpha = 1,
             fill = "black", size = 5) +
  geom_point(data = subset(Dor_ps.rda.pig.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2, colour = Class), alpha = 1,
             size = 4) +
  #geom_text(data = subset(Dor_ps.rda.pig.species.sub, Alpha == 1),aes(x = RDA1, y = RDA2,label = Genus), 
   #         nudge_y= -0.01,size=3)+
  scale_colour_manual(values = class_col) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.pig.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.pig.evals[2], 2))) +
  geom_segment(data=Dor_ps.rda.pig.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50")+
  geom_text(data=as.data.frame(Dor_ps.rda.pig.arrows*1.1),
            aes(x, y, label = substr(rownames(Dor_ps.rda.pig.arrows), 1, 1)),color="gray50", size=4)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


ggarrange(Dor_ps.rda.species.plot, Dor_ps.rda.pig.plot, labels ="AUTO",
          ncol = 2, nrow = 1, align = "hv", legend = "bottom", common.legend = TRUE)

ggsave("./figures/RDA_species.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#Effect of aquaculture
#####################################
sample_data(Dor_ps.prev_gm)$Fish_biomass_Kg_pond[is.na(sample_data(Dor_ps.prev_gm)$Fish_biomass_Kg_pond)]<- 0
sample_data(Dor_ps.prev_gm)$Food_Kg_pond[is.na(sample_data(Dor_ps.prev_gm)$Food_Kg_pond)]<- 0


Dor_ps.prev.fishponds<-Dor_ps.prev_gm %>% 
  subset_samples(
    !is.na(Fish_biomass_Kg_pond) & 
      !is.na(Food_Kg_pond))


#subset and scale parameters
#scale parameters
Fish_par.scaled <- data.frame(sample_data(Dor_ps.prev.fishponds)) %>% 
  select(Food_Kg_pond,Fish_biomass_Kg_pond) %>% 
  mutate_all(as.numeric) %>%
  mutate_if(is.numeric, scale_par)

#RDA analysis
Dor_ps.fish.rda.0 <- rda(t(otu_table(Dor_ps.prev.fishponds)) ~ ., data = Fish_par.scaled)

R2.fish<-RsquareAdj(Dor_ps.fish.rda.0)
P2.fish<-anova(Dor_ps.fish.rda.0)

#generate an RDA plot 
Dor_ps.fish.rda.scores <- vegan::scores(Dor_ps.fish.rda.0,display=c("sp","wa","lc","bp","cn"))
Dor_ps.fish.rda.sites <- data.frame(Dor_ps.fish.rda.scores$sites)
Dor_ps.fish.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.fish.rda.sites))
sample_data(Dor_ps.prev.fishponds)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.fishponds)))
sample_data(Dor_ps.prev.fishponds)$Mic.Season <- as.character(sample_data(Dor_ps.prev.fishponds)$Mic.Season )
Dor_ps.fish.rda.sites <- Dor_ps.fish.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.fishponds)) %>% 
  mutate(Mic.Season = factor(Mic.Season, levels= c("Wet","Dry")))


#Draw biplots
Dor_ps.fish.rda.arrows<- Dor_ps.fish.rda.scores$biplot*5
colnames(Dor_ps.fish.rda.arrows)<-c("x","y")
Dor_ps.fish.rda.arrows <- as.data.frame(Dor_ps.fish.rda.arrows)
Dor_ps.fish.rda.evals <- 100 * summary(Dor_ps.fish.rda.0)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.fish.plot <- ggplot() +
  geom_point(data = Dor_ps.fish.rda.sites, aes(x = RDA1, y = RDA2, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor_ps.fish.rda.sites, aes(x = RDA1, y = RDA2, colour = Mic.Season, shape = location), 
             size = 3) +
  #geom_text(data = Dor_ps.fish.rda.sites,aes(x = RDA1, y = RDA2,label = paste(Month,gsub("20","",Year))), 
  #          nudge_y= -0.8,size=5)+
  scale_colour_manual(values = c("Wet"="darkblue",
                                 "Dry"="orange")) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.fish.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.fish.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.fish.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50")+
  geom_text(data=as.data.frame(Dor_ps.fish.rda.arrows*2),
            aes(x, y, label = rownames(Dor_ps.fish.rda.arrows)),color="gray50", size = 4)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#generate an RDA plot of species
Dor_ps.fish.rda.species <- data.frame(Dor_ps.fish.rda.scores$species)
Dor_ps.fish.rda.species$ASV <- as.character(rownames(Dor_ps.fish.rda.species))
Dor_tax <- data.frame(as(tax_table(Dor_ps.prev.fishponds), "matrix"), ASV = rownames(tax_table(Dor_ps.prev.fishponds)))
Dor_ps.fish.rda.species <- Dor_ps.fish.rda.species %>%
  left_join(Dor_tax, by= "ASV") 

Dor_ps.fish.rda.species.sub <- Dor_ps.fish.rda.species %>% mutate(Class =  ifelse(Class %in% Dor_ps.ra.long.agg$Class, Class, "Other taxa")) %>% 
  mutate(Class = factor(Class,levels=c(taxa_classes,"Other taxa")),
         Alpha = case_when(ASV %in% enriched_ASV$ASV ~ 1, 
                           TRUE ~ 0.5))

#Draw biplots
Dor_ps.fish.rda.arrows<- Dor_ps.fish.rda.scores$biplot
colnames(Dor_ps.fish.rda.arrows)<-c("x","y")
Dor_ps.fish.rda.arrows <- as.data.frame(Dor_ps.fish.rda.arrows)
Dor_ps.fish.rda.evals <- 100 * summary(Dor_ps.fish.rda.0)$cont$importance[2, c("RDA1","RDA2")]

#Plot 
Dor_ps.rda.fish.species.plot <- ggplot() +
  geom_point(data = subset(Dor_ps.fish.rda.species.sub,Alpha == 0.5), aes(x = RDA1, y = RDA2), alpha = 0.2,
             fill = "gray", size = 1) +
  geom_point(data = subset(Dor_ps.fish.rda.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2), alpha = 1,
             fill = "black", size = 5) +
  geom_point(data = subset(Dor_ps.fish.rda.species.sub, Alpha == 1), aes(x = RDA1, y = RDA2, colour = Class), alpha = 1,
             size = 3) +
  #geom_text(data = subset(Dor_ps.fish.rda.species.sub, Alpha == 1),aes(x = RDA1, y = RDA2,label = Genus), 
  #          nudge_y= -0.01,size=3)+
  scale_colour_manual(values = class_col) + 
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.fish.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.fish.rda.evals[2], 2))) +
  geom_segment(data=Dor_ps.fish.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="gray50",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.fish.rda.arrows*1.1),
            aes(x, y, label = rownames(Dor_ps.fish.rda.arrows)),color="gray50",alpha=0.5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

ggarrange(Dor_ps.rda.fish.plot, Dor_ps.rda.fish.species.plot, #heights = c(2,1.2),
          ncol = 2, nrow = 1, align = "h", legend = "bottom", labels ="AUTO")

ggsave("./figures/RDA_aqua.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
