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
env.par<- c("Temp_degC", "Chl_a_mg_L", "Chl_b_mg_L",
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L","Fish_biomass_g_pond",
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
  select(Temp_degC, #Chl_a_mg_L,
         Ammonia_ug_L, NO3_NO2_N_L,TP_ug_L)

#subset pigments and mycrocysteins
pig_par.scaled<-metadata.scaled %>% 
  select(Chl_a_mg_L, Chl_b_mg_L,
         Fucoxanthin_mg_L,
         b_caroten_mg_L,Lutein_mg_L,#MC_ug_L,Dinoxanthin_mg_L,
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
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.all$CCA$eig))

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
  labs(x = sprintf("RDA1 [%s%% of explained variance]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%% of explained variance]", round(Dor_ps.rda.evals[2], 2))) +
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
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.all$CCA$eig))

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
  labs(x = sprintf("RDA1 [%s%% of explained variance]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%% of explained variance]", round(Dor_ps.rda.evals[2], 2))) +
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
#Effect of aquaculture
#####################################
#subset phys. parameters
Fish_par.scaled<-metadata.scaled %>% 
  select(Food_Kg_pond, Fish_biomass_g_pond)


#RDA analysis
Dor_ps.pr
ev.rda.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., Fish_par.scaled) # model including all variables 

Dor_ps.prev.rda.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ Food_Kg_pond+Fish_biomass_g_pond,
                         data = data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
                           mutate_at(all_of(env.par),as.numeric)) # model containing only species matrix and intercept


#generate an RDA plot 
Dor_ps.rda.scores <- vegan::scores(Dor_ps.prev.rda.0,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.sites <- data.frame(Dor_ps.rda.scores$sites)
Dor_ps.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Season <- as.character(sample_data(Dor_ps.prev.no.na)$Season )
Dor_ps.rda.sites <- Dor_ps.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate(Season = factor(Season, levels= c("Winter","Spring",
                                           "Summer","Autumn")))


#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.0$CCA$eig / sum(Dor_ps.prev.rda.all$CCA$eig))

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
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*2),
            aes(x, y, label = rownames(Dor_ps.rda.arrows)),color="black",alpha=0.5)+
  theme_bw()+
  theme(legend.position = "bottom")



