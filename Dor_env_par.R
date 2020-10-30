#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")
library(rstatix); packageVersion("rstatix")
library(PerformanceAnalytics); packageVersion("PerformanceAnalytics")

#load colors
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#subset by fraction and remove NAs in metadata
sample_data(Dor_ps.prev.vst)[sample_data(Dor_ps.prev.vst) == "#N/A"]<- NA

Dor_ps.prev.no.na <- Dor_ps.prev.vst %>% 
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

#####################################
#Check correlations between environmental parameters
#####################################

env.par<- c("Temp_degC", "Chl_a_mg_L", "Chl_b_mg_L",
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L",
            "b_caroten_mg_L","MC_ug_L","Lutein_mg_L","Zeaxanthin_mg_L","Food_Kg_pond")

#scale parameters
metadata.scaled <- data.frame(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate_at(all_of(env.par),as.numeric) %>%
  #mutate_if(is.numeric, round, 3)%>%
  mutate_if(is.numeric, scale_par)


#correlation between the env parameters
envpar_corr <- metadata.scaled %>% 
  select(env.par)%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr.pvalues<- envpar_corr %>% cor_get_pval()

#plot
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)

#remove parameters that strongly correlate to each other 
metadata.scaled_sub<-metadata.scaled %>% 
  select(all_of(env.par))%>%
  select(-c(Dinoxanthin_mg_L,
            Diatoxanthin_mg_L))
#####################################
#Variation partitioning analysis
#####################################
#RDA analysis
Dor_ps.prev.rda.all <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ ., data = metadata.scaled_sub) # model including all variables 

Dor_ps.prev.rda.0 <- rda(t(otu_table(Dor_ps.prev.no.na)) ~ 1,
                         data = metadata.scaled_sub) # model containing only species matrix and intercept
Dor_ps.prev.rda.sel.os <- ordistep(Dor_ps.prev.rda.0, scope = formula(Dor_ps.prev.rda.all), direction = 'both') #stepwise selection

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
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot*1
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


ggsave("./figures/RDA_step_all.pdf", 
       plot = Dor_ps.rda.plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)




#FL
Dor_ps.rda.scores <- vegan::scores(Dor_ps.prev.rda.all,display=c("sp","wa","lc","bp","cn"))
Dor_ps.rda.sites <- data.frame(Dor_ps.rda.scores$sites)
Dor_ps.rda.sites$Sample.ID <- as.character(rownames(Dor_ps.rda.sites))
sample_data(Dor_ps.prev.no.na)$Sample.ID <- as.character(rownames(sample_data(Dor_ps.prev.no.na)))
sample_data(Dor_ps.prev.no.na)$Season <- as.character(sample_data(Dor_ps.prev.no.na)$Season )
Dor_ps.rda.sites <- Dor_ps.rda.sites %>%
  left_join(sample_data(Dor_ps.prev.no.na)) %>% 
  mutate(Season = factor(Season, levels= c("Winter","Spring",
                                           "Summer","Autumn")))

#Draw biplots
Dor_ps.rda.arrows<- Dor_ps.rda.scores$biplot*3
colnames(Dor_ps.rda.arrows)<-c("x","y")
Dor_ps.rda.arrows <- as.data.frame(Dor_ps.rda.arrows)
Dor_ps.rda.evals <- 100 * (Dor_ps.prev.rda.sel.os$CCA$eig / sum(Dor_ps.prev.rda.sel.os$CCA$eig))



BAC_FL.rda.species <- data.frame(Dor_ps.rda.scores$species)

BAC_FL.rda.species$ASV <- as.character(rownames(BAC_FL.rda.species))
BAC_FL.rda.species$Genus <-as.character(tax_table(Dor_ps.prev.no.na)[,c("Genus")])
BAC_FL.rda.species$Class <-as.character(tax_table(Dor_ps.prev.no.na)[,c("Class")])
BAC_FL.rda.species$Phylum <-as.character(tax_table(Dor_ps.prev.no.na)[,c("Phylum")])

#BAC_FL.rda.species<- BAC_FL.rda.species[BAC_FL.rda.species$Class== "Cyanobacteriia",]

BAC_FL.rda.plot <- ggplot() +
  geom_point(data = BAC_FL.rda.species, aes(x = RDA1, y = RDA2, fill = Class), 
             shape =21, size = 4) +
  geom_text(data = BAC_FL.rda.species,aes(x = RDA1, y = RDA2,label =Genus), 
            nudge_y= -0.03,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(Dor_ps.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(Dor_ps.rda.evals[2], 2))) +
  #scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  geom_segment(data=Dor_ps.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(Dor_ps.rda.arrows*1),
            aes(x, y, label = rownames(Dor_ps.rda.arrows)),color="black",alpha=0.5)+
theme_bw()+
  theme(legend.position = "bottom")

ggsave("./figures/dada2_CCA_Res_cyano.png", 
       plot = BAC_FL.rda.plot,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



#FL
BAC_FL.rda.0 <- rda (BAC_FL.otu ~ 1, data = BAC_FL.env) # model containing only species matrix and intercept
BAC_FL.rda.sel.os <- ordistep (BAC_FL.rda.0, scope = formula (BAC_FL.rda.all), direction = 'both') #stepwise selection

#PA
BAC_PA.rda.0 <- rda (BAC_PA.otu ~ 1, data = BAC_PA.env) # model containing only species matrix and intercept
BAC_PA.rda.sel.os <- ordistep (BAC_PA.rda.0, scope = formula (BAC_PA.rda.all), direction = 'both') #stepwise selection

#variance partitioning 
varpart(BAC_FL.otu, ~ Temperature, ~ ChlA, ~ Salinity, data = BAC_FL.env[,c("Temperature", "ChlA", "Salinity")])

varpart(BAC_PA.otu, ~ Salinity, ~ dNO3, ~ Temperature, ~ ChlA,data = BAC_PA.env[,c("Salinity", "Temperature", "ChlA","dNO3")])

