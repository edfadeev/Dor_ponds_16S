#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(ggpubr); packageVersion("ggpubr")


#load colours
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Plot cyanobacterial families abundance
#####################################
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

ggsave("./figures/barplots_cyanos.pdf", 
       plot = barplots_cyanos,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))

