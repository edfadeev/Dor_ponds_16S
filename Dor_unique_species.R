#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("cowplot"); packageVersion("cowplot")
library("reshape2"); packageVersion("reshape2")
library("rstatix"); packageVersion("rstatix")
library("ggpubr"); packageVersion("ggpubr")
library("UpSetR"); packageVersion("UpSetR")

#load colors
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Unique ASVs in each pool
####################################
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#subset each pool
Dor_ps.D1<- subset_samples(Dor_ps.prev_run1, location =="D1.")
Dor_ps.D1<- prune_taxa(taxa_sums(Dor_ps.D1)>0,Dor_ps.D1)


Dor_ps.V2<- subset_samples(Dor_ps.prev_run1, location =="V2.")
Dor_ps.V2<- prune_taxa(taxa_sums(Dor_ps.V2)>0,Dor_ps.V2)

Dor_ps.Res<- subset_samples(Dor_ps.prev_run1, location =="Res.")
Dor_ps.Res<- prune_taxa(taxa_sums(Dor_ps.Res)>0,Dor_ps.Res)

#generate list of ASVs in each pool
z <- list()
z[["D1"]] <- as.character(row.names(otu_table(Dor_ps.D1)))
z[["V2"]] <- as.character(row.names(otu_table(Dor_ps.V2)))
z[["Res"]] <- as.character(row.names(otu_table(Dor_ps.Res)))

#generate overlap matrix
ASV_overlaps <- pres_abs_matrix(z)    
ASV_overlaps$ASV <- rownames(ASV_overlaps)

#split to unique ASVs of each pool
D1_unique_ASV <- ASV_overlaps %>% filter(D1 == 1, V2 == 0, Res == 0)
V2_unique_ASV <- ASV_overlaps %>% filter(D1 == 0, V2 == 1, Res == 0)
Res_unique_ASV <- ASV_overlaps %>% filter(D1 == 0, V2 == 0, Res == 1)

#extract ASVs that are present in both fishpond
fishponds_ASV<- ASV_overlaps %>% filter(Res == 0)

####################################
#Explore taxonomy of the unique ASVs
####################################
#calculate proportions
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev_run1, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)
Dor_ps.ra.long$Abundance <- Dor_ps.ra.long$Abundance*100

#fix unclassified lineages 
Dor_ps.ra.long$Class <- as.character(Dor_ps.ra.long$Class)
Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Class)] <- paste(Dor_ps.ra.long$Phylum[is.na(Dor_ps.ra.long$Class)],"uc", sep = "_")

Dor_ps.ra.long$Genus <- as.character(Dor_ps.ra.long$Genus)
Dor_ps.ra.long$Genus[is.na(Dor_ps.ra.long$Genus)] <- paste(Dor_ps.ra.long$Family[is.na(Dor_ps.ra.long$Genus)],"uc", sep = "_")

Dor_ps.ra.long$Species <- as.character(Dor_ps.ra.long$Species)
Dor_ps.ra.long$Species[is.na(Dor_ps.ra.long$Species)] <- paste(Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Species)],"uc", sep = "_")

#subset Abundance of unique ASVs in each pool and aggregate to a genus level
D1_unique_ASV_abund<- Dor_ps.ra.long %>% filter(OTU %in% D1_unique_ASV$ASV, location == "D1.") %>% 
  select(location, Month,Year,OTU,Class,Order, Family, Genus, Species, Abundance)%>%
  group_by(location, Year,Month,Class,Order, Family, Genus, Species) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  filter(Abund.total>0)

V2_unique_ASV_abund<- Dor_ps.ra.long %>% filter(OTU %in% V2_unique_ASV$ASV, location == "V2.")%>% 
  select(location, Month,Year,OTU,Class,Order, Family, Genus, Species, Abundance)%>%
  group_by(location, Year,Month,Class,Order, Family, Genus, Species) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  filter(Abund.total>0)

Res_unique_ASV_abund<- Dor_ps.ra.long %>% filter(OTU %in% Res_unique_ASV$ASV, location == "Res.")%>% 
  select(location, Month,Year,OTU,Class,Order, Family, Genus, Species, Abundance)%>%
  group_by(location, Year,Month,Class,Order, Family, Genus, Species) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  filter(Abund.total>0)

#merge together for plotting
unique_ASV_abund<- rbind(D1_unique_ASV_abund, V2_unique_ASV_abund, Res_unique_ASV_abund) %>% 
  filter(Abund.total>0.5)

#plot
unique_bar.p <- ggplot(unique_ASV_abund, aes(x = Month, y = Abund.total, fill = Class)) + 
  facet_grid(location~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol24rainbow) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (> 0.5 %) \n")+
  theme_bw()+
  theme(legend.position="bottom")


#export table
unique_species_by_pool <- rbind(D1_unique_ASV_abund, V2_unique_ASV_abund, Res_unique_ASV_abund) %>% spread(location, Abund.total)
write.csv(unique_species_by_pool, "tables/unique_species_by_pool.csv")


#explore fishponds unique ASVs
Fishpond_unique_ASV_abund<- Dor_ps.ra.long %>% filter(OTU %in% fishponds_ASV$ASV)%>% 
  select(location, Month,Year,OTU,Class,Order, Family, Genus, Species, Abundance)%>%
  group_by(location, Year,Month,Class,Order, Family, Genus, Species) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) %>% 
  filter(Abund.total>0)%>% spread(location, Abund.total)

write.csv(Fishpond_unique_ASV_abund, "tables/Fishpond_species.csv")

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
