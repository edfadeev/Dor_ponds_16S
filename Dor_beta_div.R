#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

tol21rainbow<- c("#771155", 
                 "#AA4488", 
                 "#CC99BB", 
                 "#114477", 
                 "#4477AA", 
                 "#77AADD", 
                 "#117777", 
                 "#44AAAA", 
                 "#77CCCC", 
                 "#117744", 
                 "#44AA77", 
                 "#88CCAA", 
                 "#777711", 
                 "#AAAA44", 
                 "#DDDD77", 
                 "#774411", 
                 "#AA7744", 
                 "#DDAA77", 
                 "#771122", 
                 "#AA4455", 
                 "#DD7788")
#####################################
#Plot barplots of communities
#####################################

#transform data
BAC_pruned.ra <- transform_sample_counts(Dor_ps.prev, function(x) x / sum(x))
BAC_pruned.ra.long <- psmelt(BAC_pruned.ra)
BAC_pruned.ra.long$Abundance <- BAC_pruned.ra.long$Abundance*100

#fix unclassified lineages 
BAC_pruned.ra.long$Class <- as.character(BAC_pruned.ra.long$Class)
BAC_pruned.ra.long$Class[is.na(BAC_pruned.ra.long$Class)] <- paste(BAC_pruned.ra.long$Phylum[is.na(BAC_pruned.ra.long$Class)],"uc", sep = "_")

BAC_pruned.ra.long$Order <- as.character(BAC_pruned.ra.long$Order)
BAC_pruned.ra.long$Order[is.na(BAC_pruned.ra.long$Order)] <- paste(BAC_pruned.ra.long$Class[is.na(BAC_pruned.ra.long$Order)],"uc", sep = "_")

BAC_pruned.ra.long$Species <- as.character(BAC_pruned.ra.long$Species)
BAC_pruned.ra.long$Species[is.na(BAC_pruned.ra.long$Species)] <- paste(BAC_pruned.ra.long$Class[is.na(BAC_pruned.ra.long$Species)],"uc", sep = "_")

#calculate abundance for each Class
BAC_pruned.ra.long.agg <- BAC_pruned.ra.long %>% 
  select(Month,Year,OTU,Class,Abundance)%>%
  group_by(Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

#remove below 1% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)
BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abund.total<2] <- "Other taxa"
BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

BAC_pruned.ra.long.agg$Class<- droplevels(BAC_pruned.ra.long.agg$Class)

#Plot 
barplots_total <- ggplot(BAC_pruned.ra.long.agg, aes(x = Month, y = Abund.total, fill = Class)) + 
  facet_grid(.~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol21rainbow) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom")


ggsave("./figures/dada2_Res_barplots_total.pdf", 
       plot = barplots_total,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#calculate abundance for each Class
BAC_pruned_Cyanos <- BAC_pruned.ra.long %>% 
  select(Month,Year,OTU,Class, Family, Genus, Species, Abundance)%>%
  filter(Class == "Cyanobacteriia") %>% 
  group_by(Year,Month,Family, Genus, Species) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

barplots_cyanos<- ggplot(BAC_pruned_Cyanos, aes(x = Month, y = Abund.total, fill = Family))+
  geom_col()+
  ylim(0,100)+
  facet_grid(.~Year, space= "fixed")+
  scale_fill_manual(values = tol21rainbow)+
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom")

ggsave("./figures/dada2_Res_barplots_cyanos.pdf", 
       plot = barplots_cyanos,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#PCA plot
Dor_ps.dds <- phyloseq_to_deseq2(Dor_ps.prev, ~1)
geoMeans = apply(counts(Dor_ps.dds), 1, gm_mean)
Dor_ps.dds = estimateSizeFactors(Dor_ps.dds, geoMeans = geoMeans)
Dor_ps.dds <- estimateDispersions(Dor_ps.dds)
otu.vst <- getVarianceStabilizedData(Dor_ps.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(Dor_ps.prev))

Dor_ps.prev.vst<-Dor_ps.prev
otu_table(Dor_ps.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)



Dor_ps.prev.ord <- ordinate(Dor_ps.prev.vst, method = "RDA", distance = "eucledian")
Dor_ps.prev.ord.df <- plot_ordination(Dor_ps.prev.vst, Dor_ps.prev.ord, axes = c(1,2,3),justDF = TRUE)

#extract explained variance
Dor_ps.prev.ord.evals <- 100 * (Dor_ps.prev.ord$CA$eig/ sum(Dor_ps.prev.ord$CA$eig))
Dor_ps.prev.ord.df$ID <- rownames(Dor_ps.prev.ord.df)

Dor_ps.prev.ord.p <- ggplot(data = Dor_ps.prev.ord.df, aes(x =PC1, y=PC2, colour = Season))+
  geom_point(colour = "black", size = 5)+
  geom_point(size = 4)+
  scale_colour_manual(values = tol21rainbow) + 
  stat_ellipse(data = Dor_ps.prev.ord.df, aes(x =PC1, y=PC2, group = Season, colour = Season), size = 1, type= "t")+
  geom_text(aes(label = paste(Year,"-",Month)), colour = "black", nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(Dor_ps.prev.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(Dor_ps.prev.ord.evals[2], 2)), shape = "Season")+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

ggsave("./figures/dada2_RDA_Res.pdf", 
       plot = Dor_ps.prev.ord.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#significance test
df <- as(sample_data(Dor_ps.prev.vst), "data.frame")
d <- phyloseq::distance(Dor_ps.prev.vst, "euclidean")
adonis_all <- adonis(d ~ Year + Season + Month , df)
adonis_all

