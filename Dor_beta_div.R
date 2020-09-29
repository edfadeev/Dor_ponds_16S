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

BAC_pruned.ra.long$Family <- as.character(BAC_pruned.ra.long$Family)
BAC_pruned.ra.long$Family[is.na(BAC_pruned.ra.long$Family)] <- paste(BAC_pruned.ra.long$Order[is.na(BAC_pruned.ra.long$Family)],"uc", sep = "_")

#calculate abundance for each Class
BAC_pruned.ra.long %>% select(Month,Year,OTU,Class,Abundance)%>%
  group_by(Month,Year,Class) %>%
  summarize(Abund.total= sum(Abundance)) -> BAC_pruned.ra.long.agg

#remove below 1% ra
taxa_classes <- unique(BAC_pruned.ra.long.agg$Class)
BAC_pruned.ra.long.agg$Class[BAC_pruned.ra.long.agg$Abund.total<1] <- "Other taxa"
BAC_pruned.ra.long.agg$Class <- factor(BAC_pruned.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))

#Plot 
barplots <- ggplot(BAC_pruned.ra.long.agg, aes(x = Month, y = Abund.total, fill = Class)) + 
  facet_grid(.~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = sample(col_vector)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme(legend.position="bottom")








#####################################
#PCA plot
PS107_merged.dds <- phyloseq_to_deseq2(Dor_ps.prev, ~1)
varianceStabilizingTransformation(PS107_merged.dds, blind = TRUE, fitType = "parametric")
PS107_merged.dds <- estimateSizeFactors(PS107_merged.dds)
PS107_merged.dds <- estimateDispersions(PS107_merged.dds)
otu.vst <- getVarianceStabilizedData(PS107_merged.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(Dor_ps.prev))

Dor_ps.prev.vst<-Dor_ps.prev
otu_table(Dor_ps.prev.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)



PS107.ord <- ordinate(Dor_ps.prev.vst, method = "RDA", distance = "eucledian")
PS107.ord.df <- plot_ordination(Dor_ps.prev.vst, PS107.ord, axes = c(1,2,3),justDF = TRUE)

#extract explained variance
PS107.ord.evals <- 100 * (PS107.ord$CA$eig/ sum(PS107.ord$CA$eig))
PS107.ord.df$ID <- rownames(PS107.ord.df)

PS107.ord.p <- ggplot(data = PS107.ord.df, aes(x =PC1, y=PC2, shape = Year, colour = Month))+
  geom_point(colour = "black", size = 5)+
  geom_point(size = 4)+
  scale_fill_manual(values = sample(col_vector)) + 
  #stat_ellipse(data = PS107.ord.df, aes(x =PC1, y=PC2, group = Year), size = 1, type= "norm")+
  geom_text(aes(label = Month), colour = "black", nudge_y= -0.5,  size=4)+
  labs(x = sprintf("PC1 [%s%%]", round(PS107.ord.evals[1], 2)), 
       y = sprintf("PC2 [%s%%]", round(PS107.ord.evals[2], 2)), shape = "Year")+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")


#significance test
df <- as(sample_data(Dor_ps.prev.vst), "data.frame") %>% drop_na()
Dor_ps.prev.vst_sub <- prune_samples(row.names(df), Dor_ps.prev.vst)
d <- phyloseq::distance(Dor_ps.prev.vst_sub, "euclidean")
adonis_all <- adonis(d ~ Year + Month + Chl_a + MC_ug_ml + Temp + Ammonia_ug_L + NO3_NO2 + TP_ug_L, df)
adonis_all

