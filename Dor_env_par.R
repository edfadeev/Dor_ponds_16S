#PCA plot
Dor_ps.dds <- phyloseq_to_deseq2(Dor_ps, ~1)
geoMeans = apply(counts(Dor_ps.dds), 1, gm_mean)
Dor_ps.dds = estimateSizeFactors(Dor_ps.dds, geoMeans = geoMeans)
Dor_ps.dds <- estimateDispersions(Dor_ps.dds)
otu.vst <- getVarianceStabilizedData(Dor_ps.dds)

#make sure that the dimensions of the OTU table and the DEseq object are matching
dim(otu.vst)
dim(otu_table(Dor_ps))

Dor_ps.vst<-Dor_ps
otu_table(Dor_ps.vst)<- otu_table(otu.vst, taxa_are_rows = TRUE)




sample_data(Dor_ps.vst)[sample_data(Dor_ps.vst) == "#N/A"]<- NA


#subset by fraction and remove NAs in metadata
BAC_FL.no.na <- Dor_ps.vst %>% 
  subset_samples(
    !is.na(MC_ug_L) & 
      !is.na(Fucoxanthin_mg_L) &
      !is.na(Diatoxanthin_mg_L) &
      !is.na(b_caroten_mg_L) &
      !is.na(Ammonia_ug_L) &
      !is.na(Temp_degC)&
      !is.na(Chl_a_mg_L))

##remove unobserved OTU
BAC_FL.no.na <- prune_taxa(taxa_sums(BAC_FL.no.na)>0,BAC_FL.no.na)

#extract and scale the env. parameters
BAC_FL.env <- data.frame(sample_data(BAC_FL.no.na))[c("Temp_degC", "Chl_a_mg_L", "Ammonia_ug_L", 
                                                      "Fucoxanthin_mg_L", "MC_ug_L", "Diatoxanthin_mg_L",
                                                      "b_caroten_mg_L")]%>%
  mutate_if(is.character,as.numeric)

BAC_FL.env <- as.data.frame(scale(BAC_FL.env,center = FALSE, scale = TRUE))

BAC_FL.otu <- t(otu_table(BAC_FL.no.na))


BAC_FL.rda.all <- rda (BAC_FL.otu ~ ., data = BAC_FL.env) # model including all variables 

#generate an RDA plot 
#FL
BAC_FL.rda.scores <- vegan::scores(BAC_FL.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC_FL.rda.sites <- data.frame(BAC_FL.rda.scores$sites)
BAC_FL.rda.sites$Sample.ID <- as.character(rownames(BAC_FL.rda.sites))
sample_data(BAC_FL.no.na)$Sample.ID <- as.character(rownames(sample_data(BAC_FL.no.na)))
sample_data(BAC_FL.no.na)$Season <- as.character(sample_data(BAC_FL.no.na)$Season )
BAC_FL.rda.sites <- BAC_FL.rda.sites %>%
  left_join(sample_data(BAC_FL.no.na))

#Draw biplots
BAC_FL.rda.arrows<- BAC_FL.rda.scores$biplot*5
colnames(BAC_FL.rda.arrows)<-c("x","y")
BAC_FL.rda.arrows <- as.data.frame(BAC_FL.rda.arrows)
BAC_FL.rda.evals <- 100 * (BAC_FL.rda.all$CCA$eig / sum(BAC_FL.rda.all$CCA$eig))

#Plot 
BAC_FL.rda.plot <- ggplot() +
  geom_point(data = BAC_FL.rda.sites, aes(x = RDA1, y = RDA2, fill = Season), 
             shape =21, size = 4) +
  geom_text(data = BAC_FL.rda.sites,aes(x = RDA1, y = RDA2,label = paste(Year,"-",Month)), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_FL.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_FL.rda.evals[2], 2))) +
  #scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  theme(legend.position = "none")+
  geom_segment(data=BAC_FL.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_FL.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_FL.rda.arrows)),color="black",alpha=0.5)





BAC_FL.rda.species <- data.frame(BAC_FL.rda.scores$species)

BAC_FL.rda.species$ASV <- as.character(rownames(BAC_FL.rda.species))
BAC_FL.rda.species$Genus <-as.character(tax_table(BAC_FL.no.na)[,c("Genus")])
BAC_FL.rda.species$Class <-as.character(tax_table(BAC_FL.no.na)[,c("Class")])
BAC_FL.rda.species$Phylum <-as.character(tax_table(BAC_FL.no.na)[,c("Phylum")])

BAC_FL.rda.species<- BAC_FL.rda.species[BAC_FL.rda.species$Phylum== "Cyanobacteria",]

BAC_FL.rda.plot <- ggplot() +
  geom_point(data = BAC_FL.rda.species, aes(x = RDA1, y = RDA2, fill = Phylum), 
             shape =21, size = 4) +
  geom_text(data = BAC_FL.rda.species,aes(x = RDA1, y = RDA2,label =Genus), 
            nudge_y= -0.03,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_FL.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_FL.rda.evals[2], 2))) +
  #scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  theme(legend.position = "bottom")+
  geom_segment(data=BAC_FL.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_FL.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_FL.rda.arrows)),color="black",alpha=0.5)

