#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(tidyr); packageVersion("tidyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(iNEXT); packageVersion("iNEXT")

#####################################
#Parse for Phyloseq
#####################################
ASVs_tab<- read.csv("./dada2/ASV_tab.csv", h=T, row.names = 1)

TAX_tab<- as.matrix(read.csv("./dada2/tax_tab.csv", h=T, row.names = 1))

ENV <- read.csv("./dada2/meta_table.csv", h = T, row.names = 1)

ENV<- ENV %>% mutate(Mic.Season = factor(Mic.Season, levels = c("Wet","Dry")),
                     Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                                      "May","Jun","Jul","Aug",
                                                      "Sep","Oct","Nov","Dec")))
rownames(ENV)<- paste0("X",ENV$Sample_number_dada2)
  
# Check order of samples
all.equal(rownames(ASVs_tab), rownames(TAX_tab))

#creating Phyloseq dataset
ASVs_tab <- otu_table(ASVs_tab, taxa_are_rows = TRUE)
TAX_tab <- tax_table(TAX_tab)
meta <- sample_data(ENV)
Dor_ps0 <- phyloseq(ASVs_tab, TAX_tab, meta)

#add reference sequence and replace variants with ASVs
dna <- Biostrings::DNAStringSet(taxa_names(Dor_ps0))
names(dna) <- taxa_names(Dor_ps0)
Dor_ps0 <- merge_phyloseq(Dor_ps0, dna)
taxa_names(Dor_ps0) <- paste0("ASV", seq(ntaxa(Dor_ps0)))

#remove chloroplast and mitochondria reads
Dor_ps0.chl<- subset_taxa(Dor_ps0, Order == "Chloroplast")
Dor_ps0.mit<- subset_taxa(Dor_ps0, Family == "Mitochondria")

Dor_ps0<- subset_taxa(Dor_ps0, !Order == "Chloroplast" & !Family =="Mitochondria")

#remove samples with less than 2000 sequences
Dor_ps0<-  prune_samples(sample_sums(Dor_ps0)>2000, Dor_ps0)

Dor_ps.prev<- prune_taxa(taxa_sums(Dor_ps0)>0,Dor_ps0)

saveRDS(Dor_ps.prev, "./data/Dor_ps_prev.rds")

#####################################
#Plot total number of reads and OTUs per sample
#####################################
readsumsdf <- data.frame(nreads = sort(taxa_sums(Dor_ps.prev), TRUE), sorted = 1:ntaxa(Dor_ps.prev), 
                         type = "OTUs")
readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(Dor_ps.prev), 
                                                         TRUE), sorted = 1:nsamples(Dor_ps.prev), type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
OTU_libs_overview <-  p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")+ 
  xlab("Samples")+ ylab("Number of reads")+ theme_bw()

ggsave("./figures/libs_overview.png", 
       plot = OTU_libs_overview,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Plot rarefaction
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(Dor_ps.prev)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(Dor_ps.prev), "data.frame")
meta$site <- rownames(meta)
rare$Run <- meta$Run[match(rare$site, meta$site)]
rare$Loc <- meta$location[match(rare$site, meta$site)]
rare$Year <- meta$Year[match(rare$site, meta$site)] 
rare$Month <- meta$Month[match(rare$site, meta$site)] 
rare$label <- paste(rare$Year,rare$Month, sep = "-")

rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(size =3, data= rare.point)+
  geom_text(aes(label=label), size =4, data= rare.point, colour = "black",  nudge_y = 10, nudge_x = 1000)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_grid(~Loc)+
  xlim(0,1e5)+
  #ylim(0,5000)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


ggsave("./figures/rarefactions.pdf", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))
