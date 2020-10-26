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
ASVs_tab<- read.csv("./dada2/dada2_seqtab_nochim2.txt", h=T, sep="\t")
#simplify sample names to match the metadata
colnames(ASVs_tab)<- gsub("_F_filt.fastq.gz","",colnames(ASVs_tab))

TAX_tab<- as.matrix(read.csv("./dada2/dada2_taxonomy_table.txt", h=T,sep = "\t"))

ENV <- read.csv("./dada2/sample_list.csv", sep = "," , h = T, row.names = 1, fill = T, na.strings=c("","NA"))
#remove samples that did not pass dada2 workflow
ENV<- ENV[names(ASVs_tab),] 
#add seasons
ENV<- ENV %>% 
  mutate(Season = case_when(month %in% c("12", "1", "2") ~ "Winter",
                            month %in% c("3", "4", "5") ~ "Spring",
                            month %in% c("6", "7", "8") ~ "Summer",
                            month %in% c("9", "10", "11") ~ "Autumn"),
         Month = case_when(month == 1 ~ "Jan",month == 2 ~ "Feb",month == 3 ~ "Mar",
                           month == 4 ~ "Apr",month == 5 ~ "May",month == 6 ~ "Jun",
                           month == 7 ~ "Jul",month == 8 ~ "Aug",month == 9 ~ "Sep",
                           month == 10 ~ "Oct",month == 11 ~ "Nov",month == 12 ~ "Dec")) %>% 
  mutate(Season = factor(Season, levels = c("Winter","Spring","Summer","Autumn")),
         Year = factor(year, levels = c("2013","2014","2015")),
         Month = ifelse(cod == "Res.5Jul.", "Jul", Month),#<<<---verify that with Ashraf and Sofi!!!
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Day = day) %>% 
          select(-month,year,day)

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

# subset only 2013-2014
Dor_ps0_run1<- subset_samples(Dor_ps0, Run == "1")

#Remove controls and etc
Dor_ps<- subset_samples(Dor_ps0_run1, location %in% c("Res.","V2.","D1.") & 
                          !comment %in% c("control.SP01","control.SP02","control.SP03", "duplicate batch 2","Ashraf did the PCR")&
                          !Sample_number_dada2 %in% c("95"))

Dor_ps<- prune_taxa(taxa_sums(Dor_ps)>0,Dor_ps)


#remove samples with less than 2000 sequences
Dor_ps.prev<-  prune_samples(sample_sums(Dor_ps)>5000, Dor_ps)

saveRDS(Dor_ps.prev, "./data/Dor_ps_prev.rds")

#####################################
#Plot total number of reads and OTUs per sample
#####################################
readsumsdf <- data.frame(nreads = sort(taxa_sums(Dor_ps), TRUE), sorted = 1:ntaxa(Dor_ps), 
                         type = "OTUs")
readsumsdf <- rbind(readsumsdf, data.frame(nreads = sort(sample_sums(Dor_ps), 
                                                         TRUE), sorted = 1:nsamples(Dor_ps), type = "Samples"))
title = "Total number of reads"
p <- ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
OTU_libs_overview <-  p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")+ 
  xlab("Samples")+ ylab("Number of reads")+ theme_bw()

ggsave("./figures/libs_overview_run1.pdf", 
       plot = OTU_libs_overview,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
#Plot rarefaction
####################################
iNEXT.out <- iNEXT(as.data.frame(otu_table(Dor_ps)), q=0, datatype="abundance")
rare <-fortify(iNEXT.out, type=1)

meta <- as(sample_data(Dor_ps), "data.frame")
meta$site <- rownames(meta)
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
  geom_text(aes(label=label), size =2, data= rare.point, colour = "black",  nudge_x = 10000)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_grid(~Loc)+
  #xlim(0,1e5)+
  #ylim(0,5000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/Run1_rarefactions.pdf", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)