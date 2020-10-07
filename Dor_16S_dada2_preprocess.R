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
Dor_ps <- phyloseq(ASVs_tab, TAX_tab, meta)

#remove chloroplast and mitochondria reads
Dor_ps<- subset_taxa(Dor_ps, !Order == "Chloroplast" & !Family =="Mitochondria")

#subset only Reservoir
Dor_ps<- subset_samples(Dor_ps, location == "Res." & 
                          !comment %in% c("control.SP01","control.SP02","control.SP03", "duplicate batch 2","Ashraf did the PCR")&
                          !Sample_number_dada2 %in% c("95"))

Dor_ps<- prune_taxa(taxa_sums(Dor_ps)>0,Dor_ps)

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

ggsave("./figures/dada2_Res_libs_overview.pdf", 
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
rare$Run <- meta$Run[match(rare$site, meta$site)]
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
  facet_grid(~Run)+
  #xlim(0,1e5)+
  #ylim(0,5000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/dada2_Res_rarefactions.pdf", 
       plot = rare.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#####################################
# Compute prevalence of each feature
#####################################
prevdf <-  apply(X = otu_table(Dor_ps),
                 MARGIN = ifelse(taxa_are_rows(Dor_ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})

# # Add taxonomy and total read counts to this data.frame
prevdf.tax  <-  data.frame(Prevalence = prevdf,
                           TotalAbundance = taxa_sums(Dor_ps),
                           tax_table(Dor_ps))
#summarize
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# 
# #plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(Dor_ps),color=Phylum)) +
  # # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +  theme_bw() + theme(legend.position="none")

prev_plot_phyl
ggsave("./figures/dada2_Res_prev_plot_phyl.pdf", 
       plot = prev_plot_phyl,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.05 * nsamples(Dor_ps))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
Dor_ps.prev <-  prune_taxa((prevdf > prevalenceThreshold), Dor_ps)

