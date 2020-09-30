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
OTU<- read.csv("QIIME/Dor_OTUs_SM.csv", 
               h=T, sep=";", row.names = "OUT")
TAX<- as.matrix(read.csv("QIIME/Dor_16S_tax_full.csv",
                         h=T, sep = ",", row.names = "OUT"))
ENV <- read.csv("QIIME/Dor_metadata.csv", 
                sep = "," , h = T)
#add metadata from names
ENV<- ENV %>% mutate(Month = factor(gsub("(.*)([0-9])(...)(\\.)","\\3",OUT),
                                       levels =c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")),
                    Year = gsub("(Res.)([0-9]).*","201\\2",OUT),
                    Season = case_when(Month %in% c("Dec", "Jan", "Feb") ~ "Winter",
                                       Month %in% c("Mar", "Apr", "May") ~ "Spring",
                                       Month %in% c("Jun", "Jul", "Aug") ~ "Summer",
                                       Month %in% c("Sep", "Oct", "Nov") ~ "Autumn")) %>% 
                mutate(Season = factor(Season, levels = c("Winter","Spring","Summer","Autumn")))


#add row.names
row.names(ENV)<- ENV$OUT

#extract only the samples with the metadata
OTU<- OTU %>% dplyr::select(rownames(ENV))

#extract only the overlapping OTU
OTU<- OTU[intersect(rownames(OTU), rownames(TAX)),]
TAX<- TAX[intersect(rownames(OTU), rownames(TAX)),]

# Check order of samples
all.equal(rownames(OTU), rownames(TAX))

#creating Phyloseq dataset
OTU <- otu_table(OTU, taxa_are_rows = TRUE)
TAX <- tax_table(TAX)
meta <- sample_data(ENV)
Dor_ps <- phyloseq(OTU, TAX, meta)

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

ggsave("./figures/OTU_libs_overview.pdf", 
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
rare$Year <- meta$Year[match(rare$site, meta$site)] 
rare$Month <- meta$Month[match(rare$site, meta$site)] 
rare$Type <- meta$Type[match(rare$site, meta$site)]
rare$label <- paste(rare$Year,rare$Month, sep = "-")

rare.point <- rare[which(rare$method == "observed"),]
rare.line <- rare[which(rare$method != "observed"),]
rare.line$method <- factor (rare.line$method,
                            c("interpolated", "extrapolated"),
                            c("interpolation", "extrapolation"))

rare.p <- ggplot(rare, aes(x=x, y=y, colour = site))+
  geom_line(aes(linetype = method), lwd = 0.5, data= rare.line)+
  #geom_ribbon(aes(ymin=y.lwr, ymax= y.upr, colour = NULL), alpha = 0.2)+
  geom_point(aes(shape=Year), size =3, data= rare.point)+
  geom_text(aes(label=label), size =2, data= rare.point, colour = "black",  nudge_x = 10000)+
  scale_colour_discrete(guide = FALSE)+
  labs(x = "Sample size", y = "Species richness")+
  facet_grid(~Year)+
  xlim(0,1e5)+
  ylim(0,5000)+
  theme_classic(base_size = 12)+theme(legend.position="bottom")

ggsave("./figures/rarefactions.pdf", 
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
  geom_hline(yintercept = 0.2, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) +  theme_bw() + theme(legend.position="none")

prev_plot_phyl
ggsave("./figures/prev_plot_phyl.pdf", 
       plot = prev_plot_phyl,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.2 * nsamples(Dor_ps))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
Dor_ps.prev <-  prune_taxa((prevdf > prevalenceThreshold), Dor_ps)

