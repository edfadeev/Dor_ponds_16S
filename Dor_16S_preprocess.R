#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(tidyr); packageVersion("tidyr")
library(DESeq2); packageVersion("DESeq2")
library(ggplot2); packageVersion("ggplot2")

#####################################
#Parse for Phyloseq
#####################################
OTU<- read.csv("/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DSher-projects/Dor_metaG/16S/Dor_OTU_table.csv", 
               h=T, sep=",", row.names = "OUT")
TAX<- as.matrix(read.csv("/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DSher-projects/Dor_metaG/16S/Dor_16S_tax.csv",
                         h=T, sep = ";", row.names = 1))
ENV <- read.csv("/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DSher-projects/Dor_metaG/16S/Dor_metadata.csv", 
                sep = "," , h = T)
#add metadata from names
ENV<- ENV %>% mutate(Month = factor(gsub("(.*)([0-9])(...)(\\.)","\\3",OUT),
                                       levels =c("Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")),
                    Year = gsub("(Res.)([0-9]).*","201\\2",OUT)) 
#add row.names
row.names(ENV)<- ENV$OUT

#extract only the samples with the metadata
OTU<- OTU %>% select(rownames(ENV))


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
OTU_libs_overview <-  p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

OTU_libs_overview
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
prevdf.tax.summary <- plyr::ddply(prevdf.tax, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# 
# #plot
prev_plot_phyl <- ggplot(prevdf.tax, aes(TotalAbundance, Prevalence / nsamples(Dor_ps),color=Class)) +
  # # Include a guess for parameter
  geom_hline(yintercept = 0.15, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Class) + theme(legend.position="none")

prev_plot_phyl

#  Define prevalence threshold as 5% of total samples
prevalenceThreshold <- round(0.2 * nsamples(Dor_ps))
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
Dor_ps.prev <-  prune_taxa((prevdf > prevalenceThreshold), Dor_ps)
