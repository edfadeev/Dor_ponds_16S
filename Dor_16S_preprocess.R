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
                         h=T,sep = ";", row.names = 1))
ENV <- read.csv("/Users/eduardfadeev/Google Drive (dr.eduard.fadeev@gmail.com)/DSher-projects/Dor_metaG/16S/Dor_metadata.csv", 
                sep = "," , h = T, row.names = "OUT")

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
