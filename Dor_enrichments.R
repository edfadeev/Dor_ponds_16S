#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(tidyr); packageVersion("tidyr")
library(DESeq2); packageVersion("DESeq2")
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
#Split by seasons
#####################################
Dor_winter_spring <- subset_samples(Dor_ps.prev, Season %in% c("Winter","Spring"))
Dor_winter_spring <- prune_taxa(taxa_sums(Dor_winter_spring)>0,Dor_winter_spring)


Dor_spring_summer <- subset_samples(Dor_ps.prev, Season %in% c("Spring","Summer"))
Dor_spring_summer <- prune_taxa(taxa_sums(Dor_spring_summer)>0,Dor_spring_summer)


Dor_summer_autumn <- subset_samples(Dor_ps.prev, Season %in% c("Summer","Autumn"))
Dor_summer_autumn <- prune_taxa(taxa_sums(Dor_summer_autumn)>0,Dor_summer_autumn)


Dor_autumn_winter <- subset_samples(Dor_ps.prev, Season %in% c("Winter","Autumn"))
Dor_autumn_winter <- prune_taxa(taxa_sums(Dor_autumn_winter)>0,Dor_autumn_winter)


Dor_winter_summer <- subset_samples(Dor_ps.prev, Season %in% c("Winter","Summer"))
Dor_winter_summer <- prune_taxa(taxa_sums(Dor_winter_summer)>0,Dor_winter_summer)

#####################################
#run DEseq2
#####################################
type <- c("Dor_winter_spring","Dor_spring_summer","Dor_summer_autumn","Dor_autumn_winter", "Dor_winter_summer")
deseq_res_all <- data.frame()
enriched_agg_all <- data.frame()

for (i in 1:5){
    #run DEseq
    BAC_sed.ddsMat <- phyloseq_to_deseq2(get(type[i]), ~Season)
    geoMeans = apply(counts(BAC_sed.ddsMat), 1, gm_mean)
    BAC_sed.ddsMat = estimateSizeFactors(BAC_sed.ddsMat, geoMeans = geoMeans)
    BAC_sed.ddsMat <- estimateDispersions(BAC_sed.ddsMat)
    BAC_sed.DEseq <- DESeq(BAC_sed.ddsMat, fitType="local")
    BAC_sed.DEseq.res <- results(BAC_sed.DEseq)
    
    #extract only significant OTU
    BAC_sed.DEseq.res.sig <- BAC_sed.DEseq.res[which(BAC_sed.DEseq.res$padj < 0.1), ]
    BAC_sed.DEseq.res.sig <- cbind(as(BAC_sed.DEseq.res.sig, "data.frame"),
                                   as(tax_table(get(type[i]))[rownames(BAC_sed.DEseq.res.sig), ], "matrix"))
    BAC_sed.DEseq.res.sig.sub <- BAC_sed.DEseq.res.sig[abs(BAC_sed.DEseq.res.sig$log2FoldChange)>1,]
    BAC_sed.DEseq.res.sig.sub$merging <- type[i]
    BAC_sed.DEseq.res.sig.sub$ASV <- rownames(BAC_sed.DEseq.res.sig.sub)
    rownames(BAC_sed.DEseq.res.sig.sub) <- c()
    deseq_res_all <- rbind(deseq_res_all,BAC_sed.DEseq.res.sig.sub, make.row.names = TRUE)
  }

#####################################
#Explore the results
#####################################
#separate the results by enrichment towards surface and to depth
enriched_shallow <- deseq_res_all[deseq_res_all[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Species","ASV","merging") ]
enriched_deep  <- deseq_res_all[deseq_res_all[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Species","ASV","merging")]

#Aggregate on Order level 
enriched_shallow.agg <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family+Species,
                                                       enriched_shallow, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
enriched_deep.agg <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family+Species,
                                                    enriched_deep, 
                                                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#merge data for ploting
enriched_agg <- rbind(enriched_shallow.agg,enriched_deep.agg)

#exclude taxa with less than 3 ASV
enriched_agg_top <- enriched_agg[enriched_agg$log2FoldChange.count>2,]

#order the x axis by classes
enriched_agg_top$Class <- factor(enriched_agg_top$Class, ordered = TRUE,
                                 levels= sort(unique(as.character(enriched_agg_top$Class)),decreasing=FALSE))
enriched_agg_top$Species <- factor(enriched_agg_top$Species, ordered = TRUE,
                                  levels= unique(enriched_agg_top$Species[order(enriched_agg_top$Class)]))
enriched_agg_top$merging <- factor(enriched_agg_top$merging,
                                   levels= c("Dor_winter_spring","Dor_spring_summer","Dor_summer_autumn", "Dor_autumn_winter",  "Dor_winter_summer"))

#plot
PS99_daOTU.p <- ggplot(data=enriched_agg_top,
                       aes(y=log2FoldChange.mean , x=Species, fill = Class, label = log2FoldChange.count))+ 
  #geom_point(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Order, fill = Class, label = log2FoldChange.count), size = 0, shape = 21)+
  geom_text(data=enriched_agg_top,aes(y=log2FoldChange.mean , x=Species), nudge_y= -1.5, nudge_x= 0)+
  geom_errorbar(data=enriched_agg_top,aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.2) +   
  ylab("log2foldchange")+
  #scale_y_reverse()+
  geom_point(data=enriched_agg_top[enriched_agg_top$log2FoldChange.mean<0,], aes(y=log2FoldChange.mean , x=Species, fill = Class, label = log2FoldChange.count), size = 5, shape = 24)+
  geom_point(data=enriched_agg_top[enriched_agg_top$log2FoldChange.mean>0,], aes(y=log2FoldChange.mean , x=Species, fill = Class, label = log2FoldChange.count), size = 5, shape = 25)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  theme_classic(base_size = 10)+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))+
  scale_fill_manual(values = tol21rainbow)+
  guides(shape = 22)+
  coord_flip()+
  facet_grid(.~merging)

PS99_daOTU.p

ggsave("./figures/enrichment_Speciaes.pdf", 
       plot = PS99_daOTU.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


