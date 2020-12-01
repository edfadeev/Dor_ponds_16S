#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")


#load colours
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#####################################
#ASVs enrichment test between the ponds
#####################################
Dor_ps.prev_dry<- subset_samples(Dor_ps.prev_run1, Mic.Season =="Dry")
Dor_ps.prev_dry<- prune_taxa(taxa_sums(Dor_ps.prev_dry)>0,Dor_ps.prev_dry)
sample_data(Dor_ps.prev_dry)$Type<- factor(sample_data(Dor_ps.prev_dry)$Type,
                                            levels = c("Reservoir","Fishpond"))
#run DEseq2
Dor.ddsMat <- phyloseq_to_deseq2(Dor_ps.prev_dry, ~Type)
geoMeans = apply(counts(Dor.ddsMat), 1, gm_mean)
Dor.ddsMat = estimateSizeFactors(Dor.ddsMat, geoMeans = geoMeans)
Dor.ddsMat <- estimateDispersions(Dor.ddsMat)
Dor.DEseq <- DESeq(Dor.ddsMat, fitType="local")
Dor.DEseq.res <- results(Dor.DEseq)

Dor.DEseq.res <- cbind(as(Dor.DEseq.res, "data.frame"),
                       as(tax_table(Dor_ps.prev_dry)[rownames(Dor.DEseq.res), ], "matrix"))
Dor.DEseq.res$ASV <- rownames(Dor.DEseq.res)

#extract only significant ASVs
Dor.DEseq.res.sig <- Dor.DEseq.res %>%
  filter(padj< 0.1,
         abs(log2FoldChange)>1)

#####################################
#Explore the enriched ASVs between the seasons
#####################################
#separate the results by enrichment towards surface and to depth
enriched_Res <- Dor.DEseq.res.sig[Dor.DEseq.res.sig[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","Species","ASV") ]
enriched_Fish  <- Dor.DEseq.res.sig[Dor.DEseq.res.sig[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","Species","ASV")]

#Aggregate on Genus level 
enriched_Dry.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family+Genus+Species,
                                                   enriched_Res, 
                                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

enriched_Wet.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family+Genus+Species,
                                                   enriched_Fish, 
                                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#merge data for ploting
enriched_agg <- rbind(enriched_Dry.agg,enriched_Wet.agg)

#exclude taxa with less than 3 ASV
enriched_agg_top <- enriched_agg[enriched_agg$log2FoldChange.count>2,]

#order the x axis by classes
enriched_agg_top$Class <- factor(enriched_agg_top$Class, ordered = TRUE,
                                 levels= sort(unique(as.character(enriched_agg_top$Class)),decreasing=FALSE))
enriched_agg_top$Genus <- factor(enriched_agg_top$Genus, ordered = TRUE,
                                 levels= unique(enriched_agg_top$Genus[order(enriched_agg_top$Class)]))

#plot
Dor_daASV.p <- ggplot(data=enriched_agg_top,
                      aes(y=log2FoldChange.mean , x=Genus, fill = Class, label = log2FoldChange.count))+ 
  geom_text(data=enriched_agg_top,aes(y=log2FoldChange.mean , x=Genus), size = 8, nudge_y= 0, nudge_x= -0.3)+
  geom_errorbar(data=enriched_agg_top,aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.2) +   
  ylab("log2foldchange")+
  #scale_y_reverse()+
  geom_point(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Genus, label = log2FoldChange.count), size = 9, shape = 21, fill = "black")+
  geom_point(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Genus, fill = Class, label = log2FoldChange.count), size = 7, shape = 21)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_manual(values = class_col)+
  guides(shape = 22)+
  #facet_grid(.~merging)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "none", 
        axis.title.x = element_blank())


#export list of enriched ASVs
write.csv(rbind(enriched_Dry,enriched_Wet), "tables/enriched_ASVs.csv")


#Aggregate on Family level 
enriched_Dry.Fam <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family,
                                                   enriched_Dry, 
                                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
enriched_Wet.Fam <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family,
                                                   enriched_Wet, 
                                                   FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#####################################
#Estimate abundance of the enriched ASVs
#####################################
#Plot barplots of communities
#calculate proportions
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev_run1, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)
Dor_ps.ra.long$Abundance <- Dor_ps.ra.long$Abundance*100

#fix unclassified lineages 
Dor_ps.ra.long$Class <- as.character(Dor_ps.ra.long$Class)
Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Class)] <- paste(Dor_ps.ra.long$Phylum[is.na(Dor_ps.ra.long$Class)],"uc", sep = "_")

#calculate abundance for each Class
Dor_Dry.agg <- Dor_ps.ra.long %>% 
  select(location,Mic.Season,Month,Year,OTU,Class,Abundance)%>%
  filter(Abundance>0, 
         Mic.Season == "Dry",
         OTU %in% enriched_Dry$ASV) %>% 
  group_by(Mic.Season, location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

Dor_Wet.agg <- Dor_ps.ra.long %>% 
  select(location,Mic.Season,Month,Year,OTU,Class,Abundance)%>%
  filter(Abundance>0, 
         Mic.Season == "Wet",
         OTU %in% enriched_Wet$ASV) %>% 
  group_by(Mic.Season, location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

Dor_ps.class.agg <- rbind(Dor_Dry.agg,Dor_Wet.agg) %>% 
  select(location, Mic.Season, Year, Month,Class, Abund.total) %>% 
  melt(id=c("Mic.Season","location", "Year", "Month","Class"), measure.vars = "Abund.total") %>% 
  mutate(#Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
    #                                "May","Jun","Jul","Aug",
    #                                "Sep","Oct","Nov","Dec")),
    Year = factor(Year, levels = c("2013","2014","2015")),
    location = factor(location,levels=c("Res.","V2.","D1.")))


Enriched_bar.p <- ggplot(Dor_ps.class.agg, aes(x = interaction(Month,Year), y = value, fill = Class)) + 
  facet_grid(scale= "free_x", cols = vars(Mic.Season), rows =vars(location)) +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

ggarrange(Dor_daASV.p, Enriched_bar.p, widths = c(1,1),
          ncol = 2, nrow = 1, align = "h", legend = "bottom")

ggsave("./figures/Dor_enrichment_bars.png", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#explore abundance of ASVs
Dor_Dry.genus <- Dor_ps.ra.long %>% 
  select(location,Mic.Season,Month,Year,OTU,Class,Genus,Abundance)%>%
  filter(Abundance>0, 
         Mic.Season == "Dry",
         OTU %in% enriched_Dry$ASV) %>% 
  group_by(Mic.Season, location, Year,Month,Class,Genus) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

Dor_Wet.genus <- Dor_ps.ra.long %>% 
  select(location,Mic.Season,Month,Year,OTU,Class,Genus,Abundance)%>%
  filter(Abundance>0, 
         Mic.Season == "Wet",
         OTU %in% enriched_Wet$ASV) %>% 
  group_by(Mic.Season, location, Year,Month,Class,Genus) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 




#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))