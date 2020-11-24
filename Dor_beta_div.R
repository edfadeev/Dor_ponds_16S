#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(tidyr); packageVersion("tidyr")
library(ggplot2); packageVersion("ggplot2")

#load colours
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#generate seasonal dynamics combined with bar plots
#####################################
Parameters_long<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Season, Year, Month, Temp_degC,
         #Food_Kg_pond, Fish_biomass_g_pond, O2_mg_L, pH, Ammonia_ug_L,
         NO3_NO2_N_L, TP_ug_L, MC_ug_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = "Parameters",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         location = factor(location,levels=c("Res.","V2.","D1.")))


pigments_long<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Season, Year, Month,#Chl_a_mg_L, Chl_b_mg_L,
         Diatoxanthin_mg_L,Dinoxanthin_mg_L,Fucoxanthin_mg_L,
         b_caroten_mg_L,Lutein_mg_L,Zeaxanthin_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = factor(variable),
         variable ="Pigments",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         location = factor(location,levels=c("Res.","V2.","D1.")))

chl_long <- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Season, Year, Month,Chl_a_mg_L, Chl_b_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = factor(variable),
         variable ="Chlorophyll",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         location = factor(location,levels=c("Res.","V2.","D1.")))

#merge together
Parameters_merged_df<- bind_rows(Parameters_long, pigments_long, chl_long) %>% 
  mutate(variable = factor(variable, levels = c("Temp_degC","NO3_NO2_N_L","Ammonia_ug_L","TP_ug_L",
                                                "Pigments","Chlorophyll","MC_ug_L"),
                           labels = c("Temp_degC","NO3_NO2_N_L","Ammonia_ug_L","TP_ug_L",
                                      "Pigments","Chlorophyll","MC_ug_L")),
         Family = factor(Family, levels = c(levels(pigments_long$Family),levels(chl_long$Family),"Parameters"),
                         labels = c(levels(pigments_long$Family),levels(chl_long$Family),"Parameters")))


#Plot barplots of communities
#calculate proportions
Dor_ps.ra <- transform_sample_counts(Dor_ps.prev, function(x) x / sum(x))

#melt phyloseq object
Dor_ps.ra.long <- psmelt(Dor_ps.ra)
Dor_ps.ra.long$Abundance <- Dor_ps.ra.long$Abundance*100

#fix unclassified lineages 
Dor_ps.ra.long$Class <- as.character(Dor_ps.ra.long$Class)
Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Class)] <- paste(Dor_ps.ra.long$Phylum[is.na(Dor_ps.ra.long$Class)],"uc", sep = "_")

Dor_ps.ra.long$Order <- as.character(Dor_ps.ra.long$Order)
Dor_ps.ra.long$Order[is.na(Dor_ps.ra.long$Order)] <- paste(Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Order)],"uc", sep = "_")

Dor_ps.ra.long$Species <- as.character(Dor_ps.ra.long$Species)
Dor_ps.ra.long$Species[is.na(Dor_ps.ra.long$Species)] <- paste(Dor_ps.ra.long$Class[is.na(Dor_ps.ra.long$Species)],"uc", sep = "_")

#calculate abundance for each Class
Dor_ps.ra.long.agg <- Dor_ps.ra.long %>% 
  select(location,Month,Year,OTU,Class,Abundance)%>%
  group_by(location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

#remove below 2% ra
taxa_classes <- unique(Dor_ps.ra.long.agg$Class)
Dor_ps.ra.long.agg$Class[Dor_ps.ra.long.agg$Abund.total<2] <- "Other taxa"
Dor_ps.ra.long.agg$Class <- factor(Dor_ps.ra.long.agg$Class,
                                   levels=c(taxa_classes,"Other taxa"))
Dor_ps.ra.long.agg$Class<- droplevels(Dor_ps.ra.long.agg$Class)

Dor_ps.class.agg <- Dor_ps.ra.long.agg %>% 
  select(location, Year, Month,Class, Abund.total) %>% 
  melt(id=c("location", "Year", "Month","Class"), measure.vars = "Abund.total") %>% 
  mutate(Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         Family = factor(Class),
         variable= "Bacteria",
         location = factor(location,levels=c("Res.","V2.","D1.")))


#####################################
#Generate merged plot for the reservoir
#####################################
#subset reservoir
Res_par_merged <- Parameters_merged_df %>% filter(location =="Res.")
Res_bar_class<- Dor_ps.class.agg %>% filter(location =="Res.")

#plot
Res_par.p<- ggplot()+
  geom_line(data = filter(Res_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(Res_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(Res_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(Res_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(Res_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = variable),size = 1)+
  geom_point(data = filter(Res_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = variable), size = 2)+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  guides(colour=guide_legend(title="Pigments"))+
  facet_grid(variable~Year, scales = "free_y")+
  theme_bw()+
  theme(legend.position="bottom", axis.text.x = element_blank(),axis.title.x = element_blank())

Res_bar.p <- ggplot(Res_bar_class, aes(x = Month, y = value, fill = Class)) + 
  facet_grid(variable~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom",strip.text.x = element_blank())


ggarrange(Res_par.p, Res_bar.p, heights = c(2,1.2),
          ncol = 1, nrow = 2, align = "v", legend = "bottom",
          legend.grob = do.call(rbind, c(list(get_legend(Res_bar.p),get_legend(Res_par.p)), size="first")))

ggsave("./figures/Res_overview.png", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#Generate merged plot for D1 and V2
#####################################
#subset reservoir
Pond_par_merged <- Parameters_merged_df %>% filter(location  %in% c("D1.","V2."))
Pond_bar_class<- Dor_ps.class.agg %>% filter(location %in% c("D1.","V2."))

#plot
Pond_par.p<- ggplot()+
  geom_line(data = filter(Pond_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location)),size = 1)+
  geom_point(data = filter(Pond_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), shape = location),size = 3)+
  geom_line(data = filter(Pond_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location)),size = 1)+
  geom_point(data = filter(Pond_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), shape = location),size = 3)+
  geom_line(data = filter(Pond_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = interaction(location,variable)),size = 1)+
  geom_point(data = filter(Pond_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = interaction(location,variable), shape = location), size = 3)+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  guides(colour=guide_legend(title="Pigments"), shape=guide_legend(title="Pool") )+
  facet_grid(variable~Year, scales = "free_y")+
  theme_bw()+
  theme(legend.position="none", axis.text.x = element_blank(),axis.title.x = element_blank())

Pond_bar.p <- ggplot(Pond_bar_class, aes(x = Month, y = value, fill = Class)) + 
  facet_grid(location~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom",strip.text.x = element_blank())


ggarrange(Pond_par.p, Pond_bar.p, heights = c(2,1.2),
          ncol = 1, nrow = 2, align = "v")

ggsave("./figures/Fishponds_overview.png", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#Community dissimilarities
#####################################
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#transform counts using geometric mean
Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev_run1)
#Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev)

#NMDS plot
Dor_ps.gm_mean.ord <- ordinate(Dor_ps.gm_mean, method = "NMDS", distance = "euclidean")
Dor.ord.df <- plot_ordination(Dor_ps.gm_mean, Dor_ps.gm_mean.ord, axes = c(1,2,3),justDF = TRUE)

Dor.ord.p <- ggplot()+
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = NMDS2, shape = Year), 
             fill = "black", size = 5) +
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = NMDS2, colour = Mic.Season, shape = Year), 
             size = 4) +
  geom_text(data = Dor.ord.df,aes(x = NMDS1, y = NMDS2,label = location), 
            nudge_y= -8,size=3)+
  scale_colour_manual(values = c("Wet"="darkblue",
                                 "Dry"="orange")) + 
  annotate(geom="text", x=-110, y=110, label= paste0("Stress = ", round(Dor_ps.gm_mean.ord$stress,2)),
           color="red", size = 5)+
  theme_bw()+
  theme(legend.position = "bottom")


ggsave("./figures/NMDS_2013_2014.png", 
       plot = Dor.ord.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#test grouping of samples
df <- as(sample_data(Dor_ps.gm_mean), "data.frame")
d <- phyloseq::distance(Dor_ps.gm_mean, "euclidean")
adonis_all <- adonis2(d ~ Year + Mic.Season + location , df)
adonis_all

#posthoc to check which seasons are different
groups <- df[["Mic.Season"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)

#posthoc to check which seasons are different
groups <- df[["location"]]
mod <- betadisper(d, groups)
permutest(mod)

#dispersion is different between groups
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)


#####################################
#ASVs enrichment test
#####################################
#split dataset
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
#Explore the enriched ASVs between winter and summer
#####################################
deseq_res_all<- deseq_res_all %>% filter(merging =="Dor_winter_summer")

#separate the results by enrichment towards surface and to depth
enriched_shallow <- deseq_res_all[deseq_res_all[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV","merging") ]
enriched_deep  <- deseq_res_all[deseq_res_all[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV","merging")]

#Aggregate on Genus level 
enriched_shallow.agg <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family+Genus,
                                                       enriched_shallow, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
enriched_deep.agg <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family+Genus,
                                                    enriched_deep, 
                                                    FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#merge data for ploting
enriched_agg <- rbind(enriched_shallow.agg,enriched_deep.agg)

#exclude taxa with less than 3 ASV
enriched_agg_top <- enriched_agg[enriched_agg$log2FoldChange.count>2,]

#order the x axis by classes
enriched_agg_top$Class <- factor(enriched_agg_top$Class, ordered = TRUE,
                                 levels= sort(unique(as.character(enriched_agg_top$Class)),decreasing=FALSE))
enriched_agg_top$Genus <- factor(enriched_agg_top$Genus, ordered = TRUE,
                                 levels= unique(enriched_agg_top$Genus[order(enriched_agg_top$Class)]))
enriched_agg_top$merging <- factor(enriched_agg_top$merging,
                                   levels= c("Dor_winter_spring","Dor_spring_summer","Dor_summer_autumn", "Dor_autumn_winter",  "Dor_winter_summer"))




enriched_agg_top<- enriched_agg_top[enriched_agg_top$merging %in% c("Dor_winter_summer"),]


#plot
Dor_daOTU.p <- ggplot(data=enriched_agg_top,
                       aes(y=log2FoldChange.mean , x=Genus, fill = Class, label = log2FoldChange.count))+ 
  geom_text(data=enriched_agg_top,aes(y=log2FoldChange.mean , x=Genus), nudge_y= 0, nudge_x= -0.3)+
  geom_errorbar(data=enriched_agg_top,aes(ymin = log2FoldChange.mean-log2FoldChange.se, ymax = log2FoldChange.mean +log2FoldChange.se), width = 0.2) +   
  ylab("log2foldchange")+
  #scale_y_reverse()+
  geom_point(data=enriched_agg_top, aes(y=log2FoldChange.mean , x=Genus, fill = Class, label = log2FoldChange.count), size = 5, shape = 21)+
  #geom_point(data=enriched_agg_top[enriched_agg_top$log2FoldChange.mean>0,], aes(y=log2FoldChange.mean , x=Genus, fill = Class, label = log2FoldChange.count), size = 5, shape = 25)+
  geom_hline(aes(yintercept=0), linetype="dashed")+
  scale_fill_manual(values = class_col)+
  guides(shape = 22)+
  facet_grid(.~merging)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle =90))


ggsave("./figures/Dor_enrichment.png", 
       plot = Dor_daOTU.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



#Aggregate on Family level 
enriched_winter.Fam <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family,
                                                      enriched_shallow, 
                                                      FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))
enriched_summer.Fam <-as.data.frame(as.list(aggregate(log2FoldChange~merging+Phylum+Class+Order+Family,
                                                      enriched_deep, 
                                                      FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

#####################################
#get session info and remove all objects and libraries
#####################################
sessionInfo()

rm(list = ls(all = TRUE))