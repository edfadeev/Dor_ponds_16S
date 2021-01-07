#load libraries
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(ggpubr); packageVersion("ggpubr")
#library(tidyr); packageVersion("tidyr")


#load colours
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#generate seasonal dynamics combined with bar plots
#####################################
Parameters_long<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Mic.Season, Year, Month, Temp_degC,
         Food_Kg_pond, Fish_biomass_Kg_pond,
         O2_mg_L, 
         #pH, 
         Ammonia_ug_L,
         NO3_NO2_N_L, TP_ug_L, MC_ug_L) %>% 
  mutate("N:P" = (as.numeric(NO3_NO2_N_L)/as.numeric(TP_ug_L))) %>% 
  melt(id=c("location","Mic.Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = "Parameters",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         location = factor(location,levels=c("Res.","V2.","D1.")))


pigments_long<- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Mic.Season, Year, Month,#Chl_a_mg_L, Chl_b_mg_L,
         Diatoxanthin_mg_L,Dinoxanthin_mg_L,Fucoxanthin_mg_L,
         b_caroten_mg_L,Lutein_mg_L,Zeaxanthin_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Mic.Season", "Year", "Month")) %>% 
  mutate(value =as.numeric(value),
         Family = factor(variable),
         variable ="Pigments",
         Month = factor(Month, levels = c("Jan","Feb","Mar","Apr",
                                          "May","Jun","Jul","Aug",
                                          "Sep","Oct","Nov","Dec")),
         Year = factor(Year, levels = c("2013","2014","2015")),
         location = factor(location,levels=c("Res.","V2.","D1.")))

chl_long <- as(sample_data(Dor_ps.prev),"data.frame") %>% 
  select(location, Mic.Season, Year, Month,Chl_a_mg_L, Chl_b_mg_L) %>% 
  #filter(location =="Res.") %>% 
  melt(id=c("location","Mic.Season", "Year", "Month")) %>% 
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
  mutate(variable = factor(variable, levels = c("Temp_degC","NO3_NO2_N_L","Ammonia_ug_L","TP_ug_L","N:P","Chlorophyll","Pigments","MC_ug_L","O2_mg_L", 
                                                "pH", "Food_Kg_pond", "Fish_biomass_Kg_pond"),
                           labels = c("Temp_degC","NO3_NO2_N_L","Ammonia_ug_L","TP_ug_L","N:P","Chlorophyll","Pigments","MC_ug_L","O2_mg_L", 
                                      "pH", "Food_Kg_pond", "Fish_biomass_Kg_pond")),
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
  select(Mic.Season, location,Month,Year,OTU,Class,Abundance)%>%
  group_by(Mic.Season, location, Year,Month,Class) %>%
  dplyr::summarise(Abund.total= sum(Abundance)) 

#remove below 2% ra
taxa_classes <- unique(Dor_ps.ra.long.agg$Class)
Dor_ps.ra.long.agg$Class[Dor_ps.ra.long.agg$Abund.total<2] <- "Other taxa"
Dor_ps.ra.long.agg$Class <- factor(Dor_ps.ra.long.agg$Class,
                                   levels=c(taxa_classes,"Other taxa"))
Dor_ps.ra.long.agg$Class<- droplevels(Dor_ps.ra.long.agg$Class)

Dor_ps.class.agg <- Dor_ps.ra.long.agg %>% 
  select(Mic.Season, location, Year, Month,Class, Abund.total) %>% 
  melt(id=c("Mic.Season", "location", "Year", "Month","Class"), measure.vars = "Abund.total") %>% 
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
Res_par_merged <- Parameters_merged_df %>% filter(location =="Res.", variable %in% c("Temp_degC","NO3_NO2_N_L","Ammonia_ug_L","TP_ug_L","O2_mg_L",
                                                                                     "N:P","Pigments","Chlorophyll","MC_ug_L"))
Res_bar_class<- Dor_ps.class.agg %>% filter(location =="Res.")

#plot
Res_par.p<- ggplot()+
  geom_line(data = filter(Res_par_merged, variable =="Pigments"), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(Res_par_merged, variable =="Pigments"), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(Res_par_merged, variable =="Chlorophyll"), aes(x= Month, y = value, colour = Family, group = Family),size = 1)+
  geom_point(data = filter(Res_par_merged, variable =="Chlorophyll"), aes(x= Month, y = value, colour = Family, group = Family),size = 2)+
  geom_line(data = filter(Res_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = variable),size = 1)+
  geom_point(data = filter(Res_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = variable), size = 2)+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  guides(colour=guide_legend(title="Pigments"))+
  facet_grid(variable~Year, scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),axis.text.x = element_blank(),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

Res_bar.p <- ggplot(Res_bar_class, aes(x = Month, y = value, fill = Class)) + 
  facet_grid(variable~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


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
Pond_par_merged <- Parameters_merged_df %>% filter(location  %in% c("D1.","V2."), variable %in% c("Temp_degC","NO3_NO2_N_L","TP_ug_L",#"Ammonia_ug_L","N:P","MC_ug_L","O2_mg_L",
                                                                                                 "Pigments","Chlorophyll","Food_Kg_pond", "Fish_biomass_Kg_pond"))
Pond_bar_class<- Dor_ps.class.agg %>% filter(location %in% c("D1.","V2."))

#plot
Pond_par.p<- ggplot()+
  geom_line(data = filter(Pond_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), linetype = location),size = 1)+
  geom_point(data = filter(Pond_par_merged, variable %in% c("Pigments")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), shape = location),size = 3)+
  geom_line(data = filter(Pond_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), linetype = location),size = 1)+
  geom_point(data = filter(Pond_par_merged, variable %in% c("Chlorophyll")), aes(x= Month, y = value, colour = Family, group = interaction(Family,location), shape = location),size = 3)+
  geom_line(data = filter(Pond_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = interaction(location,variable), linetype = location),size = 1)+
  geom_point(data = filter(Pond_par_merged, Family =="Parameters"), aes(x= Month, y = value, group = interaction(location,variable), shape = location), size = 3)+
  scale_fill_manual(values = cbPalette)+
  scale_colour_manual(values = cbPalette)+
  guides(colour=guide_legend(title="Pigments"), shape=guide_legend(title="Pool") )+
  facet_grid(variable~Year, scales = "free_y")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())

Pond_bar.p <- ggplot(Pond_bar_class, aes(x = Month, y = value, fill = Class)) + 
  facet_grid(location~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = class_col) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom", 
        axis.title.x = element_blank())


ggarrange(Pond_par.p, Pond_bar.p, heights = c(2,1.2),
          ncol = 1, nrow = 2, align = "v")

ggsave("./figures/Fishponds_overview.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#####################################
#Community dissimilarities 2013-2014
#####################################
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#transform counts using geometric mean
Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev_run1)
#Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev)

#NMDS plot
Dor_ps.gm_mean.ord <- ordinate(Dor_ps.gm_mean, method = "NMDS", distance = "euclidean")
Dor.ord.df <- plot_ordination(Dor_ps.gm_mean, Dor_ps.gm_mean.ord, axes = c(1,2,3),justDF = TRUE)

Dor.ord.df <- Dor.ord.df %>%  mutate(Temp_degC = as.numeric(Temp_degC),
                                     Mic.Season = factor(Mic.Season, levels =c("Dry","Wet")))

#add centroids
Dor.ord.df <- merge(Dor.ord.df,aggregate(cbind(mean.x=NMDS1,mean.y=NMDS2)~location,Dor.ord.df,mean),by="location")


Dor.ord.p <- ggplot(data = Dor.ord.df)+
  geom_point(aes(x = NMDS1, y = NMDS2, shape = location), 
             fill = "black", size = 5,alpha = 0.8) +
  #geom_point(data = Dor.ord.df, aes(x = NMDS1, y = NMDS2, colour = Mic.Season, shape = Year), 
  #           size = 7) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = Temp_degC, shape = location), 
             size = 3,alpha = 0.8) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = Comment), 
            nudge_y= -8,size=5)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  annotate(geom="text", x=-100, y=100, label= paste0("Stress = ", round(Dor_ps.gm_mean.ord$stress,2)),
           color="red", size = 5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


#plot NMDS1 vs Temperature
Dor.ord.temp.p <- ggplot()+
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = Temp_degC, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = Temp_degC, colour = Temp_degC, shape = location), 
             size = 3) +
  geom_text(data = Dor.ord.df,aes(x = NMDS1, y = Temp_degC,label = substr(Mic.Season, 1, 1)), 
            nudge_y= -0.8,size=5)+
  geom_hline(yintercept = 20, linetype = 2)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#plot NMDS2 vs Temperature
Dor.ord.temp1.p <- ggplot()+
  geom_point(data = Dor.ord.df, aes(x = NMDS2, y = Temp_degC, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor.ord.df, aes(x = NMDS2, y = Temp_degC, colour = Temp_degC, shape = location), 
             size = 3) +
  geom_text(data = Dor.ord.df,aes(x = NMDS2, y = Temp_degC, label = substr(Mic.Season, 1, 1)), 
            nudge_y= -0.8,size=5)+
  geom_hline(yintercept = 20, linetype = 2)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


ggarrange(Dor.ord.p, Dor.ord.temp.p, Dor.ord.temp1.p, widths = 3,heights = 3,labels="AUTO",
          #ncol = 3, nrow = 1, 
          align = "hv", legend = "bottom",common.legend = TRUE)


ggsave("./figures/NMDS_2013_2014.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#test grouping of samples
df <- as(sample_data(Dor_ps.gm_mean), "data.frame")
d <- phyloseq::distance(Dor_ps.gm_mean, "euclidean")
adonis_all <- adonis2(d ~ Year +Mic.Season + location  , df)
adonis_all

#posthoc to check which ponds are different
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
#Community dissimilarities 2015
#####################################
# subset only 2015
Dor_ps.prev_run2<- subset_samples(Dor_ps.prev, Run == "2")

#transform counts using geometric mean
Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev_run2)

#NMDS plot
Dor_ps.gm_mean.ord <- ordinate(Dor_ps.gm_mean, method = "NMDS", distance = "euclidean")
Dor.ord.df <- plot_ordination(Dor_ps.gm_mean, Dor_ps.gm_mean.ord, axes = c(1,2,3),justDF = TRUE)

Dor.ord.df <- Dor.ord.df %>%  mutate(Temp_degC = as.numeric(Temp_degC),
                                     Mic.Season = factor(Mic.Season, levels =c("Dry","Wet")))

#add centroids
Dor.ord.df <- merge(Dor.ord.df,aggregate(cbind(mean.x=NMDS1,mean.y=NMDS2)~location,Dor.ord.df,mean),by="location")


Dor.ord.p <- ggplot(data = Dor.ord.df)+
  geom_point(aes(x = NMDS1, y = NMDS2, shape = location), 
             fill = "black", size = 5,alpha = 0.8) +
  #geom_point(data = Dor.ord.df, aes(x = NMDS1, y = NMDS2, colour = Mic.Season, shape = Year), 
  #           size = 7) +
  geom_point(aes(x = NMDS1, y = NMDS2, colour = Temp_degC, shape = location), 
             size = 3,alpha = 0.8) +
  geom_text(aes(x = NMDS1, y = NMDS2,label = Comment), 
            nudge_y= -8,size=5)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  annotate(geom="text", x=-100, y=100, label= paste0("Stress = ", round(Dor_ps.gm_mean.ord$stress,2)),
           color="red", size = 5)+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


#plot NMDS1 vs Temperature
Dor.ord.temp.p <- ggplot()+
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = Temp_degC, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor.ord.df, aes(x = NMDS1, y = Temp_degC, colour = Temp_degC, shape = location), 
             size = 3) +
  geom_text(data = Dor.ord.df,aes(x = NMDS1, y = Temp_degC,label = substr(Mic.Season, 1, 1)), 
            nudge_y= -0.8,size=5)+
  geom_hline(yintercept = 20, linetype = 2)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

#plot NMDS2 vs Temperature
Dor.ord.temp1.p <- ggplot()+
  geom_point(data = Dor.ord.df, aes(x = NMDS2, y = Temp_degC, shape = location), 
             fill = "black", size = 5) +
  geom_point(data = Dor.ord.df, aes(x = NMDS2, y = Temp_degC, colour = Temp_degC, shape = location), 
             size = 3) +
  geom_text(data = Dor.ord.df,aes(x = NMDS2, y = Temp_degC, label = substr(Mic.Season, 1, 1)), 
            nudge_y= -0.8,size=5)+
  geom_hline(yintercept = 20, linetype = 2)+
  scale_colour_gradient(low = "blue", high = "yellow")+
  coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")


ggarrange(Dor.ord.p, Dor.ord.temp.p, Dor.ord.temp1.p, widths = 3,heights = 3,labels="AUTO",
          #ncol = 3, nrow = 1, 
          align = "hv", legend = "bottom",common.legend = TRUE)


ggsave("./figures/NMDS_2015.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

#test grouping of samples
df <- as(sample_data(Dor_ps.gm_mean), "data.frame")
d <- phyloseq::distance(Dor_ps.gm_mean, "euclidean")
adonis_all <- adonis2(d ~ Mic.Season*location  , df)
adonis_all

#####################################
#ASVs enrichment test between the seasons
#####################################
sample_data(Dor_ps.prev_run1)$Mic.Season<- factor(sample_data(Dor_ps.prev_run1)$Mic.Season,
                                                   levels = c("Dry","Wet"))
#run DEseq2
Dor.ddsMat <- phyloseq_to_deseq2(Dor_ps.prev_run1, ~Mic.Season)
geoMeans = apply(counts(Dor.ddsMat), 1, gm_mean)
Dor.ddsMat = estimateSizeFactors(Dor.ddsMat, geoMeans = geoMeans)
Dor.ddsMat <- estimateDispersions(Dor.ddsMat)
Dor.DEseq <- DESeq(Dor.ddsMat, fitType="local")
Dor.DEseq.res <- results(Dor.DEseq)

Dor.DEseq.res <- cbind(as(Dor.DEseq.res, "data.frame"),
                           as(tax_table(Dor_ps.prev_run1)[rownames(Dor.DEseq.res), ], "matrix"))
Dor.DEseq.res$ASV <- rownames(Dor.DEseq.res)

#extract only significant ASVs
Dor.DEseq.res.sig <- Dor.DEseq.res %>%
                      filter(padj< 0.1,
                             abs(log2FoldChange)>1)

#####################################
#Explore the enriched ASVs between the seasons
#####################################
#separate the results by enrichment towards surface and to depth
enriched_Dry <- Dor.DEseq.res.sig[Dor.DEseq.res.sig[, "log2FoldChange"] < 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV") ]
enriched_Wet  <- Dor.DEseq.res.sig[Dor.DEseq.res.sig[, "log2FoldChange"] > 0,c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order", "Family", "Genus","ASV")]

#Aggregate on Genus level 
enriched_Dry.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family+Genus,
                                                   enriched_Dry, 
                                                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x)), count=length(x)))))

enriched_Wet.agg <-as.data.frame(as.list(aggregate(log2FoldChange~Phylum+Class+Order+Family+Genus,
                                                   enriched_Wet, 
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