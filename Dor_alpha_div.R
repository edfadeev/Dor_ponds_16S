#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("cowplot"); packageVersion("cowplot")
library("reshape2"); packageVersion("reshape2")
library("rstatix"); packageVersion("rstatix")
library("ggpubr"); packageVersion("ggpubr")
library("venn"); packageVersion("venn")

#load colors
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Alpha diversity statistical tests
####################################
# Calculate richness
Dor_ps.alpha.div <- estimate_richness(Dor_ps.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
Dor_comm.char<- data.frame(Sample_number_dada2 = sample_data(Dor_ps.prev)$Sample_number_dada2,
                             Pool = gsub("\\.","",sample_data(Dor_ps.prev)$location),
                           Year = sample_data(Dor_ps.prev)$Year,
                             Month = sample_data(Dor_ps.prev)$Month,
                             Mic.Season = sample_data(Dor_ps.prev)$Mic.Season,
                             #Chl.sequences = sample_sums(Dor_ps.chl),
                             Sequences= sample_sums(Dor_ps.prev),
                             Observed = Dor_ps.alpha.div$Observed,
                             Chao1 = Dor_ps.alpha.div$Chao1,
                             Completness = round(100*Dor_ps.alpha.div$Observed/Dor_ps.alpha.div$Chao1, digits=2),
                             Shannon = round(Dor_ps.alpha.div$Shannon,digits=2),
                             InvSimpson = round(Dor_ps.alpha.div$InvSimpson,digits=2),
                             Evenness = round(Dor_ps.alpha.div$Shannon/log(Dor_ps.alpha.div$Observed),digits=2))

#merge with env. par.
env.par<- c("Temp_degC",  
            "Ammonia_ug_L", "NO3_NO2_N_L", "TP_ug_L","Food..Kg.","Biomass..kg.","Chl_a_mg_L", "Chl_b_mg_L",
            "Diatoxanthin_mg_L","Dinoxanthin_mg_L","Fucoxanthin_mg_L",
            "b_caroten_mg_L","MC_ug_L","Lutein_mg_L","Zeaxanthin_mg_L")

Dor_metadata<- left_join(Dor_comm.char, sample_data(Dor_ps.prev), by = c("Sample_number_dada2","Year","Month","Mic.Season")) %>% 
  dplyr::select(c(names(Dor_comm.char),env.par)) %>% 
  mutate("N:P" = as.numeric(NO3_NO2_N_L/TP_ug_L)) %>%   arrange(Pool,Year,Month)

write.csv(Dor_metadata, "./tables/Dor_alpha_table.csv")

# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#plot alpha diversity
Dor_alpha <- estimate_richness(Dor_ps.prev_run1, measures = c("Observed", "Chao1","Shannon", "Evenness"))
Dor_alpha$Evenness <- Dor_alpha$Shannon/log(Dor_alpha$Observed)

Dor_alpha <- merge_phyloseq(Dor_ps.prev_run1, sample_data(Dor_alpha))

Dor_alpha.m <- as(sample_data(Dor_alpha), "data.frame")%>%
  select(location, Year, Month, Mic.Season, Observed, Chao1, Shannon, Evenness)%>%
  melt(id.vars = c("location", "Year","Month", "Mic.Season"))

alpha.p<- ggplot(Dor_alpha.m, aes(x = Month, y = value, group = variable)) +
  labs(x = "Year", y = "Alpha diversity")+
  geom_point(aes(shape = location), size =3)+
  geom_smooth(method = loess, se = TRUE)+
  #scale_fill_manual(values =c("yellow","darkgreen"))+
  #geom_boxplot(outlier.color = NULL, notch = FALSE)+
  facet_grid(variable~Year, scales = "free")+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  coord_cartesian(clip="off")+
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("./figures/alpha_p.png", 
       plot = alpha.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

alpha_pool.p<- ggplot(Dor_alpha.m, aes (x = location, y = value, group = location, colour = Year))+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  geom_jitter(size = 3)+
  facet_wrap(variable~., scales = "free", ncol = 2)+
  theme_classic(base_size = 12)+
  #geom_signif(comparisons = list(c("D1.", "Res."),c("D1.","V2."),c("Res.","V2.")),
  #            map_signif_level=TRUE, test = "wilcox.test", color = "black")+
  theme(legend.position = "bottom")

ggsave("./figures/alpha_pools.png", 
       plot = alpha_pool.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

alpha_seasons.p<- ggplot(Dor_alpha.m, aes (x = location, y = value, group = interaction(location,Mic.Season), 
                                           colour = Mic.Season, shape = as.factor(Year)))+
  geom_boxplot(outlier.color = NA, notch = FALSE)+
  geom_jitter(size = 3)+
  facet_wrap(variable~., scales = "free", ncol = 2)+
  scale_colour_manual(values = c("Wet"="darkblue",
                                 "Dry"="orange")) + 
  #coord_fixed()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),legend.position = "bottom")

ggsave("./figures/alpha_seasons.pdf", 
       plot = alpha_seasons.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)



#####################################
#Test statistical differences between Years and Seasons
####################################
shapiro.test(sample_data(Dor_alpha)$Chao1)
#Chao1 richness did not show normal distribution (p < 0.01), thus will be analyzed using Kruskal Wallis test

kruskal.test(Chao1 ~ location, data = data.frame(sample_data(Dor_alpha)))

kruskal.test(Chao1 ~ Mic.Season, data = data.frame(sample_data(Dor_alpha)))

Chao1_Wilcox_Season <- as(sample_data(Dor_alpha),"data.frame")   %>%
  group_by(location) %>% 
  rstatix::wilcox_test(Chao1 ~ Mic.Season, p.adjust.method = "BH") %>%
  add_significance()

Chao1_Wilcox_Season <- as(sample_data(Dor_alpha),"data.frame")   %>%
  #group_by(Mic.Season) %>% 
  rstatix::wilcox_test(Chao1 ~ location, p.adjust.method = "BH") %>%
  add_significance()


#####################################
#ASVs overlap between pools
####################################
#subset each pool
Dor_ps.D1<- subset_samples(Dor_ps.prev_run1, location =="D1.")
Dor_ps.D1<- prune_taxa(taxa_sums(Dor_ps.D1)>0,Dor_ps.D1)


Dor_ps.V2<- subset_samples(Dor_ps.prev_run1, location =="V2.")
Dor_ps.V2<- prune_taxa(taxa_sums(Dor_ps.V2)>0,Dor_ps.V2)

Dor_ps.Res<- subset_samples(Dor_ps.prev_run1, location =="Res.")
Dor_ps.Res<- prune_taxa(taxa_sums(Dor_ps.Res)>0,Dor_ps.Res)

#generate list of ASVs in each pool
z <- list()
z[["D1"]] <- as.character(row.names(otu_table(Dor_ps.D1)))
z[["V2"]] <- as.character(row.names(otu_table(Dor_ps.V2)))
z[["Res"]] <- as.character(row.names(otu_table(Dor_ps.Res)))

#plot
png(file="figures/venn_pools.png",units = "cm", res = 300,
    width=30, height=30)
venn(z, snames = names(z), ilab=TRUE, zcolor = "style",
     ilcs = 1, sncs = 2)
dev.off()
