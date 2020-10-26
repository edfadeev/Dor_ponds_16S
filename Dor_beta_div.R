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
#Community dissimilarities
#####################################
#transforma counts using gemetric mean
Dor_ps.gm_mean <- phyloseq_gm_mean_trans(Dor_ps.prev)

#significance test
df <- as(sample_data(Dor_ps.gm_mean), "data.frame")
d <- phyloseq::distance(Dor_ps.gm_mean, "euclidean")
adonis_all <- adonis2(d ~ Year + Season + location + Season*Year , df)
adonis_all

#posthoc to check which seasons are different
groups <- df[["Season"]]
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
#Plot barplots of communities
#####################################
#transform data
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

#remove below 1% ra
taxa_classes <- unique(Dor_ps.ra.long.agg$Class)
Dor_ps.ra.long.agg$Class[Dor_ps.ra.long.agg$Abund.total<2] <- "Other taxa"
Dor_ps.ra.long.agg$Class <- factor(Dor_ps.ra.long.agg$Class,
                                       levels=c(taxa_classes,"Other taxa"))
Dor_ps.ra.long.agg$Class<- droplevels(Dor_ps.ra.long.agg$Class)

#reorder location
Dor_ps.ra.long.agg$location <- factor(Dor_ps.ra.long.agg$location,
                                       levels=c("Res.","V2.","D1."))

#Plot 
barplots_total <- ggplot(Dor_ps.ra.long.agg, aes(x = Month, y = Abund.total, fill = Class)) + 
  facet_grid(location~Year, space= "fixed") +
  geom_col()+
  scale_fill_manual(values = tol21rainbow) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Sequence proportions (%) \n")+
  theme_bw()+
  theme(legend.position="bottom")


ggsave("./figures/Run1_barplots_total.pdf", 
       plot = barplots_total,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)

