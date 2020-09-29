


#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("olsrr"); packageVersion("olsrr")
library("cowplot"); packageVersion("cowplot")

#set plots theme
theme_set(theme_classic())

#load colour palettes
source('./Scripts/color_palettes.R')


#####################################
#Alpha diversity statistical tests
####################################
# Calculate richness
Dor_ps.prev_alpha.div <- estimate_richness(Dor_ps.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
Dor_comm.char<- data.frame(  Sample = sample_names(Dor_ps.prev),
                             Year = sample_data(Dor_ps.prev)$Year,
                             Month = sample_data(Dor_ps.prev)$Month,
                             Sequences= sample_sums(Dor_ps.prev),
                             Observed = Dor_ps.prev_alpha.div$Observed,
                             Chao1 = Dor_ps.prev_alpha.div$Chao1,
                             Completness = round(100*Dor_ps.prev_alpha.div$Observed/Dor_ps.prev_alpha.div$Chao1, digits=2),
                             Shannon = round(Dor_ps.prev_alpha.div$Shannon,digits=2),
                             InvSimpson = round(Dor_ps.prev_alpha.div$InvSimpson,digits=2),
                             Evenness = round(Dor_ps.prev_alpha.div$Shannon/log(Dor_ps.prev_alpha.div$Observed),digits=2))

write.csv(Dor_comm.char, "./alpha_table_prev.csv")

# Create new ps object with diversity estimates added to sample_data
Dor_ps.prev_div <- merge_phyloseq(Dor_ps.prev, sample_data(Dor_comm.char))

#plot alpha diversity
Dor_alpha <- estimate_richness(Dor_ps.prev, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
Dor_alpha <- merge_phyloseq(Dor_ps.prev, sample_data(Dor_alpha))

Dor_alpha.m <- as(sample_data(Dor_alpha), "data.frame")%>%
  select(Year, Month, Observed, Chao1, Shannon, InvSimpson)%>%
  melt(id.vars = c("Year","Month"))


alpha.p<- ggplot(Dor_alpha.m, aes(x = Month, y = value, group = variable)) +
  labs(x = "Year", y = "Alpha diversity")+
  geom_point(size =3)+
  geom_smooth(method = loess, se = TRUE)+
  #scale_fill_manual(values =c("yellow","darkgreen"))+
  #geom_boxplot(outlier.color = NULL, notch = FALSE)+
  facet_wrap(variable~Year, scales = "free", ncol = 3)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  geom_vline(aes(xintercept=Inf))+
  coord_cartesian(clip="off")+
  theme_classic() +
  theme(legend.position = "bottom")

ggsave("./figures/alpha_p.pdf", 
       plot = alpha.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)
