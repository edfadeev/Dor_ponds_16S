#load libraries
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("cowplot"); packageVersion("cowplot")
library("reshape2"); packageVersion("reshape2")

#load colour palettes
source('./Scripts/color_palettes.R')

se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
#####################################
#Alpha diversity statistical tests
####################################
# Calculate richness
Dor_ps.prev_alpha.div <- estimate_richness(Dor_ps.prev, split = TRUE, measures = NULL)

#generate data set with all bacterial community characteristics
Dor_comm.char<- data.frame(  Sample = sample_names(Dor_ps.prev),
                             Year = sample_data(Dor_ps.prev)$Year,
                             Month = sample_data(Dor_ps.prev)$Month,
                             Season = sample_data(Dor_ps.prev)$Season,
                             Sequences= sample_sums(Dor_ps.prev),
                             Observed = Dor_ps.prev_alpha.div$Observed,
                             Chao1 = Dor_ps.prev_alpha.div$Chao1,
                             Completness = round(100*Dor_ps.prev_alpha.div$Observed/Dor_ps.prev_alpha.div$Chao1, digits=2),
                             Shannon = round(Dor_ps.prev_alpha.div$Shannon,digits=2),
                             InvSimpson = round(Dor_ps.prev_alpha.div$InvSimpson,digits=2),
                             Evenness = round(Dor_ps.prev_alpha.div$Shannon/log(Dor_ps.prev_alpha.div$Observed),digits=2))

write.csv(Dor_comm.char, "./alpha_table_prev.csv")

#plot alpha diversity
Dor_alpha <- estimate_richness(Dor_ps.prev, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
Dor_alpha <- merge_phyloseq(Dor_ps.prev, sample_data(Dor_alpha))

Dor_alpha.m <- as(sample_data(Dor_alpha), "data.frame")%>%
  select(Year, Month, Season, Observed, Chao1, Shannon, InvSimpson)%>%
  melt(id.vars = c("Year","Month", "Season"))


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



alpha_seasons.p<- ggplot(Dor_alpha.m, aes (x = Season, y = value, group = Season, colour = Year))+
  geom_boxplot(outlier.color = NULL, notch = FALSE)+
  geom_jitter(size = 3)+
  facet_wrap(variable~., scales = "free", ncol = 2)+
  theme_classic(base_size = 12)+
  theme(legend.position = "bottom")

ggsave("./figures/alpha_seasons.pdf", 
       plot = alpha_seasons.p,
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


#Chao1 summary
Dor_alpha.Chao1.agg <- do.call(data.frame, aggregate(Chao1~ Year+Season, as(sample_data(Dor_alpha), "data.frame"), function(x) c(mean = mean(x), se = se(x), median = median(x))))

#Shanonn summary
Dor_alpha.Shanonn.agg <- do.call(data.frame, aggregate(Shanonn~ Year+Season, as(sample_data(Dor_alpha), "data.frame"), function(x) c(mean = mean(x), se = se(x),median = median(x))))

#InvSimpson summary
Dor_alpha.InvSimpson.agg <- do.call(data.frame, aggregate(InvSimpson~ Year+Season, as(sample_data(Dor_alpha), "data.frame"), function(x) c(mean = mean(x), se = se(x),median = median(x))))

