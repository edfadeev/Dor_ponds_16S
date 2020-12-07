#load colors and functions
source("scripts/color_palletes.R")
source("scripts/extra_functions.R")

#load RDS object
Dor_ps.prev <-readRDS("data/Dor_ps_prev.rds")

#####################################
#Subset samples that have metadata
####################################
# subset only 2013-2014
Dor_ps.prev_run1<- subset_samples(Dor_ps.prev, Run == "1")

#transform to geometric mean and remove 
Dor_ps.prev_gm<- phyloseq_gm_mean_trans(Dor_ps.prev_run1)

#####################################
#plot distances distribution between fractions
#####################################
sample_data(Dor_ps.prev_gm)$SampleID <- sample_names(Dor_ps.prev_gm)
frac_distances_all <- data.frame()

frac_distances = as.matrix(phyloseq::distance(Dor_ps.prev_gm, "euclidean"))
frac_distances.m = reshape2::melt(frac_distances)

# remove self-comparisons
frac_distances.m = frac_distances.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
samples_data = as(sample_data(Dor_ps.prev_gm),"data.frame") %>%
  select(SampleID, location, Year, Month, month)

# combined distances with sample data
colnames(samples_data) = c("Var1", "location1", "Year1","Month1", "month1")
frac_distances.m = left_join(frac_distances.m, samples_data, by = "Var1")

colnames(samples_data) = c("Var2", "location2","Year2","Month2", "month2")
frac_distances.m = left_join(frac_distances.m, samples_data, by = "Var2")


#remove duplicates
frac_distances.m <- frac_distances.m%>%
  mutate(Var = map2_chr(Var1, Var2, ~toString(sort(c(.x, .y))))) %>%
  distinct(Var, .keep_all = TRUE) %>%
  select(-Var)

#extract distances between fractions by depth
frac_distances_all_type <- frac_distances.m%>%
  mutate(Jump = month1-month2,
         Date1 = factor(paste(Year1,Month1,sep ="_"), levels =c("2013_Jan","2013_Feb","2013_Mar",
                                                                   "2013_Apr","2013_May","2013_Jun",
                                                                   "2013_Jul","2013_Aug","2013_Sep",
                                                                   "2013_Oct","2013_Nov","2013_Dec",
                                                                   "2014_Jan","2014_Feb","2014_Mar",
                                                                   "2014_Apr","2014_May","2014_Jun",
                                                                   "2014_Jul","2014_Aug","2014_Sep",
                                                                   "2014_Oct","2014_Nov","2014_Dec")),
         Date2 = factor(paste(Year2,Month2,sep ="_"), levels =c("2013_Jan","2013_Feb","2013_Mar",
                                                         "2013_Apr","2013_May","2013_Jun",
                                                         "2013_Jul","2013_Aug","2013_Sep",
                                                         "2013_Oct","2013_Nov","2013_Dec",
                                                         "2014_Jan","2014_Feb","2014_Mar",
                                                         "2014_Apr","2014_May","2014_Jun",
                                                         "2014_Jul","2014_Aug","2014_Sep",
                                                         "2014_Oct","2014_Nov","2014_Dec"))) %>% 
  filter(location1 == location2,
         Date1 != Date2,
         Jump  == 1)


frac_distances_all_type$group = factor(interaction(frac_distances_all_type$Type1,
                                                   frac_distances_all_type$Region1),
                                       level= c("SRF.EGC","SRF.WSC",
                                                "EPI.EGC","EPI.WSC",
                                                "MESO.EGC","MESO.WSC",
                                                "BATHY.EGC","BATHY.WSC"))

#test for normal distribution
shapiro.test(frac_distances_all_type$value)

kruskal.test(value ~ group, data = frac_distances_all_type)

dist_Wilcox_frac <- frac_distances_all_type   %>%
  rstatix::wilcox_test(value ~ Type1, p.adjust.method = "BH") %>%
  add_significance()

type_comparisons <- list(c("SRF","EPI"),
                         c("EPI","MESO"),
                         c("MESO","BATHY"))

#plot
dist_type.p<- ggplot(subset(frac_distances_all_type, location1 =="Res."), aes(x = Date1, y = value,group = location1)) +
  labs(x = "Water layer")+
  geom_line()+
  #geom_boxplot(outlier.color = NULL, notch = FALSE, fill = "yellow", color = "black", size = 2)+
  #geom_signif(comparisons = type_comparisons, map_signif_level=TRUE, test = "wilcox.test")+
  theme_classic() +
  #facet_wrap(~location1+Year1, ncol = 2)+
  #scale_fill_manual(values = c("EGC" = "blue", "WSC"="red")) +
  theme(legend.position = "none")


#compare distances between fractions along the water column
frac_distances_all_depth<- frac_distances.m%>%
  filter(Type1 != Type2,
         Fraction1 == Fraction2)

frac_distances_all_depth$group = factor(interaction(frac_distances_all_depth$Fraction1,
                                                    frac_distances_all_depth$Region1),
                                        level= c("FL.EGC","FL.WSC",
                                                 "PA.EGC","PA.WSC"))

#test for normal distribution
shapiro.test(frac_distances_all_depth$value)

kruskal.test(value ~ group, data = frac_distances_all_depth)

dist_Wilcox_frac <- frac_distances_all_depth   %>%
  rstatix::wilcox_test(value ~ group, p.adjust.method = "BH") %>%
  add_significance()


depth_comparisons <- list(c("FL.EGC","FL.WSC"),
                          c("PA.EGC","PA.WSC"))
#plot
dist_frac.p<- ggplot(frac_distances_all_depth, aes(x = group, y = value, fill= Region1)) +
  labs(x = "Fraction")+
  geom_boxplot(outlier.color = NULL, notch = FALSE, color = "black", size = 2)+
  geom_signif(comparisons = depth_comparisons, map_signif_level=TRUE, test = "wilcox.test")+
  scale_fill_manual(values = c("EGC" = "blue", "WSC"="red")) +
  theme_classic() +
  theme(legend.position = "none")


#combined plot
plot_grid(PS99.ord.p, NA,dist_type.p, dist_frac.p, labels = c("A","","B","C"), ncol = 2, align = "h")
