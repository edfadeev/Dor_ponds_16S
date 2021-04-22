
ENV <- read.csv("./data/meta_table.csv", h = T, row.names = 1)

fish_data<- read.csv("./data/Fish_data_raw.csv", h=T) %>% 
  separate(col = "Date",
           into = c("j.day", "month", "year"), sep = "/", remove = FALSE) %>% 
  mutate(Pool = ifelse(Pool =="D1", "D1.", "V2."))



food_plots<- list()
food_df<- data.frame()

for (p in c("D1.", "V2.")) {
  for (y in c(2013:2015)) {
measure<- fish_data %>% mutate(type = "Measured") %>% filter(Pool== p, year == y)

# Estimate the rest parameters using a linear model
model.0.s <- lm(Food..Kg.~Aqua_day, data=measure)  
alpha.0.s <- coef(model.0.s)[1]
beta.0.s <- coef(model.0.s)[2]
summary(model.0.s)

predict<- ENV %>% 
  filter(location == p, Year == y, Aqua_day != "NA") %>% 
  select(Aqua_day) %>% 
  filter(!Aqua_day %in% measure$Aqua_day) %>% #no need to predict the dates that have actual measurement
  mutate(type ="Predicted", Food..Kg. = predict(model.0.s,list(Aqua_day= Aqua_day)))

######## Calculate the linear model predictions for the figure
merged<- rbind(measure[,c("Aqua_day","Food..Kg.","type")], 
               predict[,c("Aqua_day","Food..Kg.","type")]) %>% 
          mutate(Food..pred = predict(model.0.s,list(Aqua_day= Aqua_day)))

text<-paste("y =", signif(alpha.0.s,3),"+",signif(beta.0.s,3),"\u00D7 x \n",
            "R^2","=",signif(summary(model.0.s)$adj.r.squared, 3), "\n",
            "p =",signif(summary(model.0.s)$coefficients[2,4],3))

p2 <- ggplot(data = merged, aes(x=Aqua_day, y= Food..Kg.))+
  geom_point(aes(colour = type), size = 5)+
  geom_line(aes(x=Aqua_day, y= Food..pred))+
  ggtitle(paste("Food input ", p," ",y, sep=""))+
  annotate("text", label = text, x= min(merged$Aqua_day)+20, y = max(merged$Food..pred))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
        legend.position = "none", 
        axis.title.x = element_blank())



x<- paste(p,"_",y, sep ="")
#food_plots[[x]]<- add_sub(p2, text)

food_plots[[x]]<- p2

#generate a merged dataframe of all values
food_year<- merged %>% mutate(location = p, Year = y)

food_df<- rbind(food_df, food_year)
                          
  }
  }




#create fish plots
fish_plots<- list()
fish_df<- data.frame()

for (p in c("D1.", "V2.")) {
  for (y in c(2013:2015)) {
    measure<- fish_data %>% mutate(type = "Measured") %>% filter(Pool== p, year == y)
    
    # Estimate the rest parameters using a linear model
    model.0.s <- lm(Biomass..kg.~Aqua_day, data=measure)  
    alpha.0.s <- coef(model.0.s)[1]
    beta.0.s <- coef(model.0.s)[2]
    summary(model.0.s)
    
    predict<- ENV %>% 
      filter(location == p, Year == y, Aqua_day != "NA") %>% 
      select(Aqua_day) %>% 
      filter(!Aqua_day %in% measure$Aqua_day) %>% #no need to predict the dates that have actual measurement
      mutate(type ="Predicted", Biomass..kg. = predict(model.0.s,list(Aqua_day= Aqua_day)))
    
    ######## Calculate the linear model predictions for the figure
    merged<- rbind(measure[,c("Aqua_day","Biomass..kg.","type")], 
                   predict[,c("Aqua_day","Biomass..kg.","type")]) %>% 
      mutate(Biomass..pred = predict(model.0.s,list(Aqua_day= Aqua_day)))
    
    
    text<-paste("y =", signif(alpha.0.s,3),"+",signif(beta.0.s,3),"\u00D7 x \n",
                "R^2","=",signif(summary(model.0.s)$adj.r.squared, 3), "\n",
                "p =",signif(summary(model.0.s)$coefficients[2,4],3))
    
    p2 <- ggplot(data = merged, aes(x=Aqua_day, y= Biomass..kg.))+
      geom_point(aes(colour = type), size = 5)+
      geom_line(aes(x=Aqua_day, y= Biomass..pred))+
      ggtitle(paste("Fish biomass ",p,"-",y, sep=""))+
      annotate("text", label = text, x= min(merged$Aqua_day)+20, y = max(merged$Biomass..pred))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),#axis.text.x = element_blank(),
            legend.position = "none", 
            axis.title.x = element_blank())
    
    x<- paste(p,"_",y, sep ="")
    #fish_plots[[x]]<- add_sub(p2, text, x= 0)
    fish_plots[[x]]<- p2
  
    #generate a merged dataframe of all values
  fish_year<- merged %>% mutate(location = p, Year = y)
  
  fish_df<- rbind(fish_df, fish_year)
  }
  
}


plot_grid(food_plots[["D1._2013"]], food_plots[["D1._2014"]], food_plots[["D1._2015"]],
          food_plots[["V2._2013"]], food_plots[["V2._2014"]], food_plots[["V2._2015"]],
          fish_plots[["D1._2013"]], fish_plots[["D1._2014"]], fish_plots[["D1._2015"]],
          fish_plots[["V2._2013"]], fish_plots[["V2._2014"]], fish_plots[["V2._2015"]],
          ncol = 3)



ggsave("./figures/aqua_pred_lm.pdf", 
       plot = last_plot(),
       units = "cm",
       width = 30, height = 30, 
       #scale = 1,
       dpi = 300)


aqua_pred<- merge(fish_df[,c("location", "Year", "Aqua_day", "Biomass..kg.")],
            food_df[,c("location", "Year", "Aqua_day", "Food..Kg.")]) %>% 
            mutate(Biomass..kg. = ifelse(Biomass..kg.< 0, NA, Biomass..kg.), #omit negative values
                   Food..Kg. = ifelse(Food..Kg.< 0, NA,  Food..Kg.),)

ENV<- left_join(ENV, aqua_pred[,c("location", "Year", "Aqua_day", "Biomass..kg.","Food..Kg.")])

write.csv(ENV, "./data/meta_table.csv")
