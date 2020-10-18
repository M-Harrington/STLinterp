pacman::p_load(tidyverse,dplyr,devtools)
setwd("C:/Users/Matttt/Documents/STLinterp")

#######################Start Below:##################
source_url("https://github.com/M-Harrington/STLinterp/blob/master/STLinterp.R?raw=TRUE")

araria <- read_csv("araria.csv")
araria %>% head()

araria %>% ggplot(aes(time,med.welldepth))+
  geom_line(size=.75,color ="steelblue")+
  geom_point(size=2,shape=1,color="steelblue")+
  labs(title="Well Depth in Araria, Bihar", x= "Season", y="Median Well Depth (m)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
  


#example how to use STL plus to fill in missing values
araria_ts <- ts(araria$med.welldepth, frequency=4)

# which na?
missing_indicies<- which(is.na(araria_ts))

#perform stl, reconstruct missing
stlobj <- stlplus(araria_ts, s.window=4, s.degree=2, t.window=4,
                  t.degree=2,fc.window=4)

reconstruct <- seasonal(stlobj) + trend(stlobj)

#replace NAs in the dataset with interpolated values
araria_reconstructed <- araria %>%  
  mutate(med.welldepth =replace(med.welldepth,
          is.na(med.welldepth), reconstruct[missing_indicies]))


#plot seasonal and trend components
araria$season_stl <- seasonal(stlobj)
araria$trend_stl <- trend(stlobj)

araria_long <- araria %>% pivot_longer(cols=ends_with("stl"), "component")

araria_long %>% ggplot(aes(time, value, col=component))+
  geom_line(size=.75)+
  geom_point(size=2,shape=1)+
  labs(title="Estimated Season vs Trend Components with STL",
       x= "Season", y="Well Depth (m)")+
  scale_color_discrete(name = "Component", labels = c("Seasonal", "Trend"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#Show reconstructed points
names(araria_reconstructed)[3] <- "recon_welldepth" 
araria <- bind_rows(araria, araria_reconstructed[missing_indicies,])

araria %>% ggplot(aes(time,med.welldepth))+
  geom_line(size=.75,color ="steelblue")+
  geom_point(size=2,shape=1,color="steelblue")+
  geom_point(aes(x=time,y=recon_welldepth),color="red", size=2)+
  labs(title="Reconstructed points (Red) and Actual Time Series", x= "Season", y="Well Depth")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#see reconstructed time series
araria_reconstructed %>% ggplot(aes(time,recon_welldepth))+
  geom_line(size=.75,color ="steelblue")+
  geom_point(size=2,shape=1,color="steelblue")+
  labs(title="Reconstructed Well Depth for Araria, Bihar",
       x= "Season", y="Median Well Depth")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

#Looks pretty believable. 
#But how well did we actually do interpolating?  How can we quantify a "good" fit?

#check quality of prediction
params <-c(4,2,4,2,4,1) #s.window , s.degree, t.window , t.degree, fc.window, fc.degree respectively
mc_cross_val(araria_ts, grid_row = params, p=.95,k=1000)

# or
kfold_cross_val(araria_ts,grid_row=params,k = 65)

#but maybe want to test a whole bunch of hypotheses
grid_results <- STLinterp(araria_ts, s.window = c(5,6,8), 
          t.window = c(5,10,15),t.degree = c(1,2), s.degree = c(1,2), k=65,
          type="kfold", value = "grid")
grid_results %>% head

#or return just the best parameter set
STLinterp(araria_ts, s.window = c(5,6,8), 
          t.window = c(5,10,15),t.degree = c(1,2), s.degree = c(1,2), k=65,
          type="kfold", value = "best")
