infectivity <- read.csv(here("data","lee_infectivity.csv"),col.names = c("ct","y"))

infectivity %<>% 
  arrange(-ct) %>% 
  mutate(x=cumsum(y)/sum(y)) %>% 
  mutate(csum=(x-min(x))/(max(x)-min(x)))

lee_mod <- glm(csum~ct,data=infectivity,family="binomial")

#Plot the original points
# first argument is the x values, second is the y values
ggplot(infectivity)+geom_line(aes(x=ct,y=csum))


