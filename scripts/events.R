
# AK: load LFT curve to check
# setwd("~/Documents/GitHub/superspreading_testing")
# lft_curve <- read_csv("~/Documents/GitHub/pcr-profile/LFT_curve_summary.csv"); day_list <- seq(0,20,1)
# lft_curve_days <- lft_curve[match(day_list,lft_curve$days_since_infection),]
source("scripts/utils.R")
source("scripts/lee_infectivity.R")


crowd_contacts <- crossing(density=c(1:5),
               radius=c(5,10,15,20)) %>% 
      mutate(area=pi*(radius^2),
             contacts=area*density)


simul <- make_trajectories(n_cases=100,n_sims = 100) %>% 
  mutate(event_duration=2/24,
         event_start=runif(n(),start,end),
         event_end=event_start+event_duration) %>% 
  crossing(testing=c(TRUE,FALSE)) %>% 
  mutate(
    ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = event_start)),
    test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
    test_label = ifelse(testing,detector(test_p = test_p,  u = runif(n=n(),0,1)),FALSE)) %>% 
  mutate(infectiousness = map(inf_curve_func,.x=m)) %>% 
  rowwise() %>% 
  mutate(auc=MESS::auc(y=infectiousness$culture,x=infectiousness$t,from = event_start,to=event_end),
         event_auc=case_when(
           #Test positive and symptomatic, avert event transmission
           test_label & type=="asymptomatic" ~ 
                               MESS::auc(y=infectiousness$culture,x=infectiousness$t,from = event_start,to=event_end),
           #Test positive and symptomatic before event start, avert event transmission
                             test_label & type=="symptomatic"&event_start>onset_t~
                               MESS::auc(y=infectiousness$culture,x=infectiousness$t,from = event_start,to=event_end),
           #Test positive and symptomatic after event start, avert from event start onwards
                             test_label & type=="symptomatic"&event_start<onset_t~
                               MESS::auc(y=infectiousness$culture,x=infectiousness$t,from = event_start,to=event_end),
           #Test negative and symptomatic before event start, avert event transmission
                             !test_label & type=="symptomatic"&event_start>onset_t~
                               MESS::auc(y=infectiousness$culture,x=infectiousness$t,from = event_start,to=event_end),
           #Test negative and symptomatic after event, don't avert transmission
                             !test_label & type=="symptomatic"&event_start<onset_t~
                               0,
           #Test negative and asymptomatic, don't avert transmission
                             !test_label & type=="asymptomatic" ~ 
                               0)) 

#Proportion of transmission prevented
simul %>% 
  group_by(sim,testing) %>% 
  summarise(prop=sum(event_auc)/sum(auc)) %>% 
  group_by(testing) %>%
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile(.$prop,
                                             probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q)

#Ratio vs symptomatic self-isolation alone
simul %>% 
  group_by(sim,testing) %>% 
  summarise(prop=sum(event_auc)/sum(auc)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=c(testing),names_prefix = "testing_",values_from=prop) %>% 
  mutate(ratio=testing_TRUE/testing_FALSE) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$ratio,
                                              probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q)


#Crowd contacts
simul %>%
  ungroup() %>% 
  mutate(norm_auc=(event_auc-min(event_auc))/(max(event_auc)-min(event_auc))) %>% 
  group_by(sim,idx,testing) %>% 
  crossing(crowd_contacts) %>% 
  mutate(contacts=round(contacts)) %>% 
  mutate(n_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts,prob=norm_auc))) %>% View()
  group_by(density,radius,testing,contacts) %>%
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile(.$n_infected/contacts,
                                             probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q) %>% 
  ggplot(aes(x=testing,y=`50%`,ymin=`2.5%`,ymax=`97.5%`))+
  geom_pointrange()+
  facet_grid(radius~density,labeller = label_both)
  
  #Crowd contacts
n_infected <-  simul %>%
    mutate(spent_auc=auc-event_auc) %>%
    ungroup() %>% 
    mutate(norm_auc=(spent_auc-min(spent_auc))/(max(spent_auc)-min(spent_auc))) %>% 
    group_by(sim,idx,testing) %>% 
    crossing(crowd_contacts
             #%>%  filter(radius==5,density==5)
             ) %>% 
    mutate(contacts=round(contacts)) %>% 
    mutate(n_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts,prob=norm_auc))) 


n_infected %>%  group_by(radius,density,testing) %>% nest() %>% mutate(dist=map(.f=function(x){fitdist(data=x$n_infected,distr="nbinom")},.x=data)) %>% unlist(dist)
(r_k <- fitdist(n_infected%>%  group_by(radius==5,density==5) %>% filter(!testing) %>% pull(n_infected),"nbinom"))

n_infected %>% 
  ggplot(aes(x=n_infected))+geom_histogram()+facet_wrap(~testing,ncol=1)


n_infected %>% 
  mutate(prob_ss=n_infected>10) %>% 
  group_by(testing,radius,density) %>% 
  count(prob_ss) %>% 
  group_by(radius,density,testing) %>% 
  mutate(prob = prop.table(n)) %>% 
  filter(prob_ss) %>% 
  arrange(radius,density,testing) %>% 
  ggplot()+
  geom_tile(aes(x=density,y=radius,fill=prob))+
  geom_text(aes(x=density,y=radius,label=round(prob,2)),colour="white")+
  scale_fill_viridis_c(name="Probability of superspreading event\n(>10 crowd members exposed)",guide=guide_coloursteps(show.limits = F),option="rocket",begin=0.1,end=0.9)+
  facet_wrap(~testing,labeller=labeller(testing=function(x){ifelse(x,"LFT testing before event","No LFT testing before event")}),ncol = 1)+
  coord_fixed(ratio=5/25)+
  labs(x="Crowd density (people per square metre)",
       y="Dispersion distance (radius around infected person)")+
  theme_minimal()+
  theme(legend.position = "bottom")

ggsave("density_distance_ss.png",dpi=600,height = 210,width=197,units="mm")  

