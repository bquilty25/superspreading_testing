
# AK: load LFT curve to check
# setwd("~/Documents/GitHub/superspreading_testing")
# lft_curve <- read_csv("~/Documents/GitHub/pcr-profile/LFT_curve_summary.csv"); day_list <- seq(0,20,1)
# lft_curve_days <- lft_curve[match(day_list,lft_curve$days_since_infection),]
source("scripts/utils.R")
source("scripts/lee_infectivity.R")

make_trajectories(n_cases=10,n_sims = 100) %>% 
  mutate(sampling_freq=3) %>% 
  mutate(test_times=pmap(.f=test_times,list(sampling_freq=sampling_freq,onset_t=onset_t,type=type))) %>% 
  unnest(test_times) %>% 
  mutate(
    ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
    test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
    test_label = detector(test_p = test_p,  u = runif(n=n(),0,1))) %>% 
  nest(ct,test_t,test_no,test_p,test_label) %>%
  mutate(earliest_positive=map(.f=earliest_pos,.x=data)) %>%
  unnest(earliest_positive) %>% 
  select(-data) %>% 
  mutate(infectiousness = map(inf_curve_func,.x=m)) %>% 
  rowwise() %>% 
  mutate(auc=MESS::auc(y=infectiousness$culture,x=infectiousness$t),
         test_auc=ifelse(is.infinite(test_t),
                         0,
                         MESS::auc(y=infectiousness$culture,x=infectiousness$t,from=test_t)),
         symp_auc=ifelse(type=="asymptomatic",
                         0,
                         MESS::auc(y=infectiousness$culture,x=infectiousness$t,from=onset_t)),
         #test_auc=0,
         #symp_auc=0,
         trunc_auc=pmax(test_auc,symp_auc),
         #prop=(trunc_auc/auc)
         ) %>%
  #filter(type=="symptomatic") %>% 
  ungroup() %>%  
  group_by(sim) %>% 
  summarise(prop=sum(trunc_auc)/sum(auc)) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$prop*100,
                                              probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q)
  #ungroup() %>% 
  #mutate(norm_auc=((auc*(1-prop))-min((auc*(1-prop))))/(max(auc*(1-prop))-min(auc*(1-prop))))
  #mutate(norm_auc=(auc-min(auc))/(max(auc)-min(auc)))

#Load contact data
contacts_o18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_o18.csv")) %>% 
  mutate(e_repeated=e_home+e_work) %>% 
  select(-c(e_home,e_work))

contacts_u18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_u18.csv"))%>% 
  mutate(e_repeated=e_home+e_work) %>% 
  select(-c(e_home,e_work))

contacts <- bind_rows(contacts_o18,contacts_u18) 

prob_infect <- inf_curve %>% select(-u) %>% 
  unnest(infectiousness) %>% 
  mutate(contacts_casual=sample(contacts$e_other,size=n(),replace = T)) %>% 
  group_by(sim,idx,type,norm_auc,contacts_casual) %>% 
  summarise(prob_infect=1-prod(1-culture),
            n_casual_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts_casual,prob=norm_auc))) %>% 
  ungroup() %>% 
  group_by(sim,idx) %>% 
  nest() %>% 
  mutate(contacts_repeated=sample(contacts$e_repeated,size=n(),replace = T)) %>%
  unnest() %>% 
  mutate(n_repeated_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts_repeated,prob=prob_infect)),
         n_total_infected=n_repeated_infected+n_casual_infected)

(r_k <- fitdist(prob_infect$n_total_infected,"nbinom"))
ggplot(prob_infect)+geom_histogram(aes(x=n_total_infected))

#Prop responsible for 80% of transmission
propresponsible(R0=r_k$estimate[2],r_k$estimate[1],prop=0.8)



# contacts_plot <- as_tibble(contacts_sample)%>% 
#   mutate(value=ifelse(round(value)>20,">20",round(value)),
#          value=fct_relevel(value,c(as.character(seq(0,20)),">20"))) %>% 
#   ggplot(aes(x=value,
#              y=..prop..,
#              group=1))+
#   geom_bar()+
#   scale_x_discrete(breaks=c(as.character(seq(0,18,by=2)),">20"))+
#   labs(x="Number of contacts",y="Probability")+
#   theme_minimal()
# 
# auc_dist <- fitdistrplus::fitdist(inf_curve$tot_auc,"gamma")
# auc_sample <- rgamma(n=100000,shape=auc_dist$estimate[1],rate=auc_dist$estimate[2])
# 
# auc_plot <- as_tibble(auc_sample)%>% 
#   ggplot(aes(x=value),fill="grey")+
#   geom_density()+
#   lims(x=c(0,10))+
#   labs(x="Area under infectiousness curve",y="Density")+
#   theme_minimal()
# 
# R_sample <- contacts_sample*auc_sample
# (R=fitdist(round(R_sample),"nbinom"))
# 
# sec_cases_plot <- as_tibble(R_sample) %>% 
#   mutate(value=ifelse(round(value)>20,">20",round(value)),
#          value=fct_relevel(value,c(as.character(seq(0,20)),">20"))) %>% 
#   ggplot(aes(x=value,
#              y=..prop..,
#              group=1))+
#   geom_bar()+
#   labs(x="Number of secondary cases",y="Probability")+
#   theme_minimal()+
#   annotate(label=paste0("R = ",as.character(round(R$estimate[2],2)),"\n","k = ",as.character(round(R$estimate[1],2))),
#            x="17",
#            y=0.4,
#            geom="text")
# 
# (auc_plot+contacts_plot)/
#   sec_cases_plot
# 
# ggsave("results/plot.png",width=210,height=120,units="mm",dpi=600)
#        