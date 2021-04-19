
# AK: load LFT curve to check
# setwd("~/Documents/GitHub/superspreading_testing")
# lft_curve <- read_csv("~/Documents/GitHub/pcr-profile/LFT_curve_summary.csv"); day_list <- seq(0,20,1)
# lft_curve_days <- lft_curve[match(day_list,lft_curve$days_since_infection),]

# Create trajectories
inf_curve <- make_trajectories(n_cases=100,n_sims = 100) %>% 
  mutate(infectiousness = map(inf_curve_func,.x=m)) %>% 
  rowwise() %>% 
  mutate(auc=auc_wrapper(infectiousness)) %>% 
  ungroup() %>% 
  mutate(norm_auc=(auc-min(auc))/(max(auc)-min(auc)))

#Load contact data
contacts_o18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_o18.csv")) %>% 
  mutate(e_repeated=e_home+e_work) %>% 
  select(-c(e_home,e_work))

contacts_u18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_u18.csv"))%>% 
  mutate(e_repeated=e_home+e_work) %>% 
  select(-c(e_home,e_work))

contacts <- bind_rows(contacts_o18,contacts_u18) 

prob_infect <- inf_curve %>% 
  unnest(infectiousness) %>% 
  mutate(shedding=1-culture) %>% 
  mutate(contacts_casual=sample(contacts$e_other,size=n(),replace = T)) %>% 
  group_by(sim,idx,norm_auc,contacts_casual) %>% 
  summarise(prob_infect=1-prod(shedding),
            n_casual_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts_casual,prob=norm_auc))) %>% 
  ungroup() %>% 
  group_by(sim,idx) %>% 
  nest() %>% 
  mutate(contacts_repeated=sample(contacts$e_repeated,size=n(),replace = T)) %>%
  unnest() %>% 
  mutate(n_repeated_infected=rpois(n=n(),lambda=rbinom(n=n(),size=contacts_repeated,prob=prob_infect)),
         n_total_infected=n_repeated_infected+n_casual_infected)

fitdist(prob_infect$n_total_infected,"nbinom")
ggplot(prob_infect)+geom_histogram(aes(x=n_total_infected))
# 
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