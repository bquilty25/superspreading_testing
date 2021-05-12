#Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from Comix
source("scripts/utils.R")
source("scripts/lee_infectivity.R")

#Load contact data
contacts_bbc_o18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_o18.csv")) %>% 
  mutate(e_repeated=e_home+e_work,
         e_all=e_home+e_work+e_other)

contacts_bbc_u18 <- read.csv(here("2020-cov-tracing","data","contact_distributions_u18.csv"))%>% 
  mutate(e_repeated=e_home+e_work,
         e_all=e_home+e_work+e_other)

contacts_bbc <- bind_rows(contacts_bbc_o18,contacts_bbc_u18) %>% 
  mutate(time_period="pre")

contacts_comix_o18 <- read.csv(here("data","comix_contact_distributions_o18.csv")) %>% 
  mutate(e_repeated=e_home+e_work+e_school,
         e_all=e_home+e_work+e_school+e_other,
         date=as.Date(date)) 

contacts_comix_u18 <- read.csv(here("data","comix_contact_distributions_u18.csv"))%>% 
  mutate(e_repeated=e_home+e_work+e_school,
         e_all=e_home+e_work+e_school+e_other,
         date=as.Date(date)) 

contacts_comix <- bind_rows(contacts_comix_o18,contacts_comix_u18) 

comix_high_low <- contacts_comix %>%
  mutate(year = year(as.Date(date)),
         month = month(as.Date(date))) %>%
  filter(year == 2020 & month %in% c(8, 9)|
         year == 2021 & month %in% c(2, 3)) %>% 
  mutate(time_period=case_when(year == 2020 & month %in% c(8, 9) ~ "aug_sept",
                               year == 2021 & month %in% c(2, 3) ~ "feb_mar"))

#summarise number of contacts
contact_data <- comix_high_low %>% bind_rows(contacts_bbc) 

fitdist(contact_data %>%  filter(time_period=="pre") %>% pull(e_all),"nbinom")
fitdist(contact_data %>%  filter(time_period=="aug_sept") %>% pull(e_all),"nbinom")
fitdist(contact_data %>% filter(time_period=="feb_mar") %>% pull(e_all),"nbinom")

inf_curve <- make_trajectories(n_cases=10,n_sims = 100) %>% 
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
  crossing(trunc_by_onset=c(FALSE)) %>% 
  mutate(
    trunc_t=case_when(trunc_by_onset&type=="symptomatic"~onset_t,
                      trunc_by_onset&type=="asymptomatic"~NA_real_,
                      TRUE~NA_real_),
    infectiousness = pmap(inf_curve_func,.l=list(m=m,trunc_t=trunc_t)))  %>% 
  rowwise() %>% 
  mutate(sum_inf=sum(infectiousness$culture)) %>% 
  ungroup() %>% 
  mutate(norm_sum=(sum_inf-min(sum_inf))/(max(sum_inf)-min(sum_inf)))

prob_infect <- inf_curve %>% 
  select(-u) %>% 
  crossing(time_period=factor(c("pre","aug_sept","feb_mar"),
                              ordered = T,
                              levels=c("pre","aug_sept","feb_mar"))) %>% 
  mutate(contacts_repeated=case_when(time_period== "pre" ~ sample(contact_data %>%  
                                                                    filter(time_period=="pre") %>% 
                                                                    pull(e_repeated),
                                                                  size = n(),
                                                                  replace = T),
                                     time_period== "aug_sept" ~ sample(contact_data %>%  
                                                                         filter(time_period=="aug_sept") %>% 
                                                                         pull(e_repeated),
                                                                       size = n(),
                                                                       replace = T),
                                     time_period == "feb_mar" ~ sample(contact_data %>%  
                                                                         filter(time_period=="feb_mar") %>% 
                                                                         pull(e_repeated),
                                                                       size = n(),
                                                                       replace = T))) %>% 
  mutate(n_repeated_infected=rbinom(n=n(),size=contacts_repeated,prob=norm_sum)) %>% 
  unnest(infectiousness) %>%
  mutate(norm_daily=(culture/sum_inf)*norm_sum,
         ) %>% 
  mutate(contacts_casual=case_when(time_period== "pre"        ~ sample(contact_data %>%  
                                                                    filter(time_period=="pre") %>% 
                                                                    pull(e_other),
                                                                  size = n(),
                                                                  replace = T),
                                     time_period== "aug_sept" ~ sample(contact_data %>%  
                                                                         filter(time_period=="aug_sept") %>% 
                                                                         pull(e_other),
                                                                       size = n(),
                                                                       replace = T),
                                     time_period == "feb_mar" ~ sample(contact_data %>%  
                                                                         filter(time_period=="feb_mar") %>% 
                                                                         pull(e_other),
                                                                       size = n(),
                                                                       replace = T))) %>% 
  ungroup() %>% 
  mutate(n_casual_infected=rbinom(n=n(),size=contacts_casual,prob=norm_daily)) %>% 
  group_by(sim,idx,type,time_period,norm_sum,contacts_repeated,n_repeated_infected) %>% 
  summarise(contacts_casual=sum(contacts_casual),
    n_casual_infected=sum(n_casual_infected)) %>% 
  mutate(n_total_infected=n_repeated_infected+n_casual_infected) 


prob_infect %>% ggplot()+geom_histogram(aes(x=n_casual_infected))+facet_wrap(~time_period)

fitdist(prob_infect %>% filter(time_period=="pre") %>% pull(n_total_infected),"nbinom")
fitdist(prob_infect %>% filter(time_period=="aug_sept") %>% pull(n_total_infected),"nbinom")
fitdist(prob_infect %>% filter(time_period=="feb_mar") %>% pull(n_total_infected),"nbinom")

