
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
  crossing(delay=c(-1,0,1,2)) %>% 
  rowwise() %>% 
  mutate(auc=MESS::auc(y=infectiousness$culture,x=infectiousness$t),
         test_auc=ifelse(is.infinite(test_t),
                         0,
                         MESS::auc(y=infectiousness$culture,x=infectiousness$t,from=test_t)),
         symp_auc=ifelse(type=="asymptomatic",
                         0,
                         MESS::auc(y=infectiousness$culture,x=infectiousness$t,from=onset_t+delay)),
         test_auc=0,
         #symp_auc=0,
         trunc_auc=pmax(test_auc,symp_auc),
         prop=(trunc_auc/auc)
  ) %>%
  filter(type=="symptomatic") %>% 
  ungroup() %>%  
  #group_by(sim,delay) %>% 
  #summarise(prop=sum(trunc_auc)/sum(auc)) %>% 
  group_by(delay) %>% 
  nest() %>%
  mutate(Q = purrr::map(.x = data, ~quantile( .$prop*100,
                                              probs = c(0.5,0.025,0.975)))) %>%
  unnest_wider(Q)
#ungroup() %>% 
#mutate(norm_auc=((auc*(1-prop))-min((auc*(1-prop))))/(max(auc*(1-prop))-min(auc*(1-prop))))
#mutate(norm_auc=(auc-min(auc))/(max(auc)-min(auc)))
