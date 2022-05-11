#Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from Comix
source("scripts/utils.R")
source("scripts/duration.R")

N_sims <- 10000
#Make VL trajectories
traj <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE,FALSE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = N_sims,asymp_parms=asymp_fraction,variant_info=.x,browsing=F)) %>% 
  bind_rows.()

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>%
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end, interval = 1
  )))  %>%
  unnest.(infectiousness) %>%
  crossing.(
    lower_inf_thresh = c(FALSE)
  ) %>%
  mutate.(
    culture_p        = stats::predict(
      object = inf_model_choice(lower_inf_thresh),
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    infectious = rbernoulli(n = n(),
                                  p = culture_p),
    test_p           = stats::predict(
      object =  innova_mod,
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    test       = rbernoulli(n = n(),
                                  p = test_p),
    .by = c(lower_inf_thresh)
  ) %>%
  replace_na.(list(test       = FALSE,
                   infectious = FALSE)) %>% 
  select.(-c(prolif, start, end))
 
#Scenarios to investigate


key_grouping_var <- c("sim","variant","scenario_id","period","lower_inf_thresh","heterogen_vl","heterogen_contacts")

#baseline 

testing_scenarios <- traj %>% 
  filter.(heterogen_vl==T) %>% 
  select.(-m) %>% 
  crossing.(prop_self_iso_test=c(0),
            sampling_freq=c(7),
            event_size=NA) %>% 
  mutate.(self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq)) 

time_periods_of_interest <- 
  crossing(time_periods) %>% 
  filter(date_end<as.Date("2021-01-01"),period!="POLYMOD") %>% 
  #filter(period!="POLYMOD") %>% 
  #filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T))

processed_infections_baseline <- time_periods_of_interest %>% run_model(.,browsing = F)

#heterogen onoff 

testing_scenarios <- traj %>% 
  select.(-m) %>% 
  crossing.(prop_self_iso_test=c(0),
            sampling_freq=c(7),
            event_size=NA) %>% 
  mutate.(self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq)) 

time_periods_of_interest <- 
  crossing(time_periods) %>% 
  filter(date_end<as.Date("2021-01-01"),period!="POLYMOD") %>% 
  #filter(period!="POLYMOD") %>% 
  #filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T,F))

processed_infections_heterogen_on_off <- time_periods_of_interest %>% run_model

#testing

testing_scenarios <- traj %>% 
  filter.(heterogen_vl==T) %>% 
  select.(-m) %>% 
  crossing.(prop_self_iso_test=c(0,.5,1),
            sampling_freq=c(1,3,7),
            event_size=NA) %>% 
  mutate.(self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq)) 

time_periods_of_interest <- 
  crossing(time_periods) %>% 
  filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T))

processed_infections_testing <- time_periods_of_interest %>% run_model(.,browsing = F)

#event testing

testing_scenarios <- traj %>% 
  filter.(heterogen_vl==T) %>% 
  select.(-m) %>% 
  crossing.(prop_self_iso_test=c(0,.5,1),
            sampling_freq=c(NA),
            event_size=c(5,10,20,NA)) %>% 
  filter.(!(prop_self_iso_test>0&is.na(sampling_freq)&is.na(event_size))) %>% 
  mutate.(self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq)) 

time_periods_of_interest <- 
  crossing(time_periods) %>% 
  filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T))

processed_infections_events <- time_periods_of_interest %>% run_model(.,browsing = F)

#rm(repeated_infections,casual_infections)


#source("scripts/results.R")
