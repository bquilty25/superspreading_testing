#Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from Comix
source("scripts/utils.R")
source("scripts/duration.R")


#Make VL trajectories from variant characteristics and asymptomatic fraction
traj <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE,FALSE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = 10000,asymp_parms=asymp_fraction,variant_info=.x,browsing=F)) %>% 
  bind_rows.()

#Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>%
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end
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
   
scenarios <- crossing(time_periods) %>% 
  filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(TRUE,FALSE),
           #heterogen_vl=c(TRUE,FALSE),
           prop_self_iso_test=c(0,0.5,1),
           sampling_freq=c(3,7))

#individual level parameters
traj_scenarios <- crossing(traj,
                 scenarios) %>% 
  mutate.(repeated_contacts = case_when.(heterogen_contacts~sample(contact_data$e_home[contact_data$period==period],
                                     size=n(),
                                     replace=T),
                                     TRUE~round(mean(contact_data$e_home[contact_data$period==period]))),
                                     .by=c(period),
          self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq))

key_grouping_var <- c("sim","variant","scenario_id","period","lower_inf_thresh","heterogen_vl","heterogen_contacts","sampling_freq","prop_self_iso_test")

traj_scenarios_joined <- traj_ %>% 
  left_join.(traj_scenarios)  %>% 
  filter.(lower_inf_thresh==F,heterogen_vl==T,heterogen_contacts==T)

#### Generate infections of repeated (household) contacts ####

repeated_infections <- traj_scenarios_joined %>% 
  uncount.(repeated_contacts,.id="id",.remove = F) %>% 
  mutate.(hh_duration = sample(contacts_hh_duration,size=n(),replace=T),
          infected    = rbernoulli(n(),p=culture_p*hh_duration)) %>% 
  filter.(infected==T) %>% 
  slice.(min(t), .by=c(all_of(key_grouping_var),repeated_contacts,id)) %>% 
  count.(t,all_of(key_grouping_var),repeated_contacts,name = "repeated_infected")

#### Calculate casual infections ####

casual_infections <- traj_scenarios_joined %>% 
  
# Sample daily contacts
  mutate.(casual_contacts = case_when.(heterogen_contacts ~ sample(contact_data$e_other[contact_data$period==period],
                                                                  size=n(),
                                                                  replace=T),
                                       TRUE ~ round(mean(contact_data$e_other[contact_data$period==period]))),
          .by=period) %>%
  
# Testing: determine if and when testing + isolating by specified sampling frequency, adherence
  mutate.(
    test_day = ifelse((t - begin_testing) %% sampling_freq == 0 ,TRUE,FALSE),
    earliest_pos = min(t[which.max(test)&test_day]),
    .by=all_of(key_grouping_var)) %>%
  mutate.(casual_contacts=ifelse(t>=earliest_pos & self_iso_test, 0L, casual_contacts),.by=all_of(key_grouping_var)) %>% 
  uncount.(casual_contacts,.remove=F) %>% 
  
# Simulate infections
  mutate.(nhh_duration=sample(contacts_nhh_duration,size=n(),replace=T),
          infected = rbernoulli(n=n(),p = culture_p*nhh_duration)) %>% 
  count.(t,all_of(key_grouping_var),infected) %>% 
  pivot_wider.(values_from=N,names_from=infected,values_fill = 0) %>% 
  mutate.(casual_contacts=`FALSE`+`TRUE`) %>% 
  select.(everything(),"casual_infected"=`TRUE`,-`FALSE`)
  
# Join casual and repeated contacts and summarise
processed_infections <- traj_scenarios_joined %>% 
  left_join.(casual_infections) %>% 
  left_join.(repeated_infections) %>% 
  replace_na.(list(repeated_infected=0,casual_infected=0,casual_contacts=0)) %>% 
  arrange.(period,lower_inf_thresh) %>% 
  mutate.(
    total_contacts = casual_contacts+repeated_contacts,
    total_infections=casual_infected+repeated_infected)
  
rm(repeated_infections,casual_infections)


source("scripts/results.R")
