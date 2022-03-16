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
   
time_periods_of_interest <- crossing(time_periods) %>% 
  filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T))

testing_scenarios <- traj %>% 
  select.(-m) %>% 
  crossing.(prop_self_iso_test=c(0,0.5,1),
                              sampling_freq=c(3,7)) %>% 
  mutate.(self_iso_test = rbernoulli(n=n(),prop_self_iso_test),
          begin_testing = rdunif(n(),0, sampling_freq)) 

key_grouping_var <- c("sim","variant","scenario_id","period","lower_inf_thresh","heterogen_vl","heterogen_contacts")

#### Generate infections of repeated (household) contacts ####

indiv_params <- traj %>% 
  select.(-m) %>% 
  crossing.(time_periods_of_interest) %>% 
  filter.(heterogen_vl==T,heterogen_contacts==T) %>% 
  mutate.(repeated_contacts = case_when.(heterogen_contacts~sample(contact_data$e_home[contact_data$period==period],
                                     size=n(),
                                     replace=T),
                                     TRUE~round(mean(contact_data$e_home[contact_data$period==period]))),
                                     .by=period)

indiv_params_long <- indiv_params %>% 
  left_join.(traj_)

repeated_infections <- indiv_params_long %>% 
  uncount.(repeated_contacts,.id="id",.remove = F) %>% 
  mutate.(hh_duration = sample(contacts_hh_duration,size=n(),replace=T),
          infected    = rbernoulli(n(),p=culture_p*hh_duration)) %>% 
  filter.(infected==T) %>% 
  slice.(min(t), .by=c(all_of(key_grouping_var),repeated_contacts,id)) %>% 
  count.(t,all_of(key_grouping_var),repeated_contacts,name = "repeated_infected")

#### Calculate casual infections ####

casual_infections <- indiv_params_long %>% 
  
# Sample daily contacts
  mutate.(casual_contacts = case_when.(heterogen_contacts ~ sample(contact_data$e_other[contact_data$period==period],
                                                                  size=n(),
                                                                  replace=T),
                                       TRUE ~ round(mean(contact_data$e_other[contact_data$period==period]))),
          .by=period) %>%
  
# Simulate infections 
  uncount.(casual_contacts) %>% 
  mutate.(nhh_duration = sample(contacts_nhh_duration,size=n(),replace=T),
          casual_infected = rbernoulli(n=n(),p = culture_p*nhh_duration)) %>% 
  
# Testing: determine if and when testing + isolating by specified sampling frequency, adherence  
  
  left_join.(testing_scenarios) %>% 
  mutate.(
    test_day = ifelse((t - begin_testing) %% sampling_freq == 0 ,TRUE,FALSE),
    earliest_pos = min(t[which.max(test)&test_day]),
    isolating = t>=earliest_pos & self_iso_test,
    .by=c(all_of(key_grouping_var),prop_self_iso_test,sampling_freq)) %>%
  replace_na.(list(isolating=FALSE)) %>% 
 
  count.(t,all_of(key_grouping_var),prop_self_iso_test,sampling_freq,casual_infected,isolating) %>%
  filter.(isolating==F) %>% 
  pivot_wider.(values_from=N,names_from=casual_infected,values_fill = 0) %>% 
  mutate.(casual_contacts=`FALSE`+`TRUE`) %>% 
  select.(everything(),"casual_infected"=`TRUE`,-`FALSE`)
  
# Join casual and repeated contacts and summarise
processed_infections <- indiv_params_long %>% 
  right_join.(casual_infections) %>% 
  left_join.(repeated_infections) %>% 
  replace_na.(list(repeated_infected=0,casual_infected=0,casual_contacts=0)) %>% 
  arrange.(period,lower_inf_thresh) %>% 
  mutate.(
    total_contacts = casual_contacts+repeated_contacts,
    total_infections=casual_infected+repeated_infected)
  
rm(repeated_infections,casual_infections)


#source("scripts/results.R")
