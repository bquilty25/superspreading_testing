#Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from Comix
source("scripts/utils.R")
source("scripts/lee_infectivity.R")

#Load contact data
contacts_bbc_o18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_o18.csv")) 

contacts_bbc_u18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_u18.csv")) 

contacts_bbc <- bind_rows(contacts_bbc_o18, contacts_bbc_u18) %>%
  mutate(time_period = "pre")

contacts_comix_o18 <-
  read.csv(here("data", "comix_contact_distributions_o18.csv")) %>%
  mutate(
    date = as.Date(date)
  )

contacts_comix_u18 <-
  read.csv(here("data", "comix_contact_distributions_u18.csv")) %>%
  mutate(
    date = as.Date(date)
  )

contacts_comix <- bind_rows(contacts_comix_o18, contacts_comix_u18)

comix_high_low <- contacts_comix %>%
  mutate(year = year(as.Date(date)),
         month = month(as.Date(date))) %>%
  filter(year == 2020 & month %in% c(8, 9) |
           year == 2021 & month %in% c(1, 2)) %>%
  mutate(
    time_period = case_when(
      year == 2020 & month %in% c(8, 9) ~ "aug_sept",
      year == 2021 & month %in% c(1, 2) ~ "jan_feb"
    )
  )

#summarise number of contacts
contact_data <- comix_high_low %>% bind_rows(contacts_bbc)

# contact_data %>% 
#   group_by(time_period) %>% 
#   nest() %>% 
#   mutate(dists=map(.x=data,.f= .%>% 
#                      mutate(e_all = rowSums(across(c(e_home,e_school,e_work,e_other)),na.rm = T)) %>% 
#                      pull(e_all) %>% 
#                      fitdist("nbinom"))) %>% 
#   mutate(dist_params=map(dists,~tibble(mu=.x$estimate[[2]],
#                                        size=.x$estimate[[1]]))) %>% 
#   ungroup() %>% 
#   mutate(dist_ci=future_map(dists,~boot_ci(.x),.progress = T))

inf_curve <- make_trajectories(n_cases = 100, n_sims = 100) %>%
  #testing
  crossing(sampling_freq = c(NA,3)) %>%
  mutate(test_times = pmap(
    .f = test_times,
    list(
      sampling_freq = sampling_freq,
      onset_t = onset_t,
      type = type
    )
  )) %>%
  unnest(test_times) %>%
  mutate(
    ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
    test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
    test_label = detector(test_p = test_p,  u = runif(n = n(), 0, 1))
  ) %>%
  nest(ct, test_t, test_no, test_p, test_label) %>%
  mutate(earliest_positive = map(.f = earliest_pos, .x = data)) %>%
  unnest(earliest_positive) %>%
  select(-data) %>%
  mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m)))  %>%
  #proportion who self-isolate
  crossing(prop_self_iso = c(0,0.25,0.5,0.75,1)) %>%
  filter(!(prop_self_iso!=0&&type=="asymptomatic")) %>% 
  mutate(self_iso=rbinom(n=n(),size=1,prob=prop_self_iso)) %>% 
  rowwise() %>%
  mutate(sum_inf = sum(infectiousness$culture)) %>%
  ungroup() %>%
  mutate(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf)))

prob_infect <- 
  inf_curve %>%
  select(-u) %>%
  crossing(time_period = factor(
    c("pre", "aug_sept", "jan_feb"),
    ordered = T,
    levels = c("pre", "aug_sept", "jan_feb")
  )) %>%
  mutate(
    contacts_repeated_home = case_when(
      time_period == "pre" ~ sample(
        contact_data %>%
          filter(time_period == "pre") %>%
          pull(e_home),
        size = n(),
        replace = T
      ),
      time_period == "aug_sept" ~ sample(
        contact_data %>%
          filter(time_period == "aug_sept") %>%
          pull(e_home),
        size = n(),
        replace = T
      ),
      time_period == "jan_feb" ~ sample(
        contact_data %>%
          filter(time_period == "jan_feb") %>%
          pull(e_home),
        size = n(),
        replace = T
      )
    ),
    contacts_repeated_school_work = case_when(
      time_period == "pre" ~ sample(
        contact_data %>%
          filter(time_period == "pre") %>%
          pull(e_work),
        size = n(),
        replace = T
      ),
      time_period == "aug_sept" ~ sample(
        contact_data %>%
          filter(time_period == "aug_sept") %>%
          mutate(e_school_work =  e_work + e_school) %>%
          pull(e_school_work),
        size = n(),
        replace = T
      ),
      time_period == "jan_feb" ~ sample(
        contact_data %>%
          filter(time_period == "jan_feb") %>%
          mutate(e_school_work =  e_work + e_school) %>%
          pull(e_school_work),
        size = n(),
        replace = T
      )
    ),
    trunc_t=case_when(
      # if symptomatic, adhering to self isolation, and either not tested or test neg,
      # truncate at onset
      type=="symptomatic"&is.infinite(test_t)&self_iso!=0~onset_t,
      # if symptomatic, adhering to self isolation, and have onset before test, truncate at onset
      type=="symptomatic"&is.finite(test_t)&onset_t<test_t&self_iso!=0~onset_t,
      # if symptomatic, adhering to self isolation, and have onset after pos test, truncate at test
      type=="symptomatic"&is.finite(test_t)&test_t<onset_t&self_iso!=0~test_t,
      # if symptomatic, not adhering to self isolation, and have a positive test, truncate at test
      TRUE ~ test_t), 
    #if symp onset, contacts outside of home should cease; for school/work, multiply 
    #number of contacts by proportion of time since infection which occurs pre-onset
    contacts_repeated_school_work=case_when(!is.infinite(trunc_t)~ceiling(contacts_repeated_school_work*(trunc_t/(end-start))),
                                            TRUE~ceiling(contacts_repeated_school_work))
  ) %>% 
  mutate(
    contacts_repeated   = contacts_repeated_home + contacts_repeated_school_work,
    n_repeated_infected = rbinom(n = n(), 
                                      size = contacts_repeated, 
                                      prob = norm_sum)) %>%
  unnest(infectiousness) %>%
  mutate(norm_daily = (culture / sum_inf) * norm_sum) %>%
  mutate(
    contacts_casual = case_when(
      t>trunc_t  ~ 0L,
      time_period == "pre"        ~ sample(
        contact_data %>%
          filter(time_period == "pre") %>%
          pull(e_other),
        size = n(),
        replace = T
      ),
      time_period == "aug_sept" ~ sample(
        contact_data %>%
          filter(time_period == "aug_sept") %>%
          pull(e_other),
        size = n(),
        replace = T
      ),
      time_period == "jan_feb" ~ sample(
        contact_data %>%
          filter(time_period == "jan_feb") %>%
          pull(e_other),
        size = n(),
        replace = T
      )
    )
  ) %>%
  ungroup() %>%
  mutate(n_casual_infected = rbinom(n = n(), 
                                    size = contacts_casual, 
                                    prob =
                                      norm_daily)) %>%
  group_by(sim,
           idx,
           type,
           time_period,
           norm_sum,
           prop_self_iso,
           contacts_repeated,
           n_repeated_infected,
           sampling_freq) %>%
  summarise(
    contacts_casual = sum(contacts_casual),
    n_casual_infected = sum(n_casual_infected)
  ) %>%
  mutate(n_total_infected = n_repeated_infected + n_casual_infected)

dists <- prob_infect %>% 
  group_by(time_period,
         prop_self_iso,
         sampling_freq) %>% 
  nest() %>% 
  mutate(dists=map(.x=data,.f= .%>% 
                     pull(n_total_infected) %>% 
                     fitdist("nbinom"))) %>% 
  mutate(dist_params=map(dists,~tibble(mu=.x$estimate[[2]],
                                      size=.x$estimate[[1]]))) %>% 
  ungroup() %>% 
  mutate(dist_ci=future_map(dists,~boot_ci(.x),.progress = T))
  unnest_wider(dist_params) %>%  
  mutate(sampling_freq=case_when(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing"))

dists %>% 
  ungroup() %>% 
  select(-c(data,dists),"R"=mu,"k"=size) %>%
  mutate(estimate=paste("R = ",round(R,2),"k = ",round(k,2))) %>% 
  select(-c(R,k)) %>% 
  pivot_wider(names_from = time_period,values_from = estimate) %>% 
  #arrange(time_period,prop_self_iso,sampling_freq) %>% 
  base::print(n=Inf)

p1 <- prob_infect %>% 
  filter(prop_self_iso%in%c(0,0.5,1)) %>% 
  mutate(sampling_freq=case_when(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing")) %>% 
  mutate(n_total_infected=factor(ifelse(round(n_total_infected)>=10,"\u2265 10",round(n_total_infected))),
         n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
  ggplot(aes(x = n_total_infected,
             y = ..prop..,
             fill = factor(sampling_freq),
             group = 1))+
  geom_bar(width = 0.8)+
  geom_label(data = dists %>% 
               filter(prop_self_iso%in%c(0,0.5,1)),
             aes(label=paste0("R = ",round(mu,2),"\n","k = ",round(size,2)),
                 x="\u2265 10",
                 y=0.75),
             colour = "white",
             fill = rgb(0,0,0,0.4),
             hjust = 1)+
  scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
  labs(x="Number of secondary cases",y="Probability")+
  theme_minimal()+
  theme(axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        axis.line.y.left = element_line())+
  facet_rep_grid(prop_self_iso~time_period+sampling_freq, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                   prop_self_iso=function(x){scales::percent(as.numeric(x))})) + 
  scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))

ggsave(plot = p1,"results/contacts_over_time_prop_self_iso_0_1.png",width=12,height=5,units="in",dpi=400,scale=1.2)

s1 <- prob_infect %>% 
  mutate(sampling_freq=case_when(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing")) %>% 
  mutate(n_total_infected=factor(ifelse(round(n_total_infected)>=10,"\u2265 10",round(n_total_infected))),
         n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
  ggplot(aes(x = n_total_infected,
             y = ..prop..,
             fill = factor(sampling_freq),
             group = 1))+
  geom_bar(width = 0.8)+
  geom_label(data = dists,
             aes(label=paste0("R = ",round(mu,2),"\n","k = ",round(size,2)),
             x="\u2265 10",
             y=0.75),
             colour = "white",
             fill = rgb(0,0,0,0.4),
            hjust = 1)+
  scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
  labs(x="Number of secondary cases",y="Probability")+
  theme_minimal()+
  theme(axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        axis.line.y.left = element_line())+
  facet_rep_grid(prop_self_iso~time_period+sampling_freq, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                   prop_self_iso=function(x){scales::percent(as.numeric(x))})) +  
  scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))

ggsave(plot = s1,"results/contacts_over_time_prop_self_iso_full.png",width=12,height=7,units="in",dpi=400,scale=1.2)

#repeated/casual split
prob_infect %>% 
  ungroup() %>% 
  group_by(time_period,
           prop_self_iso,
           sampling_freq) %>% 
  summarise(prop_repeated = sum(n_repeated_infected) / sum(n_total_infected),
            prop_casual = sum(n_casual_infected) / sum(n_total_infected)) %>% 
  base::print(n=Inf)
