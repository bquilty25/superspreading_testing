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
           year == 2021 & month %in% c(2, 3)) %>%
  mutate(
    time_period = case_when(
      year == 2020 & month %in% c(8, 9) ~ "aug_sept",
      year == 2021 &
        month %in% c(2, 3) ~ "feb_mar"
    )
  )

#summarise number of contacts
contact_data <- comix_high_low %>% bind_rows(contacts_bbc)

fitdist(contact_data %>%  filter(time_period == "pre") %>% pull(e_all),
        "nbinom")
fitdist(contact_data %>%  filter(time_period == "aug_sept") %>% pull(e_all),
        "nbinom")
fitdist(contact_data %>% filter(time_period == "feb_mar") %>% pull(e_all),
        "nbinom")

inf_curve <- make_trajectories(n_cases = 100, n_sims = 100) %>%
  mutate(sampling_freq = 3) %>%
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
  crossing(trunc_by_onset = c(T,F)) %>%
  mutate(
    trunc_by_onset = case_when(
      trunc_by_onset & type == "symptomatic" ~ T,
      trunc_by_onset & type == "asymptomatic" ~ F,
      TRUE ~ F
    ),
    infectiousness = pmap(inf_curve_func, .l = list(m = m))
  )  %>%
  rowwise() %>%
  mutate(sum_inf = sum(infectiousness$culture)) %>%
  ungroup() %>%
  mutate(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf)))

prob_infect <- 
  inf_curve %>%
  select(-u) %>%
  crossing(time_period = factor(
    c("pre", "aug_sept", "feb_mar"),
    ordered = T,
    levels = c("pre", "aug_sept", "feb_mar")
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
      time_period == "feb_mar" ~ sample(
        contact_data %>%
          filter(time_period == "feb_mar") %>%
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
      time_period == "feb_mar" ~ sample(
        contact_data %>%
          filter(time_period == "feb_mar") %>%
          mutate(e_school_work =  e_work + e_school) %>%
          pull(e_school_work),
        size = n(),
        replace = T
      )
    ),
    #if symp onset, contacts outside of home should cease; for school/work, multiply 
    #number of contacts by proportion of infection pre-onset
    contacts_repeated_school_work=ifelse(trunc_by_onset,
                                          ceiling(contacts_repeated_school_work*(onset_t/(end-start))),
                                          contacts_repeated_school_work)
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
      trunc_by_onset&t>onset_t~0L,
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
      time_period == "feb_mar" ~ sample(
        contact_data %>%
          filter(time_period == "feb_mar") %>%
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
           trunc_by_onset,
           contacts_repeated,
           n_repeated_infected) %>%
  summarise(
    contacts_casual = sum(contacts_casual),
    n_casual_infected = sum(n_casual_infected)
  ) %>%
  mutate(n_total_infected = n_repeated_infected + n_casual_infected)

dists <- prob_infect %>% 
  group_by(time_period,
         trunc_by_onset) %>% 
  nest() %>% 
  mutate(dists=map(.x=data,.f= .%>% 
                     pull(n_total_infected) %>% 
                     fitdist("nbinom"))) %>% 
  mutate(dist_params=map(dists,~tibble(mu=.x$estimate[[2]],
                                      size=.x$estimate[[1]]))) %>% 
  unnest_wider(dist_params)


prob_infect %>% 
  mutate(n_total_infected=factor(ifelse(round(n_total_infected)>=20,"\u2265 20",round(n_total_infected))),
         n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,19)),"\u2265 20"))) %>%
  #filter(!trunc_by_onset) %>% 
  ggplot(aes(x=n_total_infected,
             y=..prop..,
             group=1))+
  geom_bar(width=0.8)+
  geom_text(data=dists,
             aes(label=paste0("R = ",as.character(round(mu,2)),"\n","k = ",as.character(round(size,2))),
             x="18",
             y=0.5),
            hjust = 1)+
  labs(x="Number of secondary cases",y="Probability")+
  theme_minimal()+
  theme(axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        axis.line.y.left = element_line())+
  facet_rep_grid(fct_rev(ifelse(trunc_by_onset,"Symptomatic self-isolation","No symptomatic self-isolation"))~time_period, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "feb_mar"="Lockdown (Comix Feb/Mar 2021)"))) + 
  scale_x_discrete(breaks=c(as.character(seq(0,18,by=2)),"\u2265 20"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,0.7),expand=c(0,0))

ggsave("results/contacts_over_time.png",width=9,height=5,units="in",dpi=400)
