#Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from Comix
source("scripts/utils.R")
source("scripts/lee_infectivity.R")

#Load contact data
contacts_bbc_o18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_o18.csv")) 

contacts_bbc_u18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_u18.csv")) 

contacts_bbc <- bind_rows(contacts_bbc_o18, contacts_bbc_u18) %>%
  mutate(time_period = "pre",
         e_school=0)

contacts_comix_o18 <-
  read.csv(here("data", "comix_contact_distributions_o18.csv")) %>%
  mutate(
    date = as.Date(date),
    age="adults"
  )

contacts_comix_u18 <-
  read.csv(here("data", "comix_contact_distributions_u18.csv")) %>%
  mutate(
    date = as.Date(date),
    age="children"
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
contact_data <- comix_high_low %>% 
  bind_rows(contacts_bbc)%>% 
  mutate(time_period=fct_relevel(time_period,"pre","aug_sept","jan_feb"),
         e_work_school=rowSums(across(c(e_work,e_school),na.rm = T)),
         e_all = rowSums(across(c(e_home,e_work_school,e_other)),na.rm = T)) %>% 
  select(-c(e_work,e_school))

contact_data_dists <- contact_data %>%
  pivot_longer(c(e_all,e_home,e_work_school,e_other)) %>% 
  mutate(name=fct_relevel(name,"e_all","e_home","e_work_school","e_other")) %>% 
  drop_na(value) %>% 
  group_by(time_period,name) %>%
  nest() %>%
  mutate(dists=map(.x=data,.f= . %>% pull(value) %>% fitdist("nbinom")),
       params=map(dists,~c(mu=.x$estimate[[2]],k=.x$estimate[[1]]))) %>% 
  unnest_wider(params) 

contact_data %>% 
  pivot_longer(c(e_all,e_home,e_work_school,e_other)) %>% 
  drop_na(value) %>% 
  mutate(name=fct_relevel(name,"e_all","e_home","e_work_school","e_other"),
         value=factor(ifelse(value>=20,"\u2265 20",value)),
         value=fct_relevel(value,c(as.character(seq(0,19)),"\u2265 20"))) %>%
  ggplot(aes(x = value,
           y = ..prop..,
           group = 1))+
  geom_bar(width = 0.8)+
  geom_label(data = contact_data_dists,
             aes(label=paste0("Mean = ", sprintf("%.2f",mu),"\nk = ",sprintf("%.2f",k)),
                 x="\u2265 20",
                 y=0.7),
             colour = "white",
             fill = rgb(0,0,0,0.4),
             hjust = 1)+
  scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
  labs(x="Number of daily contacts",y="Probability")+
  theme_minimal()+
  theme(axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        axis.line.y.left = element_line())+
  facet_rep_grid(name~time_period, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                   name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                   prop_self_iso=function(x){scales::percent(as.numeric(x))})) + 
  scale_x_discrete(breaks=c(as.character(seq(0,18,by=2)),"\u2265 20"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))

ggsave("results/contacts.png",width=12,height=7,units="in",dpi=400,scale=0.8)

traj <- make_trajectories(n_cases = 500, n_sims = 200,variant="delta")

inf_and_test <- function(traj,sampling_freq=c(NA,3)){
  #browser()
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), traj$sim[1]))
  
  traj %>% as.data.frame() %>% 
  mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m,start=start,end=end)))  %>% 
  unnest_wider(infectiousness) %>% 
  ungroup() %>%
  mutate.(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf))) %>% 
  #testing
  crossing(sampling_freq = sampling_freq) %>% 
  mutate.(test_times = pmap(
    .f = test_times,
    list(
      sampling_freq = sampling_freq,
      onset_t = onset_t,
      type = type
    )
  )) %>%
  unnest.(test_times,.drop=F) %>%
  mutate.(
    ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
    test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
    test_label = detector(test_p = test_p,  u = runif(n = n(), 0, 1))
  ) %>%
  nest(ct, test_t, test_no, test_p, test_label) %>%
  mutate.(earliest_positive = map(.f = earliest_pos, .x = data)) %>%
  unnest.(earliest_positive,.drop=F) %>%
  select.(-data)
} 

traj_ <- traj %>% 
  arrange(sim) %>% 
  group_split.(sim) %>% 
  map.(.f=inf_and_test) %>% 
  flatten() %>% 
  bind_rows.() %>% 
crossing(prop_self_iso_symp = c(0,0.25,0.5,0.75,1),
            prop_self_iso_test = c(1)) 

fst::write.fst(traj_,"results_inf_curve.fst")

sec_case_gen <- function(df){
  
  message(sprintf("\n%s", Sys.time()))
  
  df %>%
  mutate.(self_iso_symp=ifelse(type=="symptomatic",rbinom(n=n(),size=1,prob=prop_self_iso_symp),0),
            self_iso_test=rbinom(n=n(),size=1,prob=prop_self_iso_test),
            test_t = ifelse(self_iso_test==0,Inf,test_t)) %>% 
  select.(-u) %>%
  crossing(time_period = factor(
    c("pre", "aug_sept", "jan_feb"),
    ordered = T,
    levels = c("pre", "aug_sept", "jan_feb")
  )) %>%
  mutate.(
    contacts_repeated_home = case_when.(
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
    contacts_repeated_work_school = case_when(
      time_period == "pre" ~ sample(
        contact_data %>%
          filter(time_period == "pre") %>%
          pull(e_work_school),
        size = n(),
        replace = T
      ),
      time_period == "aug_sept" ~ sample(
        contact_data %>%
          filter(time_period == "aug_sept") %>%
          pull(e_work_school),
        size = n(),
        replace = T
      ),
      time_period == "jan_feb" ~ sample(
        contact_data %>%
          filter(time_period == "jan_feb") %>%
          pull(e_work_school),
        size = n(),
        replace = T
      )
    ),
    trunc_t=case_when.(
      # if symptomatic, adhering to self isolation, and either not tested or test neg,
      # truncate at onset
      type=="symptomatic"&is.infinite(test_t)&self_iso_symp!=0~onset_t,
      # if symptomatic, adhering to self isolation, and have onset before test, truncate at onset
      type=="symptomatic"&is.finite(test_t)&onset_t<test_t&self_iso_symp!=0~onset_t,
      # if symptomatic, adhering to self isolation, and have onset after pos test, truncate at test
      type=="symptomatic"&is.finite(test_t)&test_t<onset_t&self_iso_symp!=0~test_t,
      # if symptomatic, not adhering to self isolation, and have a positive test, truncate at test
      TRUE ~ test_t), 
    #if symp onset, contacts outside of home should cease; for school/work, multiply 
    #number of contacts by proportion of time since infection which occurs pre-onset
    contacts_repeated_work_school=case_when.(!is.infinite(trunc_t)~ceiling(contacts_repeated_work_school*(trunc_t/(end-start))),
                                            TRUE~ceiling(contacts_repeated_work_school))
  ) %>% 
  mutate.(
    contacts_repeated   = contacts_repeated_home + contacts_repeated_work_school,
    n_repeated_infected = rbinom(n = n(), 
                                      size = contacts_repeated, 
                                      prob = norm_sum)) %>%
  #Daily casual contacts
  unnest.(infectiousness) %>%
  mutate.(norm_daily = (culture / sum_inf) * norm_sum) %>%
  mutate.(
    contacts_casual = case_when.(
      t>trunc_t  ~ 0L,
      time_period == "pre" ~ sample(
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
  mutate.(n_casual_infected = rbinom(n = n(), 
                                    size = contacts_casual, 
                                    prob =
                                      norm_daily)) %>%
  summarise.(.by=c(sim,
                   idx,
                   type,
                   time_period,
                   norm_sum,
                   prop_self_iso_symp,
                   prop_self_iso_test,
                   contacts_repeated,
                   n_repeated_infected,
                   sampling_freq),
    contacts_casual = sum(contacts_casual),
    n_casual_infected = sum(n_casual_infected)) %>%
  mutate.(n_total_infected = n_repeated_infected + n_casual_infected) 
}

res <- inf_curve %>% 
  group_split.(sampling_freq,prop_self_iso_symp,prop_self_iso_test) %>% 
  map(~sec_case_gen(.x))

fst::write.fst(res %>% bind_rows.(),path = "results_ss.fst")

dists <- res %>% bind_rows.() %>%  
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing")) %>% 
  nest.(data=c(idx,
           type,
           norm_sum,
           contacts_repeated,
           n_repeated_infected,
          contacts_casual,
          n_casual_infected, 
          n_total_infected)) %>% 
  mutate.(dists=map(.x=data,.f= .%>% 
                     pull(n_total_infected) %>% 
                     fitdist("nbinom")),
         params=map(dists,~c(R=.x$estimate[[2]],k=.x$estimate[[1]]))) %>% 
  unnest_longer(params,values_to = "value",indices_to = "param") %>% 
  group_by(sim,time_period,
           #prop_self_iso_symp,
           #prop_self_iso_test,
           sampling_freq,
           param) %>% 
  nest() %>% 
  mutate.(mod = map(.x = data, ~approxfun(.x$value~.x$prop_self_iso_symp))) 
  #mutate.(q = map(.x = data, ~quantile(.x$value,probs=c(0.5,0.025,0.975)))) %>% 
  #unnest_wider(q) 

dists %>% 
  mutate.(pred=map(mod,~broom::augment(.x,newdata=newdata1,type.predict="response"))) %>% 
  unnest.(pred) %>% 
  nest.(data=c(sim,.fitted,.se.fit)) %>% 
  mutate.(q = map(.x = data, ~quantile(.x$.fitted,probs=c(0.5,0.025,0.975)))) %>%  
  unnest_wider(q) %>% 
  #filter.(param=="R",time_period!="jan_feb") %>% 
  ggplot(aes(x=prop_self_iso_symp,y=`50%`))+
  geom_line()+
  geom_ribbon(aes(x=prop_self_iso_symp,ymax=`97.5%`,ymin=`2.5%`),alpha=0.5)+
  #scale_y_log10()+
  facet_grid(param~sampling_freq+time_period,scales="free")

ggsave("results/smooths.png",width=210,height=120,dpi=600,units="mm")

plot_labels <- dists %>% 
  ungroup() %>% 
  select(-c(data)) %>% 
  pivot_wider(names_from=param,values_from=c(`50%`,`2.5%`,`97.5%`),names_glue = "{param}_{.value}") %>% 
  mutate(R_estimate=paste0("R<sub>c</sub> = ",sprintf("%.2f",`R_50%`), " (",sprintf("%.2f",`R_2.5%`)," - ",sprintf("%.2f",`R_97.5%`),")"),
         k_estimate = paste0("k<sub>c</sub> = ",sprintf("%.2f",`k_50%`), " (",sprintf("%.2f",`k_2.5%`)," - ",sprintf("%.2f",`k_97.5%`),")")) %>% 
  select(-c(`R_50%`, `k_50%`, `R_2.5%`, `k_2.5%`, `R_97.5%`, `k_97.5%`)) %>% 
  arrange(prop_self_iso_symp,prop_self_iso_test,time_period,sampling_freq) %>% 
  select(prop_self_iso_symp,prop_self_iso_test,time_period,sampling_freq,everything()) 

p1 <- res %>% bind_rows.() %>%  
  filter.(prop_self_iso_symp%in%c(0,0.5,1)) %>% 
  mutate.(sampling_freq=ifelse(!is.na(sampling_freq),"Testing every 3 days","No testing")) %>% 
  mutate.(n_total_infected=factor(ifelse(n_total_infected>=10,"\u2265 10",as.character(n_total_infected))),
         n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
  #as_tibble() %>% 
  ggplot(aes(x = n_total_infected,
             y = ..prop..,
             fill = factor(sampling_freq),
             group = 1))+
  geom_bar(width = 0.8)+
  geom_richtext(data = plot_labels %>% filter.(prop_self_iso_symp%in%c(0,0.5,1)),
             aes(label=paste0(R_estimate,"<br>",k_estimate),
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
  facet_rep_grid(prop_self_iso_symp~time_period+sampling_freq, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                   prop_self_iso_symp=function(x){scales::percent(as.numeric(x))})) + 
  scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))

ggsave(plot = p1,"results/sec_cases_self_iso_0_1.png",width=12,height=5,units="in",dpi=400,scale=1.15)

s1 <- res %>% bind_rows.() %>%  
  mutate.(sampling_freq=ifelse(!is.na(sampling_freq),"Testing every 3 days","No testing")) %>% 
  mutate.(n_total_infected=factor(ifelse(n_total_infected>=10,"\u2265 10",as.character(n_total_infected))),
          n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
  ggplot(aes(x = n_total_infected,
             y = ..prop..,
             fill = factor(sampling_freq),
             group = 1))+
  geom_bar(width = 0.8)+
  geom_richtext(data = plot_labels,
             aes(label=paste0(R_estimate,"<br>",k_estimate),
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
  facet_rep_grid(prop_self_iso_symp~time_period+sampling_freq, scales='free_y',
                 labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                 "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                 "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                   prop_self_iso_symp=function(x){scales::percent(as.numeric(x))})) +  
  scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))

ggsave(plot = s1,"results/sec_cases_prop_self_iso_full.png",width=12,height=9,units="in",dpi=400,scale=1.15)

#repeated/casual split
res %>% bind_rows.()  %>% 
  group_by(time_period,
           prop_self_iso,
           sampling_freq) %>% 
  summarise(prop_repeated = sum(n_repeated_infected) / sum(n_total_infected),
            prop_casual = sum(n_casual_infected) / sum(n_total_infected)) %>% 
  base::print(n=Inf)

res %>% bind_rows.() %>%  
  filter(time_period=="pre",
         prop_self_iso%in%c(0,1)) %>% 
  group_by(time_period,
           prop_self_iso,
           sampling_freq) %>%
  mutate(plus_10=n_total_infected>=10) %>% 
  count(plus_10) %>% 
  mutate(prop = prop.table(n)) %>% 
  group_by(time_period,
           prop_self_iso,
           sampling_freq) #%>% 
  # nest() %>% 
  # ungroup() %>% 
  # mutate(est = map(.x = data, ~quantile(1-.$prop,
  #                                       probs = c(0.5,0.025,0.975)))) %>% 
  # unnest_wider(est) 

