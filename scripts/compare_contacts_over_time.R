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
         e_other=rowSums(across(c(e_work,e_school,e_other),na.rm = T)),
         e_all = rowSums(across(c(e_home,e_other)),na.rm = T)) %>% 
  select(-c(e_work,e_school))

contact_data_dists <- contact_data %>%
  pivot_longer(c(e_all,e_home,e_other)) %>% 
  mutate(name=fct_relevel(name,"e_all","e_home","e_other")) %>% 
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

traj <- make_trajectories(n_cases = 500, n_sims = 200,variant=c("delta","wild"))

if(file.exists("results_traj.fst")){
  
  traj_ <- read_fst("results_traj.fst")
  
  } else {
    
  traj_ <- traj %>% 
    arrange(sim) %>% 
    group_split.(sim) %>% 
    map.(.f=inf_and_test,sampling_freq=NA) %>% 
    bind_rows.() 

write_fst(traj_,"results_traj.fst")

}

if(file.exists("results_ss.fst")){
  res <- fst::read_fst("results_ss.fst") %>% 
    crossing(prop_self_iso_symp = c(0,0.25,0.5,0.75,1),
             prop_self_iso_test = c(1))
  } else {
    
res <- traj_ %>%
  filter(is.na(sampling_freq)) %>% 
crossing(prop_self_iso_symp = 1,#c(0,0.25,0.5,0.75,1),
         prop_self_iso_test = c(1)) %>%
group_split.(sampling_freq,prop_self_iso_symp,prop_self_iso_test) %>% 
map(~sec_case_gen(.x)) %>% 
bind_rows.()

 fst::write_fst(res, "results_ss.fst")
 
}

dists <- res %>%  
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing")) %>% 
  nest.(data=c(idx,
           contacts_repeated,
           n_repeated_infected,
          contacts_casual,
          n_total_contacts,
          n_casual_infected, 
          n_total_infected)) %>% 
  mutate.(dists=map(.x=data,~hush(fitdist(.x$n_total_infected ,"nbinom"))),
         params=map(dists,~c(R=.x$estimate[[2]],k=.x$estimate[[1]]))) %>% 
  unnest_longer(params,values_to = "value",indices_to = "param") %>% 
  group_by(variant,
           time_period,
           prop_self_iso_symp,
           prop_self_iso_test,
           sampling_freq,
           param) %>% 
  summarise(x = quantile(value,probs=c(0.5,0.025,0.975)),q=c(0.5,0.025,0.975)) 

dists %>% pivot_wider(names_from = q,values_from = x) %>% 
  ggplot(aes(x=time_period,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_pointrange()+
  facet_grid(param~sampling_freq+time_period+variant,scales="free")+
  scale_y_log10()

dists %>%
  group_by(variant,
           time_period,
           sampling_freq,
           param,
           q) %>% 
  nest() %>% 
 mutate.(mod = map(.x = data, ~approxfun(y=.x$x,x=.x$prop_self_iso_symp)),
         pred= map(.x=mod, 
                           ~data.frame(prop_self_iso_symp = seq(0, 1, length.out = 101)) %>% 
                                   mutate(fitted=.x(prop_self_iso_symp)))) %>% 
  unnest.(pred) %>% mutate(q=as.character(q)) %>% pivot_wider(names_from = q,values_from = fitted) %>% 
  ggplot(aes(x=prop_self_iso_symp,y=`0.5`))+
  geom_line()+
  geom_ribbon(aes(x=prop_self_iso_symp,ymax=`0.975`,ymin=`0.025`),alpha=0.5)+
  #scale_y_log10()+
  facet_grid(param~sampling_freq+time_period+variant,scales="free")

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

