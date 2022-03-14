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
    
scenarios <- crossing(time_periods) %>% 
  filter(period%in%c("BBC Pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end))

traj_scenarios <- crossing(traj,
                 scenarios,
                 heterogen_contacts=c(TRUE,FALSE)) %>% 
  mutate.(repeated_contacts = case_when.(heterogen_contacts~sample(contact_data$e_home[contact_data$period==period],
                                     size=n(),
                                     replace=T),
                                     TRUE~round(mean(contact_data$e_home[contact_data$period==period]))),
                                     .by=c(period)
                                     )

key_grouping_var <- c(sim,variant,scenario_id,period,lower_inf_thresh,heterogen_vl,heterogen_contacts)
  
  trunc_t=case_when.(
  # if symptomatic, adhering to self isolation, and either not tested or test neg,
  # truncate at onset
  symptomatic & is.infinite(test_t) & self_iso_symp != 0 ~ onset_t,
  # if symptomatic, adhering to self isolation, and have onset before test, truncate at onset
  symptomatic &
    is.finite(test_t) & onset_t < test_t & self_iso_symp != 0 ~ onset_t,
  # if symptomatic, adhering to self isolation, and have onset after pos test, truncate at test
  symptomatic &
    is.finite(test_t) & test_t < onset_t & self_iso_symp != 0 ~ test_t,
  # if symptomatic, not adhering to self isolation, and have a positive test, truncate at test
  TRUE ~ test_t)

#calculate repeated infections
traj_scenarios_joined <- traj_ %>% 
  left_join.(traj_scenarios)  

repeated_infections <- traj_scenarios_joined %>% 
  uncount.(repeated_contacts,.id="id",.remove = F) %>% 
  mutate.(hh_duration = sample(contacts_hh_duration,size=n(),replace=T),
          infected    = rbernoulli(n(),p=culture_p*hh_duration)) %>% 
  filter.(infected==T) %>% 
  slice.(min(t), .by=c(sim,variant,scenario_id,period,lower_inf_thresh,heterogen_vl,heterogen_contacts,repeated_contacts,id)) %>% 
  count.(sim,variant,scenario_id,period,lower_inf_thresh,heterogen_vl,heterogen_contacts,t,repeated_contacts,name = "repeated_infected")

#calculate casual infections
casual_infections <- traj_scenarios_joined %>% 
  mutate.(casual_contacts = case_when.(heterogen_contacts~ sample(contact_data$e_other[contact_data$period==period],size=n(),replace=T),
                                       TRUE ~ round(mean(contact_data$e_other[contact_data$period==period]))),.by=period) %>%
  #mutate.(contacts_casual=ifelse(t>=trunc_t,0L, contacts_casual)) %>% 
  uncount.(casual_contacts,.remove=F) %>% 
  mutate.(nhh_duration=sample(contacts_nhh_duration,size=n(),replace=T),
          infected = rbernoulli(n=n(),p = culture_p*nhh_duration)) %>% 
  #filter.(infected==T) %>% 
  count.(t,sim,variant,scenario_id,period,lower_inf_thresh,heterogen_vl,heterogen_contacts,infected) %>% 
  pivot_wider.(values_from=N,names_from=infected,values_fill = 0) %>% 
  mutate.(casual_contacts=`FALSE`+`TRUE`) %>% 
  select.(everything(),"casual_infected"=`TRUE`,-`FALSE`)# %>% 
  mutate.(casual_infected=ifelse(t>5,0,casual_infected))
  
#testing


processed_infections <- traj_scenarios_joined %>% 
  left_join.(casual_infections) %>% 
  left_join.(repeated_infections) %>% 
  replace_na.(list(repeated_infected=0,casual_infected=0,casual_contacts=0)) %>% 
  arrange.(period,lower_inf_thresh) %>% 
  mutate.(
    total_contacts = casual_contacts+repeated_contacts,
    total_infections=casual_infected+repeated_infected)
  
processed_infections %>% filter.(lower_inf_thresh==F) %>% ggplot(aes(x=vl,y=total_contacts,colour=total_infections))+geom_jitter(alpha=0.5)+
  scale_y_log10()+
  scale_colour_viridis_c(trans="log10",na.value=NA,option ="turbo")+facet_grid(heterogen_vl~heterogen_contacts+period)

ggsave("contacts_infections.png")

boot_est <- processed_infections %>% 
  filter.(lower_inf_thresh==F) %>% 
  summarise.(sum_inf=sum(total_infections),.by=c(sim,variant,lower_inf_thresh,scenario_id,period,heterogen_vl,heterogen_contacts)) %>% 
  summarise.(dists=list(fitdist(sum_inf,"nbinom")$estimate %>% t()),.by=c(variant,lower_inf_thresh,scenario_id,period,heterogen_vl,heterogen_contacts)) #%>%
  # mutate.(boot_dist=map.(.x=dists, ~bootdist(f =.,bootmethod = "nonparam",parallel="snow",ncpus=8)$CI %>% 
  #                          as.data.frame() %>% rownames_to_column())) 
 
boot_est %>% unnest.(dists) %>% pivot_longer.(c(size,mu)) %>% ggplot(aes(y=value,x=period))+geom_point()+facet_grid(heterogen_vl+name~heterogen_contacts,scales="free_y",labeller=label_both)

boot_est %>% unnest.(boot_dist) %>% 
  ggplot()+
  geom_pointrange(aes(y=Median,ymin=`2.5%`,ymax=`97.5%`,x=factor(period)))+facet_wrap(~rowname,scales="free_y",ncol = 1)+lims(y=c(0,NA))

ggsave("mu_size.png")

#testing


traj_[t>=begin_testing, 
      start_iso := fifelse(any(test_label), t[which.max(test_label)], Inf),
      by=c(sim)]

sum_res <- res %>%
  summarise.(.by = c(sim,
                   idx,
                   type,
                   variant,
                   time_period,
                   prop_self_iso_symp,
                   prop_self_iso_test,
                   contacts_repeated,
                   sampling_freq),
             contacts_casual     = sum(contacts_casual),
             n_casual_infected   = sum(n_casual_infected),
             n_repeated_infected = sum(n_repeated_infected)) %>%
  mutate.(n_total_contacts       = contacts_repeated + contacts_casual,
          n_total_infected       = n_repeated_infected + n_casual_infected)

dists <- sum_res %>%  
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                 TRUE~ "No testing")) %>% 
  nest.(data=c(idx,
           contacts_repeated,
           n_repeated_infected,
          contacts_casual,
          n_total_contacts,
          n_casual_infected, 
          n_total_infected)) %>% 
  mutate.(dists=map(.x=data,~hush(fitdist(.x$n_total_infected,"nbinom"))),
         params=map(dists,~c(R=.x$estimate[[2]],k=.x$estimate[[1]]))) %>% 
  unnest_longer(params,values_to = "value",indices_to = "param") %>% 
  group_by(variant,
           time_period,
           prop_self_iso_symp,
           prop_self_iso_test,
           sampling_freq,
           param) %>% 
  summarise(x = quantile(value,probs=c(0.5,0.025,0.975)),q=c(0.5,0.025,0.975)) %>% 
  pivot_wider(names_from = q,values_from = x)

dists %>% 
  filter(prop_self_iso_symp==1,prop_self_iso_test%!in%c(0,0.5)) %>% 
  ggplot(aes(x=sampling_freq,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_hline(aes(yintercept = ifelse(param=="R",1,NA)),linetype="dashed")+
  facet_grid(fct_rev(param)~time_period,#,scales="free",
             labeller=labeller(time_period=c("BBC Pandemic"="Pre-pandemic (BBC 2017-18)",
                                             "Relaxed restrictions"="Relaxed (Comix Aug 2020)",
                                             "School reopening"="Schools reopening (Comix Sept 2020)",
                                             "Lockdown 1"="Lockdown (Comix Apr-Jul 2021)"),
                               name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                               prop_self_iso=function(x){scales::percent(as.numeric(x))}))+
  geom_linerange(aes(ymin  = `0.025`, 
                     ymax  = `0.975`,
                     colour=param,
                     group=time_period),
                 position  = position_dodge2(width = 1),
                 alpha     = 0.4,
                 size      = 3) +
  geom_point(#pch           = "-",
             size          = 2,
             position      = position_dodge2(width = 1),
             aes(y         = `0.5`,
                 colour    = param,
                 group=time_period)) +
  scale_color_manual(name = "",guide=F,
                     values = covid_pal)+
  labs(x="Asymptomatic testing frequency",
       y="Estimate")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "bottom",
        strip.placement = "outside")

ggsave("results/R_k.png",width=250,height=150,dpi=600,units="mm")

#prop infecting 0, >10, >20...
ss_dat <- sum_res %>% 
  filter(prop_self_iso_test==1,prop_self_iso_symp==1) %>%  
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                   TRUE~ "No testing")) %>% 
  group_by(idx,variant,
           time_period,
           prop_self_iso_symp,
           prop_self_iso_test,
           sampling_freq) %>% 
  summarise(n=n(),
            ss_10=sum(n_total_infected>=10),
            prop_ss_10=ss_10/n,
            ss_0=sum(n_total_infected<=0),
            prop_ss_0=ss_0/n) %>% 
  pivot_longer.(cols=c(prop_ss_0,prop_ss_10),names_to="param") %>%  
  group_by(variant, 
           time_period,
           prop_self_iso_symp, 
           prop_self_iso_test,
           sampling_freq,
           param) %>% 
  summarise(x = quantile(value,probs=c(0.5,0.025,0.975),na.rm = T),q=c(0.5,0.025,0.975)) %>% 
  pivot_wider(names_from = q,values_from = x)

ss_0_plot <- ss_dat %>% 
    filter(param=="prop_ss_0") %>% 
    ggplot(aes(x=sampling_freq,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_linerange(aes(ymin  = `0.025`, 
                     ymax  = `0.975`,
                     colour=param,
                     group=time_period),
                 position  = position_dodge2(width = 1),
                 alpha     = 0.4,
                 size      = 3) +
  geom_point(#pch           = "-",
    size          = 2,
    position      = position_dodge2(width = 1),
    aes(y         = `0.5`,
        colour    = param,
        group=time_period)) 

ss_20_plot <- ss_dat %>% 
  filter(param=="prop_ss_10") %>% 
  ggplot(aes(x=sampling_freq,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_linerange(aes(ymin  = `0.025`, 
                     ymax  = `0.975`,
                     colour=param,
                     group=time_period),
                 position  = position_dodge2(width = 1),
                 alpha     = 0.4,
                 size      = 3) +
  geom_point(#pch           = "-",
    size          = 2,
    position      = position_dodge2(width = 1),
    aes(y         = `0.5`,
        colour    = param,
        group=time_period)) +
  scale_y_log10()
    
  
(ss_0_plot/ss_20_plot)&scale_color_manual(name = "",guide=F,
                     values = covid_pal)&
  labs(x="Asymptomatic testing frequency",
       y="Proportion infecting x")&
  theme_minimal()&
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "bottom",
        strip.placement = "outside")&
  facet_grid(param~time_period,scales="free",
             labeller=labeller(param=c(prop_ss_0="0 secondary infections",
                                       prop_ss_10=">= 10 secondary infections"),
                               time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                             "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                             "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                               name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                               prop_self_iso=function(x){scales::percent(as.numeric(x))}))

ggsave("results/prop_ss.png",width=210,height=150,dpi=600,units="mm")

# SAR
sum_res %>%  
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                   TRUE~ "No testing"),
          hh_SAR=n_repeated_infected/contacts_repeated,
          nhh_SAR=n_casual_infected/contacts_casual,
          tot_SAR=n_total_infected/n_total_contacts) %>% 
  pivot_longer(cols = c("hh_SAR","nhh_SAR","tot_SAR"),values_to = "value",names_to =  "param") %>% 
  replace_na(list(value=0)) %>% 
  group_by(variant,
           time_period,
           prop_self_iso_symp,
           prop_self_iso_test,
           sampling_freq,
           param) %>% 
  summarise(x = quantile(value,probs=c(0.5,0.025,0.975),na.rm = T),q=c(0.5,0.025,0.975)) %>% 
  pivot_wider(names_from = q,values_from = x) %>% 
  ggplot(aes(x=time_period,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_pointrange()+
  facet_grid(param~sampling_freq+time_period+variant,scales="free")

#Heatmap
heatmap <- res %>% 
  filter( prop_self_iso_symp==1,prop_self_iso_test==1,is.na(sampling_freq)) %>%   
  filter(n_casual_infected>0,n_repeated_infected>0) %>%  
  ggplot(aes(x=ct,
             y=contacts_repeated+contacts_casual,
             colour=n_repeated_infected+n_casual_infected))+
  geom_point()+
  scale_y_log10()+
  scale_colour_viridis_c(trans="log",option = "turbo",breaks=c(1,10,100),limits=c(1,250), guide = guide_colourbar(nbin = 100, draw.ulim = TRUE, draw.llim = TRUE))+
  labs(x="Viral load (CT)",y="Exposed contacts",colour="Infected contacts")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "right",
        strip.placement = "outside")+
  facet_grid(~time_period,labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                          "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                          "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                            name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                            prop_self_iso=function(x){scales::percent(as.numeric(x))}))

ggsave("results/heatmap.png",width=210,height=150,dpi=600,units="mm")

hex_map <- res %>% 
  filter( prop_self_iso_symp==1,prop_self_iso_test==1,is.na(sampling_freq)) %>%   
  filter(n_casual_infected>0,n_repeated_infected>0) %>%  
  ggplot(aes(x=ct,
             y=contacts_repeated+contacts_casual))+
  geom_hex(aes(fill=stat(density)))+
  scale_y_log10()+
  scale_fill_viridis_c(option = "rocket", labels=scales::percent_format(accuracy=1L), guide = guide_colourbar(nbin = 100, draw.ulim = F, draw.llim = F))+
  labs(x="Viral load (CT)",y="Exposed contacts",fill="Density")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "right",
        strip.placement = "outside")+
  facet_grid(~time_period,labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                          "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                          "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                            name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                            prop_self_iso=function(x){scales::percent(as.numeric(x))}))

heatmap/hex_map+plot_layout(guides="collect")
ggsave("results/heatmap_hexmap.png",width=300,height=200,dpi=600,units="mm")

# Generation time
res %>% 
  filter( prop_self_iso_symp==1,prop_self_iso_test==1,is.na(sampling_freq)) %>%   
  filter(n_casual_infected>0,n_repeated_infected>0) %>%  
  pivot_longer(cols=c(n_repeated_infected,n_casual_infected)) %>% 
  ggplot(aes(x=t,
             y=value,
             fill=name))+
  geom_col()+
  scale_fill_brewer(type="qual")+
  labs(x="Viral load (CT)",y="Exposed contacts",colour="Infected contacts")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "right",
        strip.placement = "outside")+
  facet_grid(type~time_period,labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                          "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                          "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                            name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                            prop_self_iso=function(x){scales::percent(as.numeric(x))}))

# Proportion asymptomatic
res %>% 
  filter(prop_self_iso_test==1,is.na(sampling_freq)) %>%   
  mutate(onset_t=case_when(type=="asymptomatic"~Inf,
                        TRUE~onset_t),
         asymp=t<onset_t) %>%  
  filter(n_casual_infected>0,n_repeated_infected>0) %>%  
  pivot_longer(cols=c(n_repeated_infected,n_casual_infected)) %>% 
  ggplot(aes(x=t,
             y=value,
             fill=asymp))+
  geom_col()+
  scale_fill_brewer(type="qual")+
  labs(x="Viral load (CT)",y="Exposed contacts",colour="Infected contacts")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "right",
        strip.placement = "outside")+
  facet_grid(type~time_period+prop_self_iso_symp,labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                              "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                              "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                                name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                                prop_self_iso=function(x){scales::percent(as.numeric(x))}))

#  Deprecated 
# 
# plot_labels <- dists %>% 
#   ungroup() %>% 
#   select(-c(data)) %>% 
#   pivot_wider(names_from=param,values_from=c(`50%`,`2.5%`,`97.5%`),names_glue = "{param}_{.value}") %>% 
#   mutate(R_estimate=paste0("R<sub>c</sub> = ",sprintf("%.2f",`R_50%`), " (",sprintf("%.2f",`R_2.5%`)," - ",sprintf("%.2f",`R_97.5%`),")"),
#          k_estimate = paste0("k<sub>c</sub> = ",sprintf("%.2f",`k_50%`), " (",sprintf("%.2f",`k_2.5%`)," - ",sprintf("%.2f",`k_97.5%`),")")) %>% 
#   select(-c(`R_50%`, `k_50%`, `R_2.5%`, `k_2.5%`, `R_97.5%`, `k_97.5%`)) %>% 
#   arrange(prop_self_iso_symp,prop_self_iso_test,time_period,sampling_freq) %>% 
#   select(prop_self_iso_symp,prop_self_iso_test,time_period,sampling_freq,everything()) 
# 
# p1 <- res %>% bind_rows.() %>%  
#   filter.(prop_self_iso_symp%in%c(0,0.5,1)) %>% 
#   mutate.(sampling_freq=ifelse(!is.na(sampling_freq),"Testing every 3 days","No testing")) %>% 
#   mutate.(n_total_infected=factor(ifelse(n_total_infected>=10,"\u2265 10",as.character(n_total_infected))),
#          n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
#   #as_tibble() %>% 
#   ggplot(aes(x = n_total_infected,
#              y = ..prop..,
#              fill = factor(sampling_freq),
#              group = 1))+
#   geom_bar(width = 0.8)+
#   geom_richtext(data = plot_labels %>% filter.(prop_self_iso_symp%in%c(0,0.5,1)),
#              aes(label=paste0(R_estimate,"<br>",k_estimate),
#                  x="\u2265 10",
#                  y=0.75),
#              colour = "white",
#              fill = rgb(0,0,0,0.4),
#              hjust = 1)+
#   scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
#   labs(x="Number of secondary cases",y="Probability")+
#   theme_minimal()+
#   theme(axis.line.x.bottom = element_line(),
#         axis.ticks.x.bottom = element_line(),
#         axis.line.y.left = element_line())+
#   facet_rep_grid(prop_self_iso_symp~time_period+sampling_freq, scales='free_y',
#                  labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
#                                                  "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
#                                                  "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
#                                    prop_self_iso_symp=function(x){scales::percent(as.numeric(x))})) + 
#   scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
#   scale_y_continuous(limits=c(0,1),expand=c(0,0))
# 
# ggsave(plot = p1,"results/sec_cases_self_iso_0_1.png",width=12,height=5,units="in",dpi=400,scale=1.15)
# 
# s1 <- res %>% bind_rows.() %>%  
#   mutate.(sampling_freq=ifelse(!is.na(sampling_freq),"Testing every 3 days","No testing")) %>% 
#   mutate.(n_total_infected=factor(ifelse(n_total_infected>=10,"\u2265 10",as.character(n_total_infected))),
#           n_total_infected=fct_relevel(n_total_infected,c(as.character(seq(0,9)),"\u2265 10"))) %>%
#   ggplot(aes(x = n_total_infected,
#              y = ..prop..,
#              fill = factor(sampling_freq),
#              group = 1))+
#   geom_bar(width = 0.8)+
#   geom_richtext(data = plot_labels,
#              aes(label=paste0(R_estimate,"<br>",k_estimate),
#              x="\u2265 10",
#              y=0.75),
#              colour = "white",
#              fill = rgb(0,0,0,0.4),
#             hjust = 1)+
#   scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
#   labs(x="Number of secondary cases",y="Probability")+
#   theme_minimal()+
#   theme(axis.line.x.bottom = element_line(),
#         axis.ticks.x.bottom = element_line(),
#         axis.line.y.left = element_line())+
#   facet_rep_grid(prop_self_iso_symp~time_period+sampling_freq, scales='free_y',
#                  labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
#                                                  "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
#                                                  "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
#                                    prop_self_iso_symp=function(x){scales::percent(as.numeric(x))})) +  
#   scale_x_discrete(breaks=c(as.character(seq(0,8,by=2)),"\u2265 10"),expand = expansion(add=0.7))+
#   scale_y_continuous(limits=c(0,1),expand=c(0,0))
# 
# ggsave(plot = s1,"results/sec_cases_prop_self_iso_full.png",width=12,height=9,units="in",dpi=400,scale=1.15)
