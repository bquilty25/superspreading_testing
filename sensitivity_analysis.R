### Sensitivity analysis ---
#unadjusted BBC Pandemic tails
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
  #filter(period%in%c("Pre-pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  mutate(scenario_id=row_number()) %>% 
  select(-c(date_start,date_end)) %>% 
  crossing(heterogen_contacts=c(T))

processed_infections_adjusted <- run_model(testing_scenarios=testing_scenarios,contact_dat = contact_data,scenarios = time_periods_of_interest,browsing = F)

# R and K estimates over time
processed_infections_adjusted %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>%
  summarise.(dists=list(fitdist(sum_inf,"nbinom")),
             dist_means=list(fitdist(sum_inf,"nbinom")$estimate %>% enframe() %>% pivot_wider(names_from=name,values_from = value)),
             n=n(),
             ss_10=sum(sum_inf>10),
             ss_20=sum(sum_inf>20),
             ss_0=sum(sum_inf<=0),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim)) %>% 
  mutate.(
    prop_ss_10=ss_10/n*100,
    #prop_ss_20 =ss_20/n,
    prop_ss_0=ss_0/n*100) %>% 
  unnest.(dist_means) %>% 
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  filter.(variant=="wild") %>% 
  ggplot(aes(y=value,x=period,colour=name,group=name))+
  geom_point()+
  geom_line()+
  geom_text_repel(data=. %>% filter.(period=="Pre-pandemic",name=="mu"),
                  aes(x=period,y=value,label=paste0("R0 = ",round(value,1))),family="Lato",nudge_x=1)+
  geom_segment(data=rt_by_time_period %>% mutate(name="mu"),aes(x=period,xend=period,y=lower,yend=upper),alpha=0.25,size=10)+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = quad_col_pal,guide="none")+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  facet_grid2(name~.,switch="y",scales="free_y",independent = "y",
              labeller=labeller(name=c("mu"="Mean R","size"="k of R",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)")),
              axes="all",
              remove_labels = "x")+
  lims(y=c(0,NA))+
  labs(y="",
       x="Time period")+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("results/R and k over time adjusted.png",width=150,height=150,dpi=600,units="mm",bg="white")
