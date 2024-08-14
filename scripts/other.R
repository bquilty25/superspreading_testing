#sort individuals by auc, determine proportion of transmission by top 20% of auc
processed_infections_baseline %>% 
  filter(period=="Pre-pandemic") %>% 
  summarise.(sum_inf=sum(total_infections),
             sum_culture = sum(culture_p),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>% 
  #filter(sum_inf>0) %>% 
  mutate.(rank_inf=ecdf(sum_inf)(sum_inf),
          rank_culture=ecdf(sum_culture)(sum_culture),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim)) %>% 
  ggplot(aes(x=rank_culture,y=rank_inf))+
  geom_point()+
  geom_smooth()+
  lims(x=c(0,1),y=c(0,1))+
  facet_wrap2(~period)

  
processed_infections_baseline %>% 
  filter(period=="Pre-pandemic") %>% 
  summarise.(sum_inf=sum(total_infections),
             sum_culture = sum(culture_p),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>% 
  filter(sum_inf>0) %>% 
  ggplot(aes(x=sum_inf,y=sum_culture))+
  geom_density_2d_filled()+
  #lims(x=c(0,1),y=c(0,1))+
  facet_wrap2(~period)
  