pacman::p_load(qs,tidyverse)

contacts_qs <- qs::qread("data/contacts.qs")

contacts_hh_duration <- contacts_qs %>% 
  drop_na(cnt_minutes_max) %>% 
  filter(cnt_household==1) %>% 
  mutate(cnt_duration=cnt_minutes_max/1440,cnt_duration=ifelse(cnt_duration>=1,1,cnt_duration)) %>% 
  pull(cnt_duration)



contacts_nhh_duration <- contacts_qs %>% 
  drop_na(cnt_minutes_max) %>% 
  filter(cnt_household==0) %>% 
  mutate(cnt_duration=cnt_minutes_max/1440,cnt_duration=ifelse(cnt_duration>=1,1,cnt_duration)) %>% 
  pull(cnt_duration)

contacts_qs %>% 
  drop_na.(cnt_minutes_max) %>% 
  mutate.(cnt_minutes_max=ifelse(cnt_minutes_max>=1440,1440,cnt_minutes_max)) %>% 
  ggplot()+
  geom_histogram(aes(x=cnt_minutes_max/60,y=..density..,fill=factor(cnt_household)),binwidth = 1)+
  scale_x_continuous("Per-contact time (hours)",breaks = breaks_width(2))+
  scale_y_continuous("Density")+
  scale_fill_met_d(name="Signac",override.order = FALSE,direction=-1,guide=F)+
  facet_wrap(~cnt_household,labeller=labeller(cnt_household=c(`1`="Household",`0`="Non-household")))+
  plotting_theme+
  theme(legend.position = "none")

contacts_qs %>% 
  drop_na.(cnt_minutes_max) %>% 
  mutate.(cnt_minutes_max=ifelse(cnt_minutes_max>=1440,1440,cnt_minutes_max)) %>% 
  summarise.(list(enframe(quantile(cnt_minutes_max,c(0.025, 0.5, 0.975)),"quantile","duration")),.by=cnt_household) %>% 
  unnest.(V1)

ggsave("results/duration_hist.png",width=210,height=100,dpi=600,units="mm",bg="white")
