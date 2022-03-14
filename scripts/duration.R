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

