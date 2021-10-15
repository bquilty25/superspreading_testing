pacman::p_load(qs,tidyverse)

contacts_qs <- qs::qread("data/contacts.qs")

contacts_hh_duration <- contacts_qs %>% 
  filter(cnt_household==1) %>% 
  mutate(cnt_duration=cnt_minutes_max/1440) %>% 
  filter(cnt_duration<=1) %>% 
  pull(cnt_duration)


contacts_nhh_duration <- contacts_qs %>% 
  filter(cnt_household==0) %>% 
  mutate(cnt_sat=as.integer(cnt_minutes_max>60)) %>% 
  mutate(cnt_duration=cnt_minutes_max/1440) %>% 
  filter(cnt_duration<=1)%>% 
  pull(cnt_duration)

