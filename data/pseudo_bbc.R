contact_distributions_o18_full <- read_csv("data/contact_distributions_o18_full.csv") %>% rownames_to_column()

contact_distributions_o18_full %>% 
  mutate(e_all=e_home+e_work+e_other) %>% 
  summarise.(n=n(),over_250=sum(e_all>=250)) %>% mutate.(prop=over_250/n*100)
  
  ggplot(aes(x=value))+
  geom_histogram(binwidth = 1)+
  facet_wrap(~name)
