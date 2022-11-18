individuals <- processed_infections %>%
  filter.(sim%in%rdunif(10,1,N_sims),scenario_id==1,sampling_freq==7,prop_self_iso_test%in%c(0,1),variant=="wild") 

traj_plot <- individuals %>% distinct.(sim,heterogen_vl,lower_inf_thresh,variant) %>% left_join.(traj_) %>% ggplot(aes(x=t,y=vl,group=sim))+geom_line()+facet_grid(~sim,scales="free_y")

inf_plot <- individuals %>% 
  mutate.(infected_hh=hh_infected,notinfected_hh=hh_contacts-hh_infected,infected_nhh=nhh_infected,notinfected_nhh=nhh_contacts-nhh_infected) %>% 
  pivot_longer.(c(infected_hh,infected_nhh,notinfected_hh,notinfected_nhh)) %>% 
  uncount(value) %>% 
  separate.(name,into=c("infected","contact_type"),sep="_") %>% 
  #filter.(name=="infected")%>% 
  ggplot(aes(x=t,fill=infected,colour=contact_type))+
  geom_dotplot(binwidth = 1,stackgroups=TRUE,binpositions="all",dotsize=0.5)+
  scale_fill_manual(values=c(muted("red"),"grey"))+
  scale_colour_manual(values=c("black","white"))+
  scale_x_continuous(breaks=breaks_width(2))+
  scale_y_continuous(NULL, breaks = NULL)+
  facet_grid(prop_self_iso_test~sim,scales="free_y")

(traj_plot/inf_plot)&
  plotting_theme

