
#### Results ####

#inferred generation time
processed_infections_baseline %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var[-1]),
                   sampling_freq,prop_self_iso_test,t)) %>% 
  ggplot(aes(x=t,y=(sum_inf)/sum(sum_inf)))+
  geom_col(fill=bi_col_pal[2])+
  stat_function(fun=dweibull,
                args=list(shape=2.826,scale=5.665)
                ,colour=bi_col_pal[1])+
  geom_text(aes(x=Inf,y=Inf,
                hjust=1,
                vjust=1,
                label="Ferretti et al.\nWeibull shape: 2.826\nWeibull scale: 5.665"),
            family="Lato",
            colour=bi_col_pal[1])+
  scale_x_continuous("Days since first detectable by PCR (Ct<40)",
                     breaks=breaks_width(5),
                     limits = c(0,20))+
  scale_y_continuous("Proportion",labels=label_percent())+
  plotting_theme

ggsave("results/gen_time.png",width=210,height=120,dpi=600,units="mm",bg="white")

# R and K estimates over time
boot_est <- processed_infections_baseline %>% 
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

ggsave("results/R and k over time.png",width=150,height=150,dpi=600,units="mm",bg="white")

# heterogen_on_off
heterogen_plot <- processed_infections_heterogen_on_off %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>%
  summarise.(.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim),
             dist_means=list(fitdist(sum_inf,"nbinom")$estimate %>% 
                               enframe() %>% 
                               pivot_wider(names_from=name,values_from = value)
                             ),
             n=n(),
             ss_10=sum(sum_inf>10),
             ss_20=sum(sum_inf>20),
             ss_0=sum(sum_inf<=0)) %>% 
  mutate.(
    prop_ss_10=ss_10/n*100,
    prop_ss_0=ss_0/n*100) %>% 
  unnest.(dist_means) %>% 
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  filter.(variant=="wild",!(heterogen_vl==F&heterogen_contacts==F&name=="size")) %>% 
  mutate.(heterogen_label=case_when.(heterogen_vl&heterogen_contacts~"Variable viral load and contacts",
                               heterogen_vl&!heterogen_contacts~"Variable viral load, equal contacts",
                               !heterogen_vl&heterogen_contacts~"Equal viral load, variable contacts",
                               !heterogen_vl&!heterogen_contacts~"Equal viral load and equal contacts"),
    heterogen_label=fct_relevel(heterogen_label,
                                "Variable viral load and contacts",
                                "Variable viral load, equal contacts",
                                "Equal viral load, variable contacts",
                                "Equal viral load and equal contacts")) %>%
  ggplot(aes(y=value,x=period,colour=name))+
  geom_point(aes(colour=heterogen_label,group=heterogen_label))+
  geom_line(aes(colour=heterogen_label,group=heterogen_label))+
  scale_linetype(name="")+
  scale_shape(name="")+
  ggnewscale::new_scale(new_aes="linetype")+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = quad_col_pal,name="")+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  facet_grid2(name~.,
              switch="y",
              independent = "y",
              scales="free_y",
              axes="all",
              remove_labels = "x",
              labeller=labeller(name=c("mu"="Mean R","size"="k of R",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)"),
                                heterogen_vl=c("TRUE"="Variable viral load trajectory",
                                               "FALSE"="Same viral load trajectory"),
                                heterogen_contacts=c("TRUE"="Variable contacts",
                                                     "FALSE"="Same contacts"))
  )+
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits = c(0,NA)),scale_y_log10()))+
  labs(y="",#"Mean parameter value",
       x="Time period"
  )+
  plotting_theme+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.direction = "vertical",
        legend.position = "right")


ggsave("results/R and k heterogeneity.png",width=200,height=150,dpi=600,units="mm",bg="white")
ggsave("results/R and k heterogeneity.pdf",width=200,height=150,dpi=600,units="mm",bg="white")


processed_infections_baseline %>%
  #filter.(prop_self_iso_test==0,sampling_freq==3) %>%
  ggplot(aes(x=vl,y=total_contacts,colour=total_infections))+geom_jitter(alpha=0.5)+
  scale_y_log10("Daily contacts")+
  labs(x="Viral load (RNA copies/ml)")+
  scale_colour_viridis_c("Infected contacts",trans="log10",na.value=NA,option ="turbo")+
  facet_grid(prop_self_iso_test~period+sampling_freq)+
  plotting_theme

ggsave("contacts_infections.png")

processed_infections_baseline %>%
  #filter.(prop_self_iso_test==0,sampling_freq==3) %>%
  ggplot(aes(x=vl,y=total_contacts,colour=total_infections))+
  geom_jitter(alpha=0.5)+
  scale_y_log10()+
  scale_colour_viridis_c(trans="log10",na.value=NA,option="turbo")+
  facet_grid(prop_self_iso_test~period+sampling_freq)+
  plotting_theme


testing_plot <-processed_infections_testing %>% 
  summarise.(sum_inf=sum(total_infections),.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size)) %>%
  summarise.(.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size,-sim),
             dist_means=list(fitdist(sum_inf,"nbinom")$estimate %>% enframe() %>% pivot_wider(names_from=name,values_from = value)),
             n=n(),
             ss_10=sum(sum_inf>10),
             ss_0=sum(sum_inf<=0)) %>% 
  mutate.(
    prop_ss_10=ss_10/n*100,
    prop_ss_0=ss_0/n*100) %>% 
  unnest.(dist_means)  %>% 
  filter.(variant=="wild") %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  ggplot(aes(y=value,x=prop_self_iso_test*100,colour=factor(sampling_freq),group=sampling_freq))+
  #geom_point()+
  geom_line()+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = tri_col_pal)+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  facet_grid2(name~period,
              #scales="free",
              scales="free_y",
              remove_labels = "x",
              axes = "all",
              #independent = "y",
              labeller=labeller(name=c("mu"="Mean R","size"="k of R",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)")),
              switch="y"
  )+
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits = c(0,3)),
                                    scale_y_log10(),
                                    scale_y_continuous(limits = c(0,NA)),
                                    scale_y_continuous(limits = c(0,NA))))+
  lims(y=c(0,NA))+
  labs(y="",
       x="Uptake of/adherence to lateral flow testing (%)",
       colour="Testing frequency (days between tests)")+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("results/lft_impact_testing.png",width=210,height=150,dpi=600,units="mm",bg="white")

## events

events_plot <- processed_infections_events %>% 
  summarise.(sum_inf=sum(total_infections),.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size)) %>%
  summarise.(.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size,-sim),
             dist_means=list(fitdist(sum_inf,"nbinom")$estimate %>% 
                               enframe()%>% 
                               pivot_wider(names_from=name,values_from = value)),
             n=n(),
             ss_10=sum(sum_inf>10),
             ss_0=sum(sum_inf<=0)) %>% 
  mutate.(
    prop_ss_10=ss_10/n*100,
    prop_ss_0=ss_0/n*100) %>% 
  unnest.(dist_means) %>% 
  drop_na.(event_size) %>% 
  filter.(variant=="wild") %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  ggplot(aes(y=value,x=prop_self_iso_test*100,colour=factor(event_size),group=event_size))+
  #geom_point()+
  geom_line()+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = tri_col_pal)+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  facet_grid2(name~period,
              #scales="free",
              scales="free_y",
              remove_labels = "x",
              #independent = "y",
              axes = "all",
              labeller=labeller(name=c("mu"="Mean R","size"="k of R",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)")),
              switch="y"
  )+
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits = c(0,3)),
                                    scale_y_log10(),
                                    scale_y_continuous(limits = c(0,NA)),
                                    scale_y_continuous(limits = c(0,NA))))+
  #lims(y=c(0,NA))+
  labs(y="Mean parameter value",
       x="Uptake of/adherence to pre-event lateral flow testing (%)",
       colour="Minimum event size\nto prompt testing")+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("results/lft_impact_events.png",width=210,height=150,dpi=600,units="mm",bg="white")

testing_plot/events_plot+plot_annotation(tag_levels = "A")
ggsave("results/lft_plot.png",dpi=600,width=210,height=325,units="mm",bg="white")
ggsave("results/lft_plot.pdf",dpi=600,width=210,height=300,units="mm",bg="white",device = cairo_pdf)

# boot_est %>%  
#   mutate.(boot_dist=map.(.x=dists, ~bootdist(f =.,bootmethod = "nonparam",parallel="snow",ncpus=8)$CI %>% 
#                           as.data.frame() %>% rownames_to_column())) %>% 
#   unnest.(boot_dist) %>% 
#   ggplot()+
#   geom_pointrange(aes(y=Median,ymin=`2.5%`,ymax=`97.5%`,x=factor(period)))+
#   facet_wrap(~rowname,scales="free_y",ncol = 1)+lims(y=c(0,NA))
# 
# ggsave("mu_size.png")

#prop infecting 0, >10, >20...
ss_dat <- processed_infections_testing %>% 
  filter.(period%in%c("Pre-pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  summarise.(sum_inf=sum(total_infections),.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size)) %>%
  summarise.(n=n(),
             ss_10=sum(sum_inf>10),
             ss_20=sum(sum_inf>20),
             ss_0=sum(sum_inf<=0),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size,-sim)) %>% 
  mutate.(
    prop_ss_10=ss_10/n,
    prop_ss_20 =ss_20/n,
    prop_ss_0=ss_0/n)

ss_dat %>% drop_na.(sampling_freq) %>% 
  pivot_longer.(c(prop_ss_0,prop_ss_10)) %>% 
  ggplot(aes(x=prop_self_iso_test,y=value,group=sampling_freq,colour=factor(sampling_freq)))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values = tri_col_pal)+
  scale_x_continuous(labels=scales::percent,breaks=breaks_width(0.5))+
  scale_y_continuous(labels=scales::percent)+
  labs(y="Mean parameter value",
       x="Uptake of/adherence to lateral flow testing",
       #title="The relative impact of lateral flow testing is greatest when individuals have lots of contacts,\nwhen uptake is high, and when testing is frequent",
       colour="Testing frequency (days between tests)")+
  facet_rep_grid(name~period,scales="free_y",switch="y",labeller = labeller(name=c("prop_ss_0"="0 sec. inf.",
                                                                                   "prop_ss_10"=">10 sec. inf.")))+
  plotting_theme

ggsave("results/prop_ss_testing.png",width=210,height=150,dpi=600,units="mm",bg="white")

#prop infecting 0, >10, >20...
ss_dat <- processed_infections_events %>% 
  filter.(period%in%c("Pre-pandemic","Lockdown 1","Relaxed restrictions","School reopening")) %>% 
  summarise.(sum_inf=sum(total_infections),.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size)) %>%
  summarise.(n=n(),
             ss_10=sum(sum_inf>10),
             ss_20=sum(sum_inf>20),
             ss_0=sum(sum_inf<=0),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,event_size,-sim)) %>% 
  mutate.(
    prop_ss_10=ss_10/n,
    prop_ss_20 =ss_20/n,
    prop_ss_0=ss_0/n)

ss_dat %>% drop_na.(event_size) %>% 
  pivot_longer.(c(prop_ss_0,prop_ss_10)) %>% 
  ggplot(aes(x=prop_self_iso_test,y=value,group=event_size,colour=factor(event_size)))+
  geom_point()+
  geom_line()+
  scale_colour_manual(values = tri_col_pal)+
  scale_x_continuous(labels=scales::percent,breaks=breaks_width(0.5))+
  scale_y_continuous(labels=scales::percent)+
  labs(y="Mean parameter value",
       x="Uptake of/adherence to pre-event lateral flow testing",
       #title="The relative impact of lateral flow testing is greatest when individuals have lots of contacts,\nwhen uptake is high, and when testing is frequent",
       colour="Minimum number of attendees\nrequired to prompt pre-event testing")+
  facet_rep_grid(name~period,scales="free_y",switch="y",labeller = labeller(name=c("prop_ss_0"="0 sec. inf.",
                                                                                   "prop_ss_10"=">10 sec. inf.")))+
  plotting_theme

ggsave("results/prop_ss_events.png",width=210,height=150,dpi=600,units="mm",bg="white")

ss_10_plot <- ss_dat %>% 
  filter.(name=="prop_ss_10") %>% 
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
  filter.(name=="prop_ss_20") %>% 
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


(ss_10_plot/ss_20_plot)&scale_color_manual(name = "",guide=F,
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
processed_infections_baseline %>% 
  summarise.(sum_inf=sum(total_infections),
             sum_contacts=sum(total_contacts),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>%
  mutate.(sampling_freq=case_when.(sampling_freq==3~"Testing every 3 days",
                                   TRUE~ "No testing"),
          #hh_SAR=n_hh_infected/contacts_hh,
          #nhh_SAR=n_nhh_infected/contacts_nhh,
          tot_SAR=sum_inf/sum_contacts) %>% 
  pivot_longer.(cols = c("tot_SAR"),values_to = "value",names_to =  "param") %>% 
  summarise.(x = quantile(value,probs=c(0.5,0.025,0.975),na.rm = T),q=c(0.5,0.025,0.975),.by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim,param)) %>% 
  pivot_wider.(names_from = q,values_from = x) %>% 
  ggplot(aes(x=period,y=`0.5`,ymin=`0.025`,ymax=`0.975`))+
  geom_pointrange()+
  facet_grid(param~sampling_freq+variant,scales="free")

#Heatmap
heatmap <- processed_infections_baseline %>% 
  #filter.(prop_self_iso_test==0,is.na(sampling_freq)) %>%   
  #filter.(n_nhh_infected>0,n_hh_infected>0) %>%  
  ggplot(aes(x=ct,
             y=hh_contacts+nhh_contacts,
             colour=hh_infected+nhh_infected))+
  geom_point()+
  scale_y_log10()+
  scale_colour_viridis_c(trans="log",option = "turbo",breaks=c(1,10,100),limits=c(1,250), guide = guide_colourbar(nbin = 100, draw.ulim = TRUE, draw.llim = TRUE))+
  labs(x="Viral load (CT)",y="Exposed contacts",colour="Infected contacts")+
  theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        legend.position = "right",
        strip.placement = "outside")+
  facet_grid(~period,labeller=labeller(time_period=c("pre"="Pre-pandemic (BBC 2018)",
                                                     "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                                                     "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                                       name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                                       prop_self_iso=function(x){scales::percent(as.numeric(x))}))

ggsave("results/heatmap.png",width=210,height=150,dpi=600,units="mm")

hex_map <- res %>% 
  filter( prop_self_iso_symp==1,prop_self_iso_test==1,is.na(sampling_freq)) %>%   
  filter(n_nhh_infected>0,n_hh_infected>0) %>%  
  ggplot(aes(x=ct,
             y=contacts_hh+contacts_nhh))+
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
  filter(n_nhh_infected>0,n_hh_infected>0) %>%  
  pivot_longer(cols=c(n_hh_infected,n_nhh_infected)) %>% 
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
  filter(n_nhh_infected>0,n_hh_infected>0) %>%  
  pivot_longer(cols=c(n_hh_infected,n_nhh_infected)) %>% 
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
