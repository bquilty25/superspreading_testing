#### Results ####
source("scripts/utils.R")

#load simulation output
processed_infections_baseline <- qread("results/processed_infections_baseline.qs")
processed_infections_heterogen_on_off <- qread("results/processed_infections_heterogen_on_off.qs")
processed_infections_testing <- qread("results/processed_infections_testing.qs")
processed_infections_events <- qread("results/processed_infections_events.qs")
processed_infections_sens <- qread("results/processed_infections_sens.qs")

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
tic()
boot_res <- processed_infections_baseline %>% 
  #filter(period=="Pre-pandemic") %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>%
  summarise.(dists=list(bootdist(fitdist(sum_inf,"nbinom"),
                                 bootmethod="nonparam",
                                 parallel="multicore",
                                 ncpus=8)$CI%>% 
                          as.data.frame() %>% 
                          rownames_to_column(var = "name") %>% 
                          rename("lo"=`2.5%`,
                                 "hi"=`97.5%`)),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim))
toc()
qsave(boot_res,"results/R_and_k_bootstrap_ests.qs")

boot_res_sum <- boot_res %>% 
  unnest.(dists) %>%  
  filter.(variant=="wild") %>% 
  mutate.(name=as.factor(name),
          name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  as_tibble() 
write.csv(boot_res_sum,"results/R_and_k_bootstrap_ests.csv")

boot_res_sum %>% 
  ggplot(aes(y=Median,ymin=lo,ymax=hi,x=period,colour=name,fill=name,group=name))+
  geom_line()+
  geom_point()+
  geom_ribbon(alpha=0.4,colour=NA)+
  geom_segment(data=rt_by_time_period %>% mutate(name="mu"),
               aes(x=period,
                   xend=period,
                   y=lo,
                   yend=hi),
               alpha=0.25,
               size=10)+
  facet_grid2(rows=vars(name),switch="y",scales="free",
              labeller=labeller(name=c("mu"="R","size"="k",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)")),
              axes="all",
              remove_labels = "x",
  )+
  facetted_pos_scales(y=list(scale_y_continuous(limits=c(0,NA)),
                             #scale_y_continuous(),
                             #scale_y_continuous(),
                             scale_y_log10(limits=c(NA,NA),breaks=log_breaks())))+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = bi_col_pal,guide="none")+
  scale_fill_manual(values = bi_col_pal,guide="none")+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  lims(y=c(0,NA))+
  labs(y="",
       x="Time period")+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("results/Fig4 - R and k over time.png",width=150,height=100,dpi=600,units="mm",bg="white")
ggsave("results/Fig4 - R and k over time.pdf",width=150,height=100,units="mm",bg="white")

# heterogen_on_off
tic()
boot_res_heterogen <- processed_infections_heterogen_on_off %>% 
  #filter(period=="Pre-pandemic") %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test)) %>%
  summarise.(dists=list(bootdist(fitdist(sum_inf,"nbinom"),
                                 bootmethod="nonparam",
                                 parallel="multicore",
                                 ncpus=8)$CI%>% 
                          as.data.frame() %>% 
                          rownames_to_column(var = "name") %>% 
                          rename("lo"=`2.5%`,
                                 "hi"=`97.5%`)),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,-sim))
toc()
qsave(boot_res_heterogen,"results/R_and_k_bootstrap_ests_heterogen.qs")

# other_est <- tribble(~study,~xmin,~xmax,~ymin,~ymax, ~y,
#                      "Endo et al. 2020", -Inf, Inf, 0.05, 0.2, 0.1,
#                      "Rio & Althaus 2020", -Inf, Inf, 0.014, 6.95, 0.54,
#                      "Adam et al. 2020", -Inf, Inf, 0.45, 0.72, 0.58,
#                      "Laxminarayan et al. 2020", -Inf, Inf, 0.49, 0.52, 0.51)

boot_res_heterogen_sum <- boot_res_heterogen %>% 
  unnest.(dists) %>% 
  #pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0")) %>% 
  filter.(variant=="wild",!(heterogen_vl==F&heterogen_contacts==F&name=="size"),
          name%in%c("size")) %>% 
  mutate.(heterogen_label=case_when.(heterogen_vl&heterogen_contacts~"Variable viral load, overdispersed contacts",
                               heterogen_vl&!heterogen_contacts~"Variable viral load, Poisson contacts",
                               !heterogen_vl&heterogen_contacts~"Equal viral load, overdispersed contacts",
                               !heterogen_vl&!heterogen_contacts~"Equal viral load and Poisson contacts"),
    heterogen_label=fct_relevel(heterogen_label,
                                "Variable viral load, overdispersed contacts",
                                "Variable viral load, Poisson contacts",
                                "Equal viral load, overdispersed contacts",
                                "Equal viral load and Poisson contacts"),
    heterogen_vl=ifelse(heterogen_vl,"Heterogeneous viral load","Homogeneous viral load"),
    heterogen_contacts=ifelse(heterogen_contacts,"Heterogeneous contacts","Homogeneous contacts"))
write.csv(boot_res_heterogen_sum,"results/R_and_k_bootstrap_ests_heterogen.csv")

(heterogen_plot <- (boot_res_heterogen_sum %>% ggplot(aes(y=Median,ymin=lo,ymax=hi,x=period,colour=name))+
    geom_line(aes(colour=heterogen_label,fill=heterogen_label,group=heterogen_label,linetype=heterogen_label))+
    geom_point(aes(colour=heterogen_label,fill=heterogen_label,group=heterogen_label,linetype=heterogen_label))+
  geom_lineribbon(aes(colour=heterogen_label,fill=heterogen_label,group=heterogen_label,linetype=heterogen_label),
                  alpha=0.4)+
  scale_colour_manual(values = tri_col_pal)+
    scale_fill_manual(values = tri_col_pal)+  
  #scale_linetype_manual(values=c("solid","dashed","dashed"))+
  labs(y="",#"Mean parameter value",
         x="Time period",
       linetype = "",
       colour="",
       fill = ""
    )+
  facet_grid2(name~.,
              switch="y",
              independent = "y",
              scales="free_y",
              axes="all",
              remove_labels = "x",
              labeller=labeller(name=c("mu"="Mean R","size"="Overdispersion (k)",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)"),
                                heterogen_vl=c("TRUE"="Variable viral load trajectory",
                                               "FALSE"="Same viral load trajectory"),
                                heterogen_contacts=c("TRUE"="Variable contacts",
                                                     "FALSE"="Same contacts"))
  )+
  scale_y_log10()#|
    # ggplot(other_est,aes(x=study, ymin=ymin, ymax=ymax, y = y, colour=study))+
    # geom_pointrange(fatten=4, alpha=0.5)+
    # scale_y_log10(limit=c(0.01,10))
  + coord_cartesian(ylim = c(0.1,10))
  )&
    plotting_theme&
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.direction = "vertical",
          legend.position = "bottom")
  )


ggsave(heterogen_plot,file = "results/Fig5 - R and k heterogeneity.png",width=200,height=150,dpi=600,units="mm",bg="white")
ggsave(heterogen_plot,file = "results/Fig5 - R and k heterogeneity.pdf",width=200,height=150,units="mm",bg="white")

#### Sensitivity analysis ----

processed_infections_baseline %>% 
  filter.(period=="Pre-pandemic") %>% 
  mutate.(contacts="Unadjusted") %>% 
  bind_rows.(processed_infections_sens %>% 
               mutate.(contacts="Adjusted")) %>% 
  summarise.(sum_inf=sum(total_infections),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,contacts)) %>%
  summarise.(dists=list(fitdist(sum_inf,"nbinom")),
             dist_means=list(fitdist(sum_inf,"nbinom")$estimate %>% enframe() %>% pivot_wider(names_from=name,values_from = value)),
             n=n(),
             ss_10=sum(sum_inf>10),
             ss_20=sum(sum_inf>20),
             ss_0=sum(sum_inf<=0),
             .by=c(all_of(key_grouping_var),sampling_freq,prop_self_iso_test,contacts,-sim)) %>% 
  mutate.(
    prop_ss_10=ss_10/n*100,
    #prop_ss_20 =ss_20/n,
    prop_ss_0=ss_0/n*100) %>% 
  unnest.(dist_means) %>% 
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>% 
  mutate.(name=fct_relevel(name,"mu","size","prop_ss_10","prop_ss_0"),
          contacts=fct_relevel(contacts,"Unadjusted")) %>% 
  filter.(variant=="wild") %>% 
  ggplot(aes(y=value,x=contacts,colour=name,group=name,shape=contacts))+
  geom_point(show.legend = F)+
  # geom_text_repel(data=. %>% filter.(period=="Pre-pandemic",name=="mu"),
  #                 aes(x=contacts,y=value,label=paste0("R0 = ",round(value,1))),family="Lato",
  #                 nudge_y = -0.2)+
  geom_hline(aes(linetype=name,yintercept=1),colour=quad_col_pal[1])+
  scale_colour_manual(values = quad_col_pal,guide="none")+
  scale_linetype_manual(values=c("dashed",NA,NA,NA),guide="none")+
  facet_grid2(~name,switch="y",scales="free_y",independent = "y",
              labeller=labeller(name=c("mu"="R","size"="k",
                                       "prop_ss_0"= "Proportion infecting\n 0 others (%)",
                                       "prop_ss_10"="Proportion infecting\n over 10 others (%)")),
              axes="all",
              remove_labels = "x")+
  lims(y=c(0,NA))+
  labs(y="",
       x="Contact data tail adjustment")+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))
  
ggsave("results/R_and_k_contact_data_sensitivity.png",width=200,height=100,dpi=600,units="mm",bg="white")

processed_infections_baseline %>%
  #filter.(prop_self_iso_test==0,sampling_freq==3) %>%
  ggplot(aes(x=vl,y=total_contacts,colour=total_infections))+geom_jitter(alpha=0.5)+
  scale_y_log10("Daily contacts")+
  labs(x="Viral load (RNA copies/ml)")+
  scale_colour_viridis_c("Infected contacts",trans="log10",na.value=NA,option ="turbo")+
  facet_grid(prop_self_iso_test~period+sampling_freq)+
  plotting_theme

# ggsave("contacts_infections.png")

#### LFT testing ----

## regular testing ----

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
              labeller=labeller(name=c("mu"="R","size"="k",
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
  theme(axis.text.x=element_text(angle = 0, vjust = 1, hjust=1))

ggsave("results/lft_impact_testing.png",width=210,height=150,dpi=600,units="mm",bg="white")

## events ----

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
              labeller=labeller(name=c("mu"="R","size"="k",
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
  theme(axis.text.x=element_text(angle = 0, vjust = 1, hjust=1))

ggsave("results/lft_impact_events.png",width=210,height=150,dpi=600,units="mm",bg="white")
ggsave("results/lft_impact_events.pdf",width=210,height=150,units="mm",bg="white")


testing_plot/events_plot+plot_annotation(tag_levels = "A")
ggsave("results/Fig6 - lft_plot.png",dpi=600,width=210,height=325,units="mm",bg="white")
ggsave("results/Fig6 - lft_plot.pdf",width=210,height=300,units="mm",bg="white")