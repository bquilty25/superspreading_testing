
plot_dat <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen=c(TRUE,FALSE)) %>% 
  group_split.(variant,heterogen) %>% 
  map.(~make_trajectories(n_sims = 100,asymp_parms=asymp_fraction,variant_info=.x,browsing = F)) %>% 
  bind_rows.()%>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end,interval=0.1
  )))  %>%
  unnest.(infectiousness) %>%
  crossing.(lower_inf_thresh = c(TRUE, FALSE)) %>%
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

log_plot <- plot_dat %>% 
  #mutate.(vl=10^(vl)) %>% 
  ggplot()+
  geom_line(aes(x=t,
                y=vl,
                group=sim,
                colour=test_p),
            alpha=0.2
  )+
  scale_colour_gradient2(name="Probability of detection by LFT",midpoint=0.5,low=muted("blue"),high=muted("red"))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  scale_y_continuous(name="log10 RNA copies/ml")+facet_wrap(~heterogen)
  # scale_y_reverse(# Features of the first axis
  #   limits=c(40,NA),
  #   name = "Cycle\nthreshold")+
  #coord_cartesian(ylim=c(40,NA))+
  #plotting_theme+
  # theme( strip.background = element_blank(),
  #        panel.grid.minor.x = element_blank(),
  #        panel.grid.major.x = element_blank(),
  #        #strip.text.x = element_blank(),
  #        #axis.title.x = element_blank()
  # )
#+coord_cartesian(ylim=c(NA,1000000))
#facet_wrap(~type,nrow = 2,labeller=labeller(type=capitalize))

ggsave("results/log_plot.png",dpi=600,width=210,height=150,units="mm")


cont_plot <- plot_dat %>% 
  #mutate.(vl=10^(vl)) %>% 
  ggplot()+
  geom_line(aes(x=t,
                y=ct,
                group=sim,
                colour=test_p),
            alpha=0.2
  )+
  # geom_point(data=. %>% filter(type=="symptomatic"),
  #            aes(x=onset_t,y=onset_ct,group=interaction(sim,idx,variant),colour=onset_p),pch=16,alpha=0.5)+
  scale_colour_gradient2(name="Probability of detection by LFT",midpoint=0.5,low=muted("blue"),high=muted("red"))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  scale_y_reverse(name="RNA copies/ml",sec.axis = sec_axis(trans=convert_Ct_logGEML))+
  # scale_y_reverse(# Features of the first axis
  #   limits=c(40,NA),
  #   name = "Cycle\nthreshold")+
  #coord_cartesian(ylim=c(40,NA))+
  # plotting_theme+
  # theme( strip.background = element_blank(),
  #        panel.grid.minor.x = element_blank(),
  #        panel.grid.major.x = element_blank()
  #        #strip.text.x = element_blank(),
  #        #axis.title.x = element_blank()
  # )+ 
  facet_zoom(ylim=c(c(0,1e6)))
#+coord_cartesian(ylim=c(NA,1000000))
  #facet_wrap(~type,nrow = 2,labeller=labeller(type=capitalize))


ggsave("results/cont_plot.png",dpi=600,width=210,height=150,units="mm")
