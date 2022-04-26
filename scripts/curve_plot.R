
plot_dat <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE,FALSE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = 1000,asymp_parms=asymp_fraction,variant_info=.x,browsing = F)) %>% 
  bind_rows.()%>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end,interval=0.1
  )))  %>%
  unnest.(infectiousness) %>%
  crossing.(lower_inf_thresh = c(FALSE)) %>%
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
  mutate.(vl=10^(vl)) %>% 
  ggplot()+
  geom_line(data=. %>% filter.(heterogen_vl==T),
            aes(x=t,
                y=vl,
                group=sim,
                colour=culture_p),
            alpha=0.05
  )+
  geom_line(data=. %>% filter.(heterogen_vl==F),
            aes(x=t,
                y=vl,
                group=sim,
                colour=culture_p),
            alpha=0.1,size=1
  )+
  scale_colour_viridis_c(name="Probability of being infectious",option="inferno",begin = 0.2,end=0.8,
                         guide=guide_colorsteps(barwidth=unit(5,"cm")))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  scale_y_log10(name="RNA copies/ml",labels=label_log())+
  coord_cartesian(ylim=c(10^3.5,NA))+
  plotting_theme

#ggsave("results/log_plot.png",dpi=600,width=210,height=150,units="mm",bg="white")

cv_plot <- plot_dat %>% 
  filter.(heterogen_vl==T,t%%1==0,lower_inf_thresh==F) %>% 
  summarise.(cv=cv(culture_p),.by=t) %>% 
  ggplot(aes(x=t,y=cv))+
  geom_point(colour="#2E4C6D")+
  geom_line(colour="#2E4C6D")+
  lims(y=c(0,NA))+
  labs(x="Days since first detectable by PCR (Ct<40)",y="Coefficient of variation\nin infectiousness")+
  plotting_theme


auc_plot_lab <- plot_dat %>% 
  filter.(heterogen_vl==T,t%%1==0,lower_inf_thresh==F) %>% 
  summarise.(sum_inf=sum(culture_p),.by=sim) %>% 
  summarise.(cv=cv(sum_inf)) %>% pull(cv)


auc_plot <- plot_dat %>% 
  filter.(heterogen_vl==T) %>% 
  summarise.(sum_inf=sum(culture_p),.by=sim) %>% 
  ggplot(aes(x=sum_inf))+
  geom_density(fill="#FC997C",colour="#2E4C6D",alpha=0.1)+
  geom_text(aes(x=Inf,y=Inf),label=paste("CV = ",round(auc_plot_lab,digits = 2)),size=6,hjust=1.1,vjust=1.5,colour="#2E4C6D")+
  lims(y=c(0,NA))+
  labs(x="Area under infectiousness curve (AU)",y="Density")+
  plotting_theme

log_plot/(auc_plot+cv_plot)+plot_annotation(tag_levels = "A")
ggsave("results/log_and_cv_plot.png",dpi=600,width=210,height=150,units="mm",bg="white")


