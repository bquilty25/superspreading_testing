
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
                #colour=tri_col_pal[2]),
                #colour=zoo::rollmean(culture_p,4,fill=0)
                #,
            alpha=0.05
  )+
  geom_line(data=. %>% filter.(heterogen_vl==F),
            aes(x=t,
                y=vl,
                group=sim,
                colour=culture_p),
                #colour=tri_col_pal[2],
                #colour=zoo::rollmean(culture_p,4,fill=0)
                
            #alpha=0.1,
            size=1
  )+
  scale_colour_gradient(high=bi_col_pal[2],low=bi_col_pal[1])+
  # scale_colour_viridis_c(name="Probability of infectiousness",option="inferno",begin = 0.2,end=0.8,limits=c(0,1),
  #                        guide=guide_colorsteps(barwidth=unit(3,"cm"),barheight=unit(0.5,"cm"),show.limits = T))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  scale_y_log10(name="RNA copies/ml",labels=label_log())+
  coord_cartesian(xlim=c(NA,20),ylim=c(10^3.5,NA))+
  plotting_theme+guides(colour="none")

ggsave("results/log_plot.png",dpi=600,width=210,height=150,units="mm",bg="white")

# cv_plot <- plot_dat %>% 
#   filter.(heterogen_vl==T,t%%1==0,lower_inf_thresh==F) %>% 
#   summarise.(cv=cv(culture_p),.by=t) %>% 
#   ggplot(aes(x=t,y=cv))+
#   geom_point(colour="#2E4C6D")+
#   geom_line(colour="#2E4C6D")+
#   lims(y=c(0,NA))+
#   labs(x="Days since first detectable by PCR (Ct<40)",y="Coefficient of variation\nin infectiousness")+
#   plotting_theme

# 
# auc_plot_lab <- plot_dat %>% 
#   filter.(heterogen_vl==T,t%%1==0,lower_inf_thresh==F) %>% 
#   summarise.(sum_inf=sum(culture_p),.by=sim) %>% 
#   summarise.(cv=cv(sum_inf)) %>% pull(cv)


auc_plot <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = 1000,asymp_parms=asymp_fraction,variant_info=.x,browsing = F)) %>% 
  bind_rows.()%>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end,interval=1
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
  select.(-c(prolif, start, end)) %>% 
  summarise.(sum_inf=sum(culture_p),.by=sim) %>% 
  ggplot(aes(x=sum_inf))+
  geom_density(fill=bi_col_pal[2],colour=bi_col_pal[2],alpha=0.25,adjust=2)+
  #geom_text(aes(x=Inf,y=Inf),label=paste("CV = ",round(auc_plot_lab,digits = 2)),size=6,hjust=1.1,vjust=1.5,colour="#2E4C6D")+
  lims(y=c(0,NA))+
  labs(x="Area under infectiousness curve (AU)",y="Density")+
  plotting_theme

#log_plot/(auc_plot+cv_plot)+plot_annotation(tag_levels = "A")
#ggsave("results/log_and_cv_plot.png",dpi=600,width=210,height=150,units="mm",bg="white")


# jitter_plot <- plot_dat %>% 
#   mutate.(vl=10^(vl)) %>% 
#   ggplot()+
#   geom_jitter(data=. %>% filter.(heterogen_vl==T),
#             aes(x=t,
#                 y=culture_p,
#                 group=sim,
#                 colour=culture_p),
#             alpha=0.05
#   )+
#   scale_colour_viridis_c(name="Probability of infectiousness",option="inferno",begin = 0.2,end=0.8,
#                          guide=guide_colorsteps(barwidth=unit(5,"cm"),barheight=unit(1,"cm")))+
#   scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
#  # scale_y_log10(name="RNA copies/ml",labels=label_log())+
#   #coord_cartesian(ylim=c(10^3.5,NA))+
#   plotting_theme

# dens_plot <- plot_dat %>% 
#   mutate.(vl=10^(vl)) %>% 
#   ggplot()+
#   geom_density(data=. %>% filter.(heterogen_vl==T),aes(x=culture_p)
#   )+
#   coord_flip()+
#   scale_colour_viridis_c(name="Probability of infectiousness",option="inferno",begin = 0.2,end=0.8,
#                          guide=guide_colorsteps(barwidth=unit(5,"cm")))+
#   scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
#   # scale_y_log10(name="RNA copies/ml",labels=label_log())+
#   #coord_cartesian(ylim=c(10^3.5,NA))+
#   plotting_theme

#dens_plot+jitter_plot

violin_plot <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = 1000,asymp_parms=asymp_fraction,variant_info=.x,browsing = F)) %>% 
  bind_rows.()%>% 
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end,interval=1
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
  select.(-c(prolif, start, end))  %>% 
  mutate.(median_p=mean(culture_p),.by=c(t,heterogen_vl)) %>% 
  filter.(t<20) %>% 
  ggplot()+
  geom_violin(data=. %>% filter.(heterogen_vl==T),
                          aes(x=t,
                              y=culture_p,
                              group=t,
                             fill=median_p
                              ),colour="white",scale="width",draw_quantiles = c(0.25, 0.5, 0.75))+
  # ggdist::stat_slab(data=. %>% filter.(heterogen_vl==T),
  #             aes(x=factor(t),
  #                 y=culture_p,
  #                 group=t,
  #                fill=median_p
  #                 ), #fill=bi_col_pal[2],
  #             side="both",
  #             normalize="groups",
  #             breaks=seq(0,1,by=0.05)
  # )+
  scale_fill_gradient(high=bi_col_pal[2],low=bi_col_pal[1],guide="none")+
  # scale_fill_viridis_c(guide="none",option="inferno",begin = 0.2,end=0.8#,limits=c(0,1)
  #                        #guide=guide_colorsteps(barwidth=unit(5,"cm"))
  #                      )+
  labs(x="Days since first detectable by PCR (Ct<40)",y="Probability of infectiousness")+
  coord_cartesian(xlim=c(NA,20))+
  #scale_x_discrete(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  #scale_y_log10(name="RNA copies/ml",labels=label_log())+
  #coord_cartesian(ylim=c(10^3.5,NA))+
  plotting_theme


inf_plot <- pickering %>% 
  mutate(vl=10^vl) %>% 
  ggplot()+
  geom_smooth(aes(x=vl,y=culture),
              method="glm",
              method.args=list(family="binomial"),
              colour=tri_col_pal[1],
              fill=tri_col_pal[1],
              alpha=0.25)+
  scale_x_log10(name="RNA copies/ml",labels=label_log())+
  ylab("Probability of infectiousness")+
  plotting_theme

((log_plot|inf_plot)/(violin_plot|auc_plot))+plot_annotation(tag_levels = "A")

ggsave("results/log_and_violin_plot.png",dpi=600,width=210,height=120,units="mm",bg="white")

