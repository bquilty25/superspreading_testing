source("scripts/utils.R")

n_sims <- 5000

traj <- vl_params %>% 
  filter.(variant%in%c("wild")) %>%
  mutate.(variant=fct_drop(variant)) %>% 
  crossing(heterogen_vl=c(TRUE,FALSE)) %>% 
  group_split.(variant,heterogen_vl) %>% 
  map.(~make_trajectories(n_sims = n_sims,asymp_parms=asymp_fraction,variant_info=.x,browsing = F)) %>% 
  bind_rows.()

infctsnss_params <- generate_params(culture_mod, n_sims) %>%
  as_tidytable(.) %>%
  rename.(beta0 = "(Intercept)", beta1 = vl) %>%
  mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

plot_dat <- traj %>% 
  mutate.(infectivity = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end,interval=1
  )))  %>%
  unnest.(infectivity) %>%
  crossing.(lower_inf_thresh = c(FALSE)) %>%
  mutate.(
    culture_p = culture_prob(vl, beta0, beta1),
    infectious = rbernoulli(n = n(),
                            p = culture_p),
    test_p = stats::predict(
      object =  innova_mod,
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    test = rbernoulli(n = n(),
                      p = test_p),
    .by = c(lower_inf_thresh)
  ) %>%
  replace_na.(list(test       = FALSE,
                   infectious = FALSE)) %>% 
  select.(-c(prolif, start, end)) 

plot_dat1 <- plot_dat %>% 
  filter.(sim <= 1000)

log_plot <- plot_dat1 %>% 
  filter.(heterogen_vl==T) %>%
  mutate.(vl=10^(vl)) %>% 
  ggplot()+
  geom_line(data=. %>% filter.(heterogen_vl==T),
            aes(x=t,
                y=vl,
                group=sim,
                #colour=culture_p
                ),
                colour=tri_col_pal[1],
                #colour=zoo::rollmean(culture_p,4,fill=0)
                #,
            alpha=0.05
  )+
  # geom_line(data=. %>% filter.(heterogen_vl==F) %>% filter.(sim==1),
  #           aes(x=t,
  #               y=vl,
  #               group=sim,
  #               #colour=culture_p
  #               ),
  #               colour=tri_col_pal[1],
  #               #colour=zoo::rollmean(culture_p,4,fill=0)
  #               
  #           #alpha=0.1,
  #           size=1
  # )+
  geom_hline(yintercept=1.6e7,linetype="dashed",hjust=1)+
  #scale_colour_gradient(high=bi_col_pal[2],low=bi_col_pal[1])+
  # scale_colour_viridis_c(name="Probability of infectivity",option="inferno",begin = 0.2,end=0.8,limits=c(0,1),
  #                        guide=guide_colorsteps(barwidth=unit(3,"cm"),barheight=unit(0.5,"cm"),show.limits = T))+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  scale_y_log10(name="RNA copies/ml",labels=label_log())+
  coord_cartesian(xlim=c(NA,20),ylim=c(10^3.5,NA))+
  plotting_theme+guides(colour="none")+
  labs(caption="Kissler et al. 2021")

ggsave("results/log_plot.png",dpi=600,width=210,height=150,units="mm",bg="white")

# Average number of days spent infectious
plot_dat %>% 
  summarise.(n_inf=sum(infectious==T),.by=c(sim,heterogen_vl)) %>% 
  summarise.(q=list(quibble2(n_inf,c(0.025,0.5,0.975))),.by=heterogen_vl) %>% 
  unnest.(q)

days_inf_plot <- plot_dat %>% 
  summarise.(n_inf=sum(infectious==T),.by=c(sim,heterogen_vl)) %>% 
  filter.(heterogen_vl==T) %>% 
  ggplot()+
  geom_bar(aes(x=n_inf,y=..count../sum(..count..)),fill=bi_col_pal[2])+
  ylab("Probability")+
  scale_x_continuous(name="Number of days individuals infected contacts",
                     breaks = breaks_width(1),
                     expand = c(0.04,0.04))+
  plotting_theme+guides(colour="none")


culture_plot <- plot_dat1 %>%
  filter.(heterogen_vl==T) %>% 
  ggplot()+
  geom_line(data=. %>% filter.(heterogen_vl==T),
            aes(x=t,
                y=culture_p,
                group=sim,
                #colour=culture_p
            ),
            colour=bi_col_pal[2],
            #colour=zoo::rollmean(culture_p,4,fill=0)
            #,
            alpha=0.05
  )+
  ylab("Probability of infectivity")+
  scale_colour_gradient(high=bi_col_pal[2],low=bi_col_pal[1],guide="none")+
  scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  #scale_y_log10(name="RNA copies/ml",labels=label_log())+
  #coord_cartesian(xlim=c(NA,20),ylim=c(10^3.5,NA))+
  plotting_theme+guides(colour="none")#+
  #labs(title="Inferred infectivity curves"
       #,subtitle = "Duration of P(infectivity)>0.5:\n1.9 days (95% CI: 0, 5.8)")


auc_dat <- plot_dat %>%
  filter.(heterogen_vl==T) %>%
  summarise.(sum_inf=sum(culture_p),.by=sim) 

q95 <- auc_dat[,quantile(sum_inf, c(0.025, 0.975))]
q95[2]/q95[1]

auc_gamma_params <- fitdist(auc_dat[,sum_inf], "gamma")$estimate

auc_plot <- auc_dat %>% 
  ggplot(aes(x=sum_inf))+
  geom_density(fill=bi_col_pal[1],colour=bi_col_pal[1],alpha=0.25,adjust=2)+
  #geom_text(aes(x=Inf,y=Inf),label=paste("CV = ",round(auc_plot_lab,digits = 2)),size=6,hjust=1.1,vjust=1.5,colour="#2E4C6D")+
  lims(y=c(0,NA))+
  labs(x="Area under infectivity curve (AU)",y="Probability density")+
  plotting_theme

(days_inf_plot/auc_plot)+plot_annotation(tag_levels =  "A")
ggsave("results/days_inf_and_auc_plot.png",dpi=600,width=210,height=100,units="mm",bg="white")
ggsave("results/days_inf_and_auc_plot.pdf",dpi=600,width=210,height=100,units="mm",bg="white")

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
#   scale_colour_viridis_c(name="Probability of infectivity",option="inferno",begin = 0.2,end=0.8,
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
#   scale_colour_viridis_c(name="Probability of infectivity",option="inferno",begin = 0.2,end=0.8,
#                          guide=guide_colorsteps(barwidth=unit(5,"cm")))+
#   scale_x_continuous(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
#   # scale_y_log10(name="RNA copies/ml",labels=label_log())+
#   #coord_cartesian(ylim=c(10^3.5,NA))+
#   plotting_theme

#dens_plot+jitter_plot

violin_plot <- plot_dat1  %>% 
  filter.(heterogen_vl==T) %>%
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
  labs(title="Daily variation in infectivity",x="Days since first detectable by PCR (Ct<40)",y="Probability of infectivity")+
  coord_cartesian(xlim=c(NA,20))+
  #scale_x_discrete(name="Days since first detectable by PCR (Ct<40)",breaks = breaks_width(5))+
  #scale_y_log10(name="RNA copies/ml",labels=label_log())+
  #coord_cartesian(ylim=c(10^3.5,NA))+
  plotting_theme

prob_culture <- infctsnss_params %>%
  crossing(ct=seq(10, 40, by=0.5)) %>%
  mutate.(vl=convert_Ct_logGEML(ct), 
          culture_p = culture_prob(vl, beta0, beta1)) %>%
  arrange.(sim)

inf_plot <- prob_culture %>%
  mutate.(vl=10^vl) %>%
  group_by(vl) %>%
  summarise(
    median_p = median(culture_p, na.rm = TRUE),
    lower = quantile(culture_p, 0.025, na.rm = TRUE),
    upper = quantile(culture_p, 0.975, na.rm = TRUE)
  ) %>%
  ggplot(aes(x = vl)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = bi_col_pal[1], alpha = 0.3) +
  geom_line(aes(y = median_p), colour = bi_col_pal[1], size = 1) +
  geom_vline(xintercept=1.6e7,linetype="dashed")+
  scale_x_log10(name="RNA copies/ml",labels=label_log(),limits=c(10^3.5,NA)) +
  scale_colour_gradient(high=bi_col_pal[2],low=bi_col_pal[1],guide = "none") +
  labs(y="Probability of\nculturing virus",
       caption="Pickering et al. 2021") +
  plotting_theme

# inf_plot <- pickering %>% 
#   mutate(vl=10^vl) %>% 
#   ggplot()+
#   stat_smooth(aes(x=vl,y=culture,
#                   colour=..y..),
#                   method="glm",
#                   method.args=list(family="binomial"),
#                   #alpha=0.25,
#               size=2,
#               geom="line")+
#   geom_vline(xintercept=1.6e7,linetype="dashed")+
#   scale_x_log10(name="RNA copies/ml",labels=label_log(),limits=c(10^3.5,NA))+
#   scale_colour_gradient(high=bi_col_pal[2],low=bi_col_pal[1],guide = "none")+
#   labs(y="Probability of\nculturing virus",
#        caption="Pickering et al. 2021",
#        #subtitle=expression(paste("VL at 50% = 1.6 x ",10^7," RNA copies/ml (Ct = 21.8)"))
#        )+
#   plotting_theme


(log_plot|inf_plot|culture_plot)/(auc_plot|days_inf_plot)+plot_annotation(tag_levels = "A")

#((log_plot|inf_plot)/(violin_plot|auc_plot))+plot_annotation(tag_levels = "A")

ggsave("results/Fig3 - viral_load_combined_plot.png",dpi=600,width=300,height=150,units="mm",bg="white")
ggsave("results/Fig3 - viral_load_combined_plot.pdf",dpi=600,width=300,height=150,units="mm",bg="white")

