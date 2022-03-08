
make_trajectories(variant=c("wild"),n_sims = 10,n_cases=5) %>% 
  mutate(pred = pmap(.l = list(m=m,start=start,end=end), 
                     ~with(list(...),data.frame(x = seq(start, end, length.out = 31)) %>% 
                             mutate(ct=m(x),
                                    test_p=stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct))))),
         onset_ct=pmap_dbl(.l=list(m=m,onset_t=onset_t),
                           ~with(list(...),m(onset_t))),
         onset_p=stats::predict(innova_mod, type = "response", newdata = data.frame(ct = onset_ct)))%>% 
  unnest(pred) %>% 
  ggplot()+
  geom_line(aes(x=x,
                y=ct,
                group=interaction(sim,idx,variant),
                colour=test_p),
            alpha=0.5
  )+
  geom_point(data=. %>% filter(type=="symptomatic"),
             aes(x=onset_t,y=onset_ct,group=interaction(sim,idx,variant),colour=onset_p),pch=16,alpha=0.5)+
  scale_colour_viridis_c(name="Probability of culture (infectiousness)",option="mako",begin=0.2,end=0.8)+
  scale_x_continuous(name="Days since exposure",breaks = breaks_width(1))+
  scale_y_reverse(# Features of the first axis
    limits=c(40,NA),
    name = "Cycle\nthreshold")+
  coord_cartesian(ylim=c(40,NA))+
  plotting_theme+
  theme( strip.background = element_blank(),
         panel.grid.minor.x = element_blank(),
         panel.grid.major.x = element_blank(),
         #strip.text.x = element_blank(),
         #axis.title.x = element_blank()
  )+
  facet_wrap(~type,nrow = 2,labeller=labeller(type=capitalize))


ggsave("results/piecewise_linear.png",dpi=600,width=210,height=150,units="mm")
