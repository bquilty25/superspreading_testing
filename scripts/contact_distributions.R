
# Run code and estimate mean and k of contact distributions (SLOW)
if(!file.exists("results/contact_dat_params.qs")){
tic()
contact_data_dists <- contact_data %>% 
  #pivot_longer.(c(e_all,e_home,e_other)) %>% 
  #mutate.(name=fct_relevel(name,"e_all","e_home","e_other")) %>% 
  #drop_na.(value) %>% 
  #nest_by.(period) %>%
  summarise.(#dists=map.(.x=data,.f= . %>% pull.(e_all) %>% fitdist("nbinom")),
          dist_means=list(fitdist(e_all,"nbinom")$estimate %>% t()),.by=period
         #boot_dist=map.(.x=dists, ~bootdist(f =.,bootmethod = "nonparam",parallel="snow",ncpus=4)$CI %>% 
                          #as.data.frame() %>% 
                          #rownames_to_column)
         )

qs::qsave(contact_data_dists,"results/contact_dat_params.qs")

toc()
} else {
  contact_data_dists <- qs::qread("results/contact_dat_params.qs")
}


         #params=map.(dists,~c(mu=.x$estimate[[2]],k=.x$estimate[[1]]))
         #) %>% 
  #unnest_wider(params) 

contact_data_dists %>% 
  unnest.(dist_means) %>% 
  filter.(period!="POLYMOD") %>% 
  pivot_longer.(c(size,mu)) %>% 
  ggplot(aes(y=value,x=period))+
  geom_point(aes(group=name,colour=name))+
  geom_line(aes(group=name,colour=name))+
  facet_rep_grid(name~.,scales="free_y",switch="y",labeller=labeller(name=c("mu"="Mean","size"="k")))+
  scale_colour_manual(values = bi_col_pal,guide="none")+
  scale_linetype_manual(values=c("dashed",NA),guide="none")+
  lims(y=c(0,NA))+
  labs(y="Mean parameter value",
       x="Time period",
       #title="The relative impact of lateral flow testing is greatest when individuals have lots of contacts,\nwhen uptake is high, and when testing is frequent")+
  )+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))

ggsave("results/contact_dists.png",width=210,height=150,dpi=600,units="mm",bg="white")


contact_data %>% 
  pivot_longer(c(e_all,e_home,e_other)) %>% 
  drop_na(value) %>% 
  mutate(name=fct_relevel(name,"e_all","e_home","e_other"),
         value=factor(ifelse(value>=20,"\u2265 20",value)),
         value=fct_relevel(value,c(as.character(seq(0,19)),"\u2265 20"))) %>%
  ggplot(aes(x = value,
             y = ..prop..,
             group = 1))+
  geom_bar(width = 0.8)+
  geom_label(data = contact_data_dists %>% unnest(boot_dist) %>% select(-c(`2.5%`,`97.5%`)) %>% pivot_wider(values_from = Median,names_from = rowname),
             aes(label=paste0("Mean = ", sprintf("%.2f",mu),"\nk = ",sprintf("%.2f",size)),
                 x="\u2265 20",
                 y=0.7),
             colour = "white",
             fill = rgb(0,0,0,0.4),
             hjust = 1)+
  scale_fill_brewer(type="qual",guide=F,palette = "Set2")+
  labs(x="Number of daily contacts",y="Probability")+
  theme_minimal()+
  theme(axis.line.x.bottom = element_line(),
        axis.ticks.x.bottom = element_line(),
        axis.line.y.left = element_line())+
  facet_rep_grid(name~period, scales='free_y', switch="y",
                 #labeller=labeller(period=c("pre"="Pre-pandemic (BBC 2018)",
                 #                                "aug_sept"="Relaxed (Comix Aug/Sept 2020)",
                 #                                "jan_feb"="Lockdown (Comix Jan/Feb 2021)"),
                 #                   name=c("e_all"="All","e_home"="Home","e_work_school"="Work/School","e_other"="Other"),
                 #                  prop_self_iso=function(x){scales::percent(as.numeric(x))})
  ) + 
  scale_x_discrete(breaks=c(as.character(seq(0,18,by=2)),"\u2265 20"),expand = expansion(add=0.7))+
  scale_y_continuous(limits=c(0,1),expand=expansion(add=0.01))

ggsave("results/contacts.png",width=12,height=7,units="in",dpi=400,scale=0.8)
