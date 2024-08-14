#create figure 2
source("scripts/utils.r")

colour_pal <- c("#17877b", "#D7402B", "#055a8c", "#daa520", "#20bdcc", "#010f5b", "#d72638")

(line_plot <- contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  summarise.(n=n(),
             #zero = sum(e_all==0),
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             over_50=sum(e_all>50),
             over_100=sum(e_all>100),
             over_200=sum(e_all>200),
             .by=c(period)) %>% 
  pivot_longer.(c(over_5:over_200)) %>% 
  mutate.(name=fct_relevel(name,
                          "over_5",
                          "over_10",
                          "over_20",
                          "over_50",
                          "over_100",
                          "over_200")) %>% 
  mutate.(hi_lo=if_else(extract_numeric(name)>=50,"hi","lo"))  %>% 
  rowwise %>% 
  mutate(tst = list(broom::tidy(prop.test(value, n, conf.level=0.95)))) %>%
  tidyr::unnest(tst) %>% 
  ggplot()+
  # geom_texthline(data=. %>% filter(period=="Pre-pandemic"),
  #                aes(yintercept=estimate*100,
  #                    colour=fct_relevel(name,"over_5"),
  #                    label=paste0("Over ", extract_numeric(name), " pre-pandemic")
  #                    ),
  #                size=2.5,
  #                hjust=1.00,
  #                linetype="dashed",
  #                show.legend=F)+
    geom_line(data=. %>% filter(period%!in%c("POLYMOD")),
               aes(x=period,
                   group=name,
                   y=estimate*100,
                   colour=fct_relevel(name,"over_5")),
               #size=0.5
    )+
  geom_point(data=. %>% filter(period%!in%c("POLYMOD")),
            aes(x=period,
                group=name,
                y=estimate*100,
                colour=fct_relevel(name,"over_5")),
            size=0.5
            )+
  geom_linerange(data=. %>% filter(period%!in%c("POLYMOD")),
            aes(x=period,
                group=name,
                ymin=conf.low*100,
                ymax=conf.high*100,
                #alpha=0.15,
                colour=fct_relevel(name,"over_5"),
                fill=fct_relevel(name,"over_5")),
            #size=0.5,
            show.legend = F)+
  facet_grid2(fct_rev(hi_lo)~.,scales="free_y",axes="all",remove_labels = "x")+
  scale_y_continuous("Percentage of participants (%)"#,expand = expansion(mult = c(0,0.05))
                     )+
  scale_x_discrete(expand = expansion(0,c(0.5,0.5)))+
  guides(colour = guide_legend(nrow = 1))+
    scale_colour_brewer(name="Reported daily contacts", 
                        palette = "Set2",
                        labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
                                 "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
    scale_fill_brewer(name="Daily number of contacts", 
                      palette = "Set2",
                      labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
                               "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  # scale_colour_manual(name="Daily number of contacts", 
  #                     values=rev(colour_pal),
  #                     labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
  #                              "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  # scale_fill_manual(name="Daily number of contacts", 
  #                     values=rev(colour_pal),
  #                     labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
  #                              "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  plotting_theme+
  labs(x="")+
  theme(strip.text = element_blank())+
  ggpubr::rotate_x_text(angle=45)
)

ggsave("results/high_n_contacts4.png",width=210,height=150,dpi=600,units="mm",bg="white")

dot_plot <-   contact_data %>% 
  filter.(period%in%c("Pre-pandemic","1st lockdown","School reopening")) %>% 
  pivot_longer.(cols=c(e_home,e_other,e_all)) %>%
  mutate.(ecdf_x=ecdf(value)(value),.by=c(name,period)) %>% 
  ggplot()+
  geom_point(aes(x=value,y=1-ecdf_x,colour=period),alpha=0.5)+
  facet_wrap2(~name,labeller = labeller(name=c("e_all"="All contacts",
                                               "e_home"="Household contacts",
                                               "e_other"="Out of household contacts")),
              axes="all")+
  scale_x_continuous(trans="pseudo_log",breaks = c(0,1,10,100,1000),expand = expansion(0,0))+
  scale_y_continuous(trans="log10",labels=label_percent())+
  scale_colour_manual(name="Time period",values=colour_pal)+
  plotting_theme+
  labs(x="Reported daily contacts",y=str_wrap("Percentage of participants reporting at least X contacts (%)",35))

# dot_plot <-   contact_data %>% 
#   filter.(period%in%c("Pre-pandemic","1st lockdown","School reopening")) %>% 
#   pivot_longer(cols=c(e_home,e_other,e_all)) %>%
#   summarise(n=n(),
#             .by=c(value,name,period))%>% 
#   mutate(prop=n/sum(n),.by=name,period) %>% 
#   ggplot()+
#   geom_point(aes(x=value,y=prop*100,colour=period),alpha=0.3)+
#   facet_grid2(~name,labeller = labeller(name=c("e_all"="All contacts",
#                                                "e_home"="Household contacts",
#                                                "e_other"="Out of household contacts")),
#               axes="all")+
#   scale_x_continuous(trans="pseudo_log",breaks = c(0,1,10,100,1000),expand = expansion(0,0))+
#   scale_y_continuous(trans="log10",labels=label_number(accuracy = 0.01))+
#   scale_colour_manual(name="Time period",values=colour_pal)+
#   plotting_theme+
#   labs(x="Number of daily reported contacts",y="Percentage of participants (%)")


ggsave("results/contacts_ecdf.png",width=210,height=120,dpi=600,units="mm",bg="white")

nbinom_plot <- contact_data %>%  
  pivot_longer.(c(e_all,e_home,e_other),names_to = "contact_type") %>% 
  drop_na.(value) %>% 
  summarise.(dist_means=list(fitdist(value,"nbinom")$estimate %>% enframe()),.by=c(period,contact_type)) %>% 
  unnest.(dist_means) %>% 
  filter.(period!="POLYMOD",
          contact_type=="e_all") %>% 
  ggplot(aes(y=value,x=period))+
  geom_point(aes(group=name,colour=name))+
  geom_line(aes(group=name,colour=name))+
  scale_colour_manual(values = bi_col_pal,guide="none")+
  scale_linetype_manual(values=c("dashed",NA),guide="none")+
  lims(y=c(0,NA))+
  labs(y="",
       x="",
  )+
  plotting_theme+
  theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  facet_grid2(name~.,switch="y",scales="free_y",labeller=labeller(name=c("mu"="Mean\nreported daily contacts","size"="Overdispersion (k) of\nreported daily contacts")),axes="all",remove_labels = "x")+
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits=c(0,NA),expand = expansion(c(0,0.1))),scale_y_log10(limits=c(0.25,3))))

line_plot/dot_plot/nbinom_plot+plot_annotation(tag_levels = "A")+plot_layout(heights=c(1.5,1,1.5))
ggsave("results/Fig2 - contacts.png",width=210,height=320,dpi=600,units="mm",bg="white")
ggsave("results/Fig2 - contacts.pdf",width=210,height=320,dpi=600,units="mm",bg="white")

# for supplement
dot_plot_all <- contact_data %>% 
  filter(period%!in%c("POLYMOD")) %>% 
  pivot_longer(cols=c(e_home,e_other)) %>%
  summarise(n=n(),
            .by=c(value,name,period))%>% 
  mutate(prop=n/sum(n),.by=name,period) %>% 
  ggplot()+
  geom_point(aes(x=value,y=prop*100),alpha=0.3)+
  facet_grid2(period~name,labeller = labeller(name=c("e_all"="All contacts",
                                               "e_home"="Household contacts",
                                               "e_other"="Out of household contacts")),
              axes="all")+
  scale_x_continuous(trans="pseudo_log",breaks = c(0,1,10,100,1000),expand = expansion(0,0))+
  scale_y_continuous(trans="log10",labels=label_number(accuracy = 0.01))+
  scale_colour_manual(name="Time period",values=colour_pal)+
  plotting_theme+
  labs(x="Number of daily reported contacts",y="Percentage of participants (%)")

    
dot_plot_all <- contact_data %>% 
    filter.(period%!in%c("POLYMOD")) %>% 
    pivot_longer.(cols=c(e_home,e_other)) %>%
    mutate.(ecdf_x=ecdf(value)(value),.by=c(name,period)) %>% 
    ggplot()+
    geom_point(aes(x=value,y=1-ecdf_x,group=period),alpha=0.5)+
    ggh4x::facet_grid2(period~name,labeller = labeller(name=c("e_all"="All contacts",
                                                 "e_home"="Household contacts",
                                                 "e_other"="Out of household contacts")
                                                 ),
                       axes= "all"
                )+
    scale_x_continuous(trans="pseudo_log",breaks = c(0,1,10,100,1000),expand = expansion(0,0))+
    scale_y_continuous(trans="log10",labels=label_percent())+
    scale_colour_manual(name="Time period",values=colour_pal)+
    plotting_theme+
    labs(x="Number of reported daily contacts",y=str_wrap("Percentage of participants reporting at least X contacts (%)",35))
  
  ggsave("results/contacts_ecdf_all1.png",width=210,height=350,dpi=600,units="mm",bg="white")
  ggsave("results/contacts_ecdf_all1.pdf",width=210,height=350,dpi=600,units="mm",bg="white")
  
  