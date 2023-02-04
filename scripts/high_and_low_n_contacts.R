cols <- c("#17877b", "#D7402B", "#055a8c", "#daa520", "#20bdcc", "#010f5b", "#d72638")


design <- 
"AA
BC"

dist_plot <- contact_data %>% 
  filter(period%!in%c("POLYMOD")) %>%
  pivot_longer(cols=c(e_home,e_other,e_all))%>% 
  ggplot(aes(y=value,
             x=period,
             fill=fct_rev(name)))+
  scale_fill_manual(values=cols)+
  stat_histinterval(slab_type="histogram",
                    side="both",
                    normalize="groups",
                    size=0.25,
            breaks=function(x){seq(min(x),max(x),by=1)},
            alpha=0.75
            #breaks=c(-Inf,0,1,2,3,4,5,6,7,8,9,20,50,Inf)
            )+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,20,50,100))+
  #scale_y_continuous(labels = percent)+
  coord_cartesian(ylim = c(NA,50))+
  facet_manual(~name,
              labeller = labeller(name=c("e_all"="All contacts",
                                         "e_home"="Household contacts",
                                         "e_other"="Out of household contacts")),
              scales="free_y",
              axes="all",
              design=design)+
  guides(fill="none")+
  labs(y="Number of reported daily contacts",x="")+
  plotting_theme+
  ggpubr::rotate_x_text(angle=45)

ggsave("results/contact_dist.png",width=210,height=175,dpi=600,units="mm",bg="white")


line_plot <- contact_data %>% 
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
  mutate(name=fct_relevel(name,
                          "over_5",
                          "over_10",
                          "over_20",
                          "over_50",
                          "over_100",
                          "over_200")) %>% 
  mutate.(prop_over=value/n,
          #date_character=factor(date_yw,ordered = T),
          #survey=ifelse.(date_y>="2020","Pandemic\n(Comix)","Pre-pandemic\n(Pre-pandemic)"),
          #survey=fct_rev(survey),
          hi_lo=if_else(extract_numeric(name)>=50,"hi","lo"))  %>% 
  rowwise %>% 
  mutate(tst = list(broom::tidy(prop.test(value, n, conf.level=0.95)))) %>%
  tidyr::unnest(tst) %>% 
  ggplot()+
  geom_texthline(data=. %>% filter(period=="Pre-pandemic"),
                 aes(yintercept=prop_over,
                     colour=fct_relevel(name,"over_5"),
                     label=paste0("Over ", extract_numeric(name), " pre-pandemic")
                     ),
                 size=2.5,
                 hjust=1.00,
                 linetype="dashed",
                 show.legend=F)+
  geom_line(data=. %>% filter(period%!in%c("Pre-pandemic","POLYMOD")),
            aes(x=period,
                group=name,
                y=prop_over,
                colour=fct_relevel(name,"over_5")),
            size=0.5)+
  geom_ribbon(data=. %>% filter(period%!in%c("Pre-pandemic","POLYMOD")),
            aes(x=period,
                group=name,
                ymin=conf.low,
                ymax=conf.high,
                alpha=0.25,
                fill=fct_relevel(name,"over_5")),
            size=0.5,show.legend = F)+
  # geom_rect(data=time_periods %>% filter(period!="POLYMOD"),
  #           aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  facet_grid2(fct_rev(hi_lo)~.,scales="free_y",axes="all")+
  #scale_x_yearmonth("", breaks="1 month",expand=expansion(),date_labels="%b '%y",limits=c(as.Date("2020-03-15"),NA))+
  scale_y_continuous("% of participants",labels = percent#,expand = expansion(mult = c(0,0.05))
                     )+
  scale_x_discrete(expand = expansion(0,c(0.5,1.5)))+
  guides(colour = guide_legend(nrow = 1))+
  scale_colour_manual(name="Daily number of contacts", 
                      values=cols,
                      labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
                               "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  scale_fill_manual(name="Daily number of contacts", 
                      values=cols,
                      labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20",
                               "over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  plotting_theme+
  labs(x="")+
  theme(strip.text = element_blank())+
  ggpubr::rotate_x_text(angle=45)

ggsave("results/high_n_contacts4.png",width=210,height=150,dpi=600,units="mm",bg="white")

dist_plot/line_plot+plot_annotation(tag_levels = "A")

ggsave("results/contacts_combined.png",width=210,height=300,dpi=600,units="mm",bg="white")
