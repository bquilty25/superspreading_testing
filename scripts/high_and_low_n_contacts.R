source("scripts/utils.R")


pre_pandemic_high<- contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(period=="BBC Pandemic") %>% 
  summarise.(n=n(),
             over_50=sum(e_all>50),
             over_100=sum(e_all>100),
             over_200=sum(e_all>200),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(over_50:over_200)) %>% 
  mutate.(prop_over=value/n,survey="BBC Pandemic") %>% 
  ggplot()+
  geom_point(aes(x=as.Date("2018-01-01"),y=prop_over,colour=fct_relevel(name,"over_50")),size=2,show.legend = F,pch=18)+
  geom_rect(data=time_periods %>% filter(period=="BBC Pandemic"),aes(xmin=-Inf,xmax=Inf,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  #geom_text(data=time_periods %>% filter(period=="BBC Pandemic"),aes(x=as.Date("2018-01-01"),y=0.03,label=period),vjust=1,size=2,colour="#2E4C6D")+
  facet_wrap(~survey)+
  coord_cartesian(xlim=c(as.Date("2017-12-01"),as.Date("2018-02-01")))+
  scale_x_yearmonth("", breaks="1 year",expand=expansion(),date_labels="%Y")+
  scale_y_continuous("% of participants",labels = percent,limits=c(0,0.03),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="Daily number of contacts", values=met.brewer("Java",3,override.order =T),labels=c("over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  scale_fill_manual("Time period",values=met.brewer("Java",6,override.order = F))+
  plotting_theme

pre_pandemic_low <- contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(period=="BBC Pandemic") %>% 
  summarise.(n=n(),
             zero = sum(e_all==0),
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(zero:over_20)) %>% 
  mutate.(prop_over=value/n,survey="BBC Pandemic")  %>% 
  ggplot()+
  geom_point(aes(x=as.Date("2018-01-01"),y=prop_over,colour=fct_relevel(name,"over_5")),size=2,show.legend = F,pch=18)+
  geom_rect(data=time_periods %>% filter(period=="BBC Pandemic"),aes(xmin=-Inf,xmax=Inf,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  #geom_text(data=time_periods %>% filter(period=="BBC Pandemic"),aes(x=as.Date("2018-01-01"),y=0.75,label=period),vjust=1,size=2,colour="#2E4C6D")+
  facet_wrap(~survey)+
  coord_cartesian(xlim=c(as.Date("2017-12-01"),as.Date("2018-02-01")))+
  scale_x_yearmonth("", breaks="1 year",expand=expansion(),date_labels="%Y")+
  scale_y_continuous("% of participants",labels = percent,limits = c(0,0.75),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="Daily number of contacts", values=met.brewer("Hokusai3",4,override.order =T),labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20"))+
  scale_fill_manual("Time period",values=met.brewer("Java",6,override.order = F))+
  plotting_theme

high_n_plot <-contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(date>"2020-05-17",date<"2021-01-01") %>% 
  summarise.(n=n(),
             
            over_50=sum(e_all>50),
            over_100=sum(e_all>100),
            over_200=sum(e_all>200),
            .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(over_50:over_200)) %>% 
  mutate.(prop_over=value/n,survey="Comix")  %>% 
  ggplot()+
  geom_point(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_50")),size=0.8)+
  geom_line(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_50")),size=0.8)+
  geom_rect(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")),aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")) %>% mutate(period=fct_recode(period, "Lockdown\n2 easing"="Lockdown 2 easing")),aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=0.03,label=period),vjust=1,size=2,colour="#2E4C6D")+
  facet_wrap(~survey)+
  scale_x_yearmonth("", breaks="1 month",expand=expansion(0,0),date_labels="%b '%y")+
  scale_y_continuous("",labels = percent,limits=c(0,0.03),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="", values=met.brewer("Java",3,override.order =T),labels=c("over_50"= "Over 50", "over_100"= "Over 100", "over_200"="Over 200"))+
  scale_fill_manual("Time period",values=met.brewer("Java",6,override.order = F))+
  plotting_theme

low_n_plot <-contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(date>"2020-01-01",date<"2021-01-01") %>% 
  summarise.(n=n(),
             zero = sum(e_all==0),
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(zero:over_20)) %>% 
  mutate.(prop_over=value/n,survey="Comix")  %>% 
  ggplot()+
  geom_point(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_5")),size=0.8)+
  geom_line(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_5")),size=0.8)+
  geom_rect(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")),aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")) %>% mutate(period=fct_recode(period, "Lockdown\n2 easing"="Lockdown 2 easing")),aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=0.75,label=period),vjust=1,size=2,colour="#2E4C6D")+
  facet_wrap(~survey)+
  scale_x_yearmonth("", breaks="1 month",expand=expansion(),date_labels="%b '%y")+
  scale_y_continuous("",labels = percent,limits = c(0,0.75),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="Daily number of contacts", values=met.brewer("Hokusai3",4,override.order =T),labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20"))+
  plotting_theme

((pre_pandemic_low/pre_pandemic_high)|(low_n_plot/high_n_plot))+plot_layout(widths=c(1,10),guides="collect")&theme(legend.position = "bottom")

ggsave("results/high_n_contacts.png",width=300,height=150,dpi=600,units="mm",bg="white")


contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(period!="POLYMOD") %>% 
  summarise.(n=n(),
             zero=sum(e_all==0),
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             over_50=sum(e_all>50),
             over_100=sum(e_all>100),
             over_200=sum(e_all>200),
             .by=c(period,date_start,date_end)) %>% 
  pivot_longer.(c(zero:over_200)) %>% 
  mutate.(prop_over=value/n) %>% 
  select.(-value) %>% 
  pivot_wider.(names_from = name,values_from = prop_over,names_sort = T) %>% 
  select.(period,date_start, date_end, n, zero, over_5, over_10, over_20, over_50, over_100, over_200) %>% 
  mutate.(across.(c(zero:over_200),function(x)formatC(x*100,digits = 1,format="f")),
          date_start=format(date_start,"%d %b %Y"),
          date_end=format(date_end,"%d %b %Y")) %>% 
  htmlTable::htmlTable(rnames = F,header = c("Time period","From", "To", "N","Zero", "Over 5", "Over 10", "Over 20", "Over 50", "Over 100", "Over 200"))

contact_data %>% 
  mutate.(date_start=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date_start),
          date_end=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date_end),
          date=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date),
    date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(period!="POLYMOD") %>% 
  pivot_longer.(c(e_all)) %>% 
  summarise.(n=n(),
             zero=sum(value==0),
             over_20=sum(value>20),
             .by=c(date_yw,date_y,name)) %>% 
  pivot_longer.(c(zero:over_20)) %>% 
  mutate.(prop_over=value/n,survey=ifelse.(date_y>="2020","Pandemic\n(Comix)","Pre-pandemic\n(BBC Pandemic)"),
          survey=fct_rev(survey)) %>% 
  select.(-value) %>% 
  ggplot()+
  geom_point(aes(x=date_yw,y=prop_over,colour=fct_relevel(name.1,"over_20")))+
  geom_line(aes(x=date_yw,y=prop_over,colour=fct_relevel(name.1,"over_20")),size=0.8)+
  geom_rect(data=time_periods %>% filter(period!="POLYMOD"),
            aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_vline(data=time_periods %>% filter(period!="POLYMOD"),
            aes(xintercept=date_start,ymax=Inf,ymin=-Inf),alpha=0.1,fill=bi_col_pal[1])+
  geom_vline(data=time_periods %>% filter(period!="POLYMOD"),
             aes(xintercept=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill=bi_col_pal[1])+
  geom_text(data=time_periods %>% filter(period!="POLYMOD"), #%>% 
              #mutate(period=fct_recode(period, 
              #                         "Lockdown\n2 easing"="Lockdown 2 easing",
              #                         "Lockdown\n3 + schools" = "Lockdown 3 + schools"
              #                         )),
            aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=0.3,label=str_wrap(period,10)),vjust=1,hjust=0.5,size=2,colour="#2E4C6D")+
  facet_grid2(~survey,scales="free",space="free")+
  facetted_pos_scales(x=list(
    scale_x_yearmonth("", breaks="1 year",expand=expansion(mult = 20),date_labels="%Y",limits=c(as.Date("2018-01-01"),as.Date("2018-01-01"))),
    scale_x_yearmonth("", breaks="1 month",expand=expansion(),date_labels="%b '%y",limits=c(as.Date("2020-03-15"),NA))))+
  scale_y_continuous("% of participants",labels = percent,limits = c(0,0.30),expand = expansion(mult = c(0,0.05)))+
  # scale_colour_manual(name="Daily number of contacts", values=met.brewer("Hokusai3",2,override.order =F),labels=c( "over_20"=">20","zero"="0"))+
  scale_colour_manual(name="Daily number of contacts", values=bi_col_pal,labels=c( "over_20"=">20","zero"="0"))+
  plotting_theme


ggsave("results/high_n_contacts2.png",width=350,height=150,dpi=600,units="mm",bg="white")


contact_data %>% 
  filter(period!="POLYMOD") %>%  
  mutate.(date_start=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date_start),
                    date_end=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date_end),
                    date=ifelse(period=="BBC Pandemic",as.Date("2018-01-01"),date),
    pandemic=ifelse(period!="BBC Pandemic",
                          "Pandemic\n(Comix, 2020)",
                          "Pre-pandemic\n(BBC Pandemic, 2017)"), 
    date_yw=yearweek(date),
          date_y=year(date),
          date_ym=yearmonth(date),
          period=if_else(!(str_detect(period,"BBC Pandemic")),
                                          paste("Comix: ",period," (",date_ym,")"),
                                          paste0(period," (",date_ym,")"))) %>% 
  mutate.(n_cat=cut(e_all,
                    breaks=c(-Inf,1,5,10,20,50,100,Inf),
                    labels=c("0","1-5","6-10","11-20","21-50","51-100","100+"),
                    right=F)) %>% 
  ggplot()+
  geom_bar(aes(x=date_yw,fill=fct_rev(n_cat),y=after_stat(count/sum(count))),position = "fill")+
  geom_segment(data=time_periods %>% filter(period!="POLYMOD"),
             aes(x=date_start,xend=date_end,yend=1.1,y=1.1),alpha=0.1,fill=bi_col_pal[1])+
  geom_text(data=time_periods %>% filter(period!="POLYMOD",period!="BBC Pandemic"), #%>% 
            #mutate(period=fct_recode(period, 
            #                         "Lockdown\n2 easing"="Lockdown 2 easing",
            #                         "Lockdown\n3 + schools" = "Lockdown 3 + schools"
            #                         )),
            aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=1.1,label=str_wrap(period,10)),vjust=1,hjust=0.5,size=2,colour="#2E4C6D")+
  scale_fill_met_d(palette="Hokusai1",name="Number of daily reported contacts")+
  rotate_x_text(angle=45)+
  labs(x="",
       y="Proportion of participants")+
  facet_wrap2(~fct_rev(pandemic),scales="free_x")+
  facetted_pos_scales(x=list(
    scale_x_yearmonth("", breaks="1 year",expand=expansion(mult=5),date_labels="%Y",limits=c(as.Date("2018-01-01"),as.Date("2018-01-01"))),
    scale_x_yearmonth("", breaks="1 month",expand=expansion(),date_labels="%b '%y",limits=c(as.Date("2020-03-15"),NA))))+
  force_panelsizes(cols = c(0.2,0.8))+
  plotting_theme+
  guides(
    #reverse color order (higher value on top)
    fill = guide_colorbar(reverse = TRUE))

ggsave("results/high_n_contacts3.png",width=210,height=150,dpi=600,units="mm",bg="white")

