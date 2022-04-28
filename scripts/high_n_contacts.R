source("scripts/utils.R")


pre_pandemic_high <- contact_data %>% 
  mutate.(date_yw=yearweek(date),
          date_y=year(date)) %>% 
  filter.(period=="BBC Pandemic") %>% 
  summarise.(n=n(),
             over_50=sum(e_all>50),
             over_100=sum(e_all>100),
             over_200=sum(e_all>200),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(over_50:over_200)) %>% 
  mutate.(prop_over=value/n)  %>% 
  ggplot()+
  geom_point(aes(x=as.Date("2018-01-01"),y=prop_over,colour=fct_relevel(name,"over_50")),size=2,show.legend = F,pch=18)+
  geom_rect(data=time_periods %>% filter(period=="BBC Pandemic"),aes(xmin=-Inf,xmax=Inf,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(period=="BBC Pandemic"),aes(x=as.Date("2018-01-01"),y=0.03,label=period),vjust=1,size=2,colour="#2E4C6D")+
  #facet_wrap(~date_y,scales="free_x")+
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
             
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(over_5:over_20)) %>% 
  mutate.(prop_over=value/n)  %>% 
  ggplot()+
  geom_point(aes(x=as.Date("2018-01-01"),y=prop_over,colour=fct_relevel(name,"over_5")),size=2,show.legend = F,pch=18)+
  geom_rect(data=time_periods %>% filter(period=="BBC Pandemic"),aes(xmin=-Inf,xmax=Inf,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(period=="BBC Pandemic"),aes(x=as.Date("2018-01-01"),y=0.75,label=period),vjust=1,size=2,colour="#2E4C6D")+
  #facet_wrap(~date_y,scales="free_x")+
  coord_cartesian(xlim=c(as.Date("2017-12-01"),as.Date("2018-02-01")))+
  scale_x_yearmonth("", breaks="1 year",expand=expansion(),date_labels="%Y")+
  scale_y_continuous("% of participants",labels = percent,limits = c(0,0.75),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="Daily number of contacts", values=met.brewer("Hokusai3",3,override.order =T),labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20"))+
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
  mutate.(prop_over=value/n)  %>% 
  ggplot()+
  geom_point(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_50")),size=0.8)+
  geom_line(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_50")),size=0.8)+
  geom_rect(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")),aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")) %>% mutate(period=fct_recode(period, "Lockdown\n2 easing"="Lockdown 2 easing")),aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=0.03,label=period),vjust=1,size=2,colour="#2E4C6D")+
  #facet_wrap(~date_y,scales="free_x")+
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
             
             over_5=sum(e_all>5),
             over_10=sum(e_all>10),
             over_20=sum(e_all>20),
             .by=c(date_yw,date_y)) %>% 
  pivot_longer.(c(over_5:over_20)) %>% 
  mutate.(prop_over=value/n)  %>% 
  ggplot()+
  geom_point(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_50")),size=0.8)+
  geom_line(aes(x=date_yw,y=prop_over,colour=fct_relevel(name,"over_5")),size=0.8)+
  geom_rect(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")),aes(xmin=date_start,xmax=date_end,ymax=Inf,ymin=-Inf),alpha=0.1,fill="#FC997C")+
  geom_text(data=time_periods %>% filter(date_start>as.Date("2020-01-01"),date_start<as.Date("2021-01-01")) %>% mutate(period=fct_recode(period, "Lockdown\n2 easing"="Lockdown 2 easing")),aes(x=(as.numeric(date_end)+as.numeric(date_start))/2,y=0.75,label=period),vjust=1,size=2,colour="#2E4C6D")+
  #facet_wrap(~date_y,scales="free_x")+
  scale_x_yearmonth("", breaks="1 month",expand=expansion(),date_labels="%b '%y")+
  scale_y_continuous("",labels = percent,limits = c(0,0.75),expand = expansion(mult = c(0,0.05)))+
  scale_colour_manual(name="Daily number of contacts", values=met.brewer("Hokusai3",3,override.order =T),labels=c("over_5"= "Over 5", "over_10"= "Over 10", "over_20"="Over 20"))+
  plotting_theme

((pre_pandemic_low/pre_pandemic_high)|(low_n_plot/high_n_plot))+plot_layout(widths=c(1,10),guides="collect")&theme(legend.position = "bottom")

ggsave("results/high_n_contacts.png",width=220,height=150,dpi=600,units="mm",bg="white")

