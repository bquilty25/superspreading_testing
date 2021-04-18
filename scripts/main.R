# Load required packages scripts
pacman::p_load("fitdistrplus","tidyverse","patchwork")

# Load data
inf_curve <- read.csv("data/auc.csv")

contacts_o18 <- read.csv("2020-cov-tracing/data/contact_distributions_o18.csv") %>% 
  pivot_longer(1:3) %>% 
  mutate(age="o18")

contacts_u18 <- read.csv("2020-cov-tracing/data/contact_distributions_u18.csv")%>% 
  pivot_longer(1:3) %>% 
  mutate(age="u18")

contacts <- bind_rows(contacts_o18,contacts_u18)

# estimate parameters and take product
#fitdistrplus::descdist(contacts$value,boot = 1000)
contacts_dist <- fitdist(contacts$value,"nbinom")
contacts_sample <- rnbinom(n=100000,size=contacts_dist$estimate[1],mu=contacts_dist$estimate[2])

contacts_plot <- as_tibble(contacts_sample)%>% 
  mutate(value=ifelse(round(value)>20,">20",round(value)),
         value=fct_relevel(value,c(as.character(seq(0,20)),">20"))) %>% 
  ggplot(aes(x=value,
             y=..prop..,
             group=1))+
  geom_bar()+
  scale_x_discrete(breaks=c(as.character(seq(0,18,by=2)),">20"))+
  labs(x="Number of contacts",y="Probability")+
  theme_minimal()

#fitdistrplus::descdist(inf_curve$tot_auc,boot = 1000)
auc_dist <- fitdistrplus::fitdist(inf_curve$tot_auc,"gamma")
auc_sample <- rgamma(n=100000,shape=auc_dist$estimate[1],rate=auc_dist$estimate[2])

auc_plot <- as_tibble(auc_sample)%>% 
  ggplot(aes(x=value),fill="grey")+
  geom_density()+
  lims(x=c(0,10))+
  labs(x="Area under infectiousness curve",y="Density")+
  theme_minimal()

R_sample <- contacts_sample*auc_sample
(R=fitdist(round(R_sample),"nbinom"))

sec_cases_plot <- as_tibble(R_sample) %>% 
  mutate(value=ifelse(round(value)>20,">20",round(value)),
         value=fct_relevel(value,c(as.character(seq(0,20)),">20"))) %>% 
  ggplot(aes(x=value,
             y=..prop..,
             group=1))+
  geom_bar()+
  labs(x="Number of secondary cases",y="Probability")+
  theme_minimal()+
  annotate(label=paste0("R = ",as.character(round(R$estimate[2],2)),"\n","k = ",as.character(round(R$estimate[1],2))),x="17",y=0.5,geom="text")

(auc_plot+contacts_plot)/
  sec_cases_plot

ggsave("results/plot.png",width=210,height=120,units="mm",dpi=600)
       