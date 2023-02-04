# Load required packages scripts
pacman::p_load(
  "qs",
  "ggdist",
  "fitdistrplus",
  "EnvStats",
  "tidyverse",
  "patchwork",
  "here",
  "rriskDistributions",
  "rms",
  "DescTools",
  "MESS",
  "lubridate",
  "lemon",
  "boot",
  "furrr",
  "data.table",
  "ggtext",
  "fst",
  "extraDistr",
  "emdbook",
  "colorspace",
  "fuzzyjoin",
  "qs",
  "ggpubr",
  "bench",
  "tictoc",
  "naniar",
  "scales",
  "ggforce",
  "RGeode",
  "tsibble",
  "MetBrewer",
  "ggrepel",
  "ggh4x",
  "lemon",
  "geomtextpath"
)

if(packageVersion("tidytable")!="0.8.0"){
remotes::install_version("tidytable", version = "0.8")
}else{  
library(tidytable)
}

seed <- 1000
set.seed(seed)

# plotting options
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

`%!in%` = Negate(`%in%`)

plotting_theme <- theme_minimal(base_family = "Lato")+
  theme(axis.ticks = element_line(colour="#2E4C6D"),
        axis.title = element_text(colour="#2E4C6D"),
        axis.text = element_text(colour="#2E4C6D"),
        strip.text = element_text(colour="#2E4C6D"),
        axis.line.x = element_line(colour="#2E4C6D"),
        axis.line.y = element_line(colour="#2E4C6D"),
        #panel.border = element_rect(fill=NA,colour="#2E4C6D"),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.placement = "outside",
        axis.line = element_line(colour="#2E4C6D"),
        line = element_line(colour="#2E4C6D"),
        text = element_text(colour="#2E4C6D",family = "Lato"))

bi_col_pal <- c("#396EB0","#FC997C")
tri_col_pal <- c("#396EB0","#DADDFC","#FC997C")
quad_col_pal <- c("#2E4C6D","#396EB0","#DADDFC","#FC997C")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

### Load viral load data ----
kissler_dat <- read_csv("CtTrajectories_B117/output/shared_params_df.csv") %>% 
  select("alpha_peakvl" = dpmeanB, 
         "wild_peakvl" = dpmeanW,
         "alpha_prolif" = wpmeanB, 
         "wild_prolif" = wpmeanW, 
         "alpha_clear" = wrmeanB, 
         "wild_clear" = wrmeanW,
         "peakvl_sd" = dpsd,  
         "prolif_sd" = wpsd,  
         "clear_sd" = wrsd)

kissler_dat_means <- kissler_dat %>% 
  select(-c(peakvl_sd:clear_sd)) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(mean=median(value)) %>% 
  separate(name,sep = "_",into=c("variant","param"))

kissler_dat_sd <- kissler_dat %>% 
  select(c(peakvl_sd:clear_sd)) %>% 
  pivot_longer(everything(),names_to = "param") %>% 
  group_by(param) %>% 
  summarise(sd=median(value)) %>% 
  mutate(param=str_extract(param, "[^_]+"))

kissler_dat_est <- kissler_dat_means %>% 
  left_join(kissler_dat_sd) %>% 
  pivot_wider(names_from=param,values_from = c(mean,sd))

hay_dat <- read_csv(file="data/shared_params_df.csv") %>% 
  select("omicron_peakvl" = dpmeanB_trans, 
         "delta_peakvl" = dpmeanW_trans,
         "omicron_prolif" = wpmeanB_trans, 
         "delta_prolif" = wpmeanW_trans, 
         "omicron_clear" = wrmeanB_trans, 
         "delta_clear" = wrmeanW_trans,
         "peakvl_sd" = dpsd,  
         "prolif_sd" = wpsd,  
         "clear_sd" = wrsd)

hay_dat_means <- hay_dat %>% 
  select(-c(peakvl_sd:clear_sd)) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(mean=median(value)) %>% 
  separate(name,sep = "_",into=c("variant","param"))

hay_dat_sd <- hay_dat %>% 
  select(c(peakvl_sd:clear_sd)) %>% 
  pivot_longer(everything(),names_to = "param") %>% 
  group_by(param) %>% 
  summarise(sd=median(value)) %>% 
  mutate(param=str_extract(param, "[^_]+"))

hay_dat_est <- hay_dat_means %>% 
  left_join(hay_dat_sd) %>% 
  pivot_wider(names_from=param,values_from = c(mean,sd))

vl_params <- bind_rows(kissler_dat_est,hay_dat_est) %>% 
  mutate(variant=as_factor(variant))

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
  out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
  return(out) 
}

#Define time periods of interest
time_periods <- tribble(~idx,~period,~date_start,~date_end,
                        -1, "POLYMOD",             as_date("01/01/2008",format="%d/%m/%Y"), as_date("01/01/2008",format="%d/%m/%Y"),
                        0, "Pre-pandemic",         as_date("01/09/2017",format="%d/%m/%Y"), as_date("01/12/2018",format="%d/%m/%Y"),
                        1, "1st Lockdown",           as_date("23/03/2020",format="%d/%m/%Y"), as_date("03/06/2020",format="%d/%m/%Y"),
                        2, "1st Lockdown easing",    as_date("04/06/2020",format="%d/%m/%Y"), as_date("29/07/2020",format="%d/%m/%Y"),
                        3, "Relaxed restrictions", as_date("30/07/2020",format="%d/%m/%Y"), as_date("03/09/2020",format="%d/%m/%Y"),
                        4, "School reopening",     as_date("04/09/2020",format="%d/%m/%Y"), as_date("24/10/2020",format="%d/%m/%Y"),
                        5, "2nd Lockdown",           as_date("05/11/2020",format="%d/%m/%Y"), as_date("02/12/2020",format="%d/%m/%Y"),
                        6, "2nd Lockdown easing",    as_date("03/12/2020",format="%d/%m/%Y"), as_date("19/12/2020",format="%d/%m/%Y"),
                        7, "3rd Lockdown",           as_date("05/01/2021",format="%d/%m/%Y"), as_date("07/03/2021",format="%d/%m/%Y"),
                        8, "3rd Lockdown + schools", as_date("08/03/2021",format="%d/%m/%Y"), as_date("31/03/2021",format="%d/%m/%Y"),
                        9, "Step 2 + schools",     as_date("16/04/2021",format="%d/%m/%Y"), as_date("16/05/2021",format="%d/%m/%Y")) %>% 
  mutate(period=factor(period,levels = period))

#### Load contact data ----
contacts_polymod <- 
  read.csv(here("data","2008_Mossong_POLYMOD_contact_common.csv")) %>% 
  pivot_longer(cols=c(cnt_home,cnt_school,cnt_work,cnt_transport,cnt_leisure,cnt_otherplace)) %>% 
  filter(value) %>% 
  select(-value) %>% 
  count(name,part_id) %>% 
  pivot_wider(names_from = name,values_from = n) %>% 
  mutate(e_home=cnt_home,
         e_other=rowSums(across(c(cnt_work,cnt_school,cnt_transport,cnt_leisure,cnt_otherplace)),na.rm=T)) %>% 
  select(part_id,e_home,e_other) %>% 
  complete(part_id=full_seq(part_id,1),fill=list(e_home=0,e_other=0)) %>% 
  mutate(date=as_date("01/01/2008",format="%d/%m/%Y"))

contacts_bbc_o18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_o18.csv")) 

contacts_bbc_u18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_u18.csv")) 

contacts_bbc <- bind_rows(contacts_bbc_o18, contacts_bbc_u18) %>%
  mutate(date=as_date("01/09/2017",format="%d/%m/%Y"),
         e_school=0)

# Clean and load Comix data
source("scripts/comix_clean.R")

#summarise number of contacts
contact_data <- contacts_bbc %>% 
  mutate(e_other=rowSums(across(c(e_work,e_school,e_other)),na.rm = T)) %>% 
  select(date,e_home,e_other) %>% 
  bind_rows(contacts_polymod) %>% 
  mutate(e_all = rowSums(across(c(e_home,e_other)),na.rm = T),
         date=as_date(date),
         part_id=as.character(part_id))%>%
  bind_rows(contacts_comix) %>% 
  fuzzyjoin::fuzzy_inner_join(time_periods,
                              by=c("date"="date_start","date"="date_end"),
                              match_fun=list(`>=`,`<=`))

#### impute out of HH values > 250 for Pre-pandemic by fitting distribution to values from non-lockdown periods ----

# Calculate proportion over 250 by time period
contact_data %>%
  filter(period %in% c("Relaxed restrictions","School reopening","Step 2 + schools")) %>% 
  summarise.(n=n(),over_250=sum(e_other>=250)) %>% 
  mutate.(prop=over_250/n)

#0.00160 or 0.16% over 250

# assume distribution of high contacts is exponential and fit distribution
dist_over_250 <- contact_data %>% 
  filter(period %in% c("Relaxed restrictions","School reopening","Step 2 + schools")) %>% 
  filter(e_other>=250)  %>% 
  pull(e_other) %>% 
  fitdistr(.,"exponential")

#simulate individuals with high numbers of contacts for Pre-pandemic
dat_append <- data.frame(e_other=round(rexptr(n=0.0016*1.0016*nrow(contact_data %>% 
                                                                     filter(period=="Pre-pandemic")),
                                              lambda=dist_over_250$estimate[1],
                                              range = c(250,Inf))),
                         e_home=sample(size=0.0016*1.0016*nrow(contact_data %>% 
                                                                 filter(period=="Pre-pandemic")),
                                       x=contact_data %>% 
                                         filter(period=="Pre-pandemic") %>% 
                                         pull(e_home))) %>% 
  mutate(e_all=e_home+e_other,period="Pre-pandemic",idx=3)

#append to data
contact_data_adjusted <- contact_data %>% bind_rows(dat_append)


##### KCL ANALYSIS ----
pickering <- readxl::read_xlsx(here::here("data","pickering_dat.xlsx")) %>% 
  select(-c(`Viral Growth`,...7,...8)) %>% 
  rename("culture"=...6) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)ifelse(x=="ND",NA,x)) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)case_when(x%in%c(0.5,1,2)~1,
                                         is.na(x)~NA_real_,
                                         TRUE~0)) %>% 
  mutate(id=row_number()) %>% 
  rename(ct=`Ct N1`) %>% 
  mutate(vl=(-(ct-44.34)/3.134))

innova_mod <- glm(Innova~vl,
                  data=pickering,
                  family="binomial") 

innova_higher_mod <- glm(Innova~vl,
                         data=pickering %>% 
                           mutate(vl=vl+2.5),family="binomial") 

#sensitivity analysis on test probability 
test_model_choice <- function(boolean){
  if(boolean){
    innova_higher_mod
  }else{
    innova_mod
  }
}

culture_mod <- glm(culture~vl,data=pickering,family="binomial") 

#sensitivity analysis for infectiousness as culture prob or lft prob
inf_model_choice <- function(boolean){
  #browser()
  if(boolean){
    innova_mod
  }else{
    culture_mod
  }
}

# Historical Rt estimates
rt <- readxl::read_excel("data/221123_R_and_growth_rate_time_series_for_publication_v1.0.xlsx",
                         range = "Table1_-_R!B10:D52",
                         col_names = c("date", "lower", "upper"))

rt_by_time_period <- rt %>% filter(date<as.Date("2021-01-01")) %>%
  fuzzyjoin::fuzzy_inner_join(time_periods,
                              by=c("date"="date_start","date"="date_end"),
                              match_fun=list(`>=`,`<=`)) %>%
  group_by(across(-c(date,lower,upper))) %>% 
  summarise(lower=mean(lower),upper=mean(upper)) %>% 
  ungroup() %>% 
  filter(period!="Lockdown 1") #very minimal overlap with lockdown 1 period (Rt starts 29/5/2020)

# Create viral load trajectories for a given number of sims
make_trajectories <- function(
    n_sims = 100,
    asymp_parms = asymp_fraction,
    variant_info, 
    browsing = FALSE
){
  
  if (browsing) browser()
  
  set.seed(seed)
  #simulate CT trajectories
  
  inf <- rbbinom(n = n_sims,
                 size=1,
                 alpha = asymp_parms$shape1,
                 beta  = asymp_parms$shape2) %>%
    as_tidytable() %>%
    rename.("asymptomatic"=x) %>%
    mutate.(sim=row_number.(),
            asymptomatic=as.logical(asymptomatic))
  
  traj <- inf %>% 
    crossing.(start=0) %>% 
    crossing(variant_info) %>% 
    mutate.(
      prolif=case_when.(heterogen_vl~rnormTrunc(n=n(),mean=mean_prolif,
                                                sd=sd_prolif,min = 1,max=14),
                        TRUE~median(rnormTrunc(n=n(),mean=mean_prolif,
                                               sd=sd_prolif,min = 1,max=14))),
      clear=case_when.(heterogen_vl~rnormTrunc(n=n(), mean=mean_clear, 
                                               sd=sd_clear, min = 1,max=30),
                       TRUE~median(rnormTrunc(n=n(), mean=mean_clear, 
                                              sd=sd_clear, min = 1,max=30))),
      end=prolif+clear,
      onset_t=prolif+rnorm(n=n(),mean = 2,sd=1.5)
    ) %>%
    select.(-c(mean_prolif, sd_prolif, mean_clear, sd_clear,clear)) %>%
    pivot_longer.(cols = -c(sim,variant,onset_t, asymptomatic, heterogen_vl,
                            mean_peakvl,sd_peakvl),
                  values_to = "x") %>%
    mutate.(y=case_when.(name=="start" ~ 40,
                        name=="end"   ~ 40,
                        name=="prolif"~case_when.(heterogen_vl~rnormTrunc(n=n(),
                                                                          mean=mean_peakvl,
                                                                          sd=sd_peakvl,min=0,max=40),
                                                  TRUE~median(rnormTrunc(n=n(),
                                                                         mean=mean_peakvl,
                                                                         sd=sd_peakvl,min=0,max=40))))) %>% 
    select.(-c(mean_peakvl,sd_peakvl))
  
  
  models <- traj %>%
    nest.(data = -c(sim,variant,onset_t,asymptomatic,heterogen_vl)) %>%  
    mutate.(
      # Perform approxfun on each set of points
      m  = map.(data, ~approxfun(x=.x$x,y=.x$y))) 
  
  #cannot pivot wider with "m" column - extract and rejoin
  x_model <- models %>% 
    select.(-data)
  
  models <- models %>% 
    select.(-m) %>% 
    unnest.(data,.drop=F) %>%  
    select.(-c(y)) %>% 
    pivot_wider.(names_from=name,values_from = x) %>% 
    left_join.(x_model) %>%
    select.(c(sim, variant, heterogen_vl, onset_t, prolif, start, end, m)) %>% 
    arrange.(sim)
  
}

inf_curve_func <- function(m,start=0,end=30,interval=1,trunc_t){
  #browser()
  x <- tidytable(t=seq(start,end,by=interval)) %>% 
    mutate.(ct=m(t),
            vl=convert_Ct_logGEML(ct))
  
  return(infectiousness=x)
}

calc_sensitivity <- function(model, x){
  #browser()
  if(!is.na(x)){
    s <- model(x)
  } else {
    s <- NA_real_
  }
  
  return(s)
}

propresponsible=function(R0,k,prop){
  qm1=qnbinom(1-prop,k+1,mu=R0*(k+1)/k)
  remq=1-prop-pnbinom(qm1-1,k+1,mu=R0*(k+1)/k)
  remx=remq/dnbinom(qm1,k+1,mu=R0*(k+1)/k)
  q=qm1+1
  1-pnbinom(q-1,k,mu=R0)-dnbinom(q,k,mu=R0)*remx
}

sample_filter <- function(condition,df,col,n){
  sample(df %>% filter.(period==condition) %>% pull.(col),size=n,replace=T)
}

mean_filter <- function(condition,df,col){
  mean(df %>% filter.(period==condition) %>% pull.(col))
}

#### Main Model ----
run_model <- function(scenarios, browsing=F){
  
  if(browsing){browser()}
  
  #### Generate infections of hh (household) contacts ####
  indiv_params <- traj %>%  
    select.(-m) %>%
    crossing.(heterogen_contacts = unique(scenarios$heterogen_contacts),
              period=unique(scenarios %>% 
                              mutate.(period=fct_drop(period)) %>% 
                              pull.(period))) %>%   
    mutate.(hh_contacts=ifelse(heterogen_contacts,
                               sample_filter(condition = period, 
                                             df = contact_data_adjusted, 
                                             col="e_home", 
                                             n=n()),
                               round(mean_filter(condition = period, 
                                                 df = contact_data_adjusted, 
                                                 col="e_home"))),
            .by=c(period)) 
  
  indiv_params_long <- indiv_params %>% 
    left_join.(traj_)
  
  #simulate infections (and keep first instance)
  hh_infections <- indiv_params_long %>% 
    uncount.(hh_contacts,.id="id",.remove = F) %>% 
    mutate.(hh_duration = case_when.(heterogen_contacts ~ sample(contacts_hh_duration,
                                                                 size=n(),replace=T),
                                     TRUE               ~ median(contacts_hh_duration)),
            infected    = rbernoulli(n(),p=culture_p*hh_duration)) %>% 
    filter.(infected==T) %>% 
    slice.(min(t), .by=c(all_of(key_grouping_var),hh_contacts,id)) %>% 
    count.(t,all_of(key_grouping_var),hh_contacts,name = "hh_infected") %>% 
    arrange.(sim)
  
  #### Calculate nhh infections ####
  
  nhh_infections <- indiv_params_long %>% 
    
    # Sample daily contacts
    mutate.(nhh_contacts = ifelse(heterogen_contacts,
                                  sample_filter(condition = period,
                                                df=contact_data_adjusted,
                                                col="e_other",n=n()),
                                  round(mean_filter(condition = period,
                                                    df=contact_data_adjusted,
                                                    col="e_other"))),
            .by=c(period)) %>% 
    
    # Simulate infections 
    uncount.(nhh_contacts,.remove = F) %>% 
    mutate.(nhh_duration = case_when.(heterogen_contacts ~ sample(contacts_nhh_duration,
                                                                  size=n(),replace=T),
                                      TRUE               ~ median(contacts_nhh_duration)),
            nhh_infected = rbernoulli(n=n(),p = culture_p*nhh_duration)) %>% 
    summarise.(nhh_infected=sum(nhh_infected),.by=c(t,all_of(key_grouping_var),nhh_contacts,test)) %>% 
    
    # Testing: determine if and when testing + isolating by specified sampling frequency, adherence  
    right_join.(testing_scenarios) %>% 
    mutate.(
      test_day = case_when.((t - begin_testing) %% sampling_freq == 0 ~ TRUE,
                            nhh_contacts > event_size ~ TRUE,
                            TRUE ~ FALSE)) %>% 
    mutate.(
      earliest_pos = min(t[test&test_day]),
      test_iso = t>=earliest_pos & self_iso_test,
      .by=c(all_of(key_grouping_var),prop_self_iso_test,self_iso_test,begin_testing,sampling_freq,event_size)) %>%
    filter.(test_iso==F) %>% 
    select.(everything(),-test_iso,-test,-earliest_pos,-test_day) 
  
  # Join nhh and hh contacts and summarise
  processed_infections <- indiv_params_long %>% 
    right_join.(testing_scenarios) %>% 
    left_join.(hh_infections) %>% 
    left_join.(nhh_infections) %>% 
    replace_na.(list(hh_infected=0,nhh_infected=0,nhh_contacts=0)) %>% 
    arrange.(period,lower_inf_thresh) %>% 
    mutate.(
      total_contacts = nhh_contacts+hh_contacts,
      total_infections=nhh_infected+hh_infected)
}

#function to calculate the proportion above or below a defined threshold
prop_n <- function(df, threshold=10, col=e_all, op=">="){
  browser()
  df %>% 
    summarise.(n=n(),
               s = sum(match.fun(op)({{col}},{{threshold}})),
               "prop_{{threshold}}":= s/n) %>% 
    rename.("n_{{threshold}}":= s)
}

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24,  0.38),
  p = c(0.025, 0.975), 
  show.output = F, plot = F) %>%
  as.list

approx_sd <- function(x1, x2){
  (x2-x1) / (qnorm(0.95) - qnorm(0.05) )
}

#bootstrap confidence interval function
boot_ci <- function(x,nrep=100) {
  
  trueval <- tibble(param=c("mu","size"),
                    mean=c(x$estimate[[2]],
                           x$estimate[[1]])) 
  
  ci <- bootdist(f=x,niter = nrep)$CI %>% 
    as.data.frame() %>% 
    select(-Median) %>% 
    rownames_to_column("param")
  
  left_join(trueval,ci)
}

#run code quietly
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

# logarithmic spaced sequence
# blatantly stolen from library("emdbook"), because need only this
lseq <- function(from=1, to=100000, length.out=6) {
  
  exp(seq(log(from), log(to), length.out = length.out))
}

#quantile function
quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}
