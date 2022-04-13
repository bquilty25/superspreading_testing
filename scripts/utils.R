# Load required packages scripts
pacman::p_load(
  "fitdistrplus",
  "EnvStats",
  "tidyverse",
  "patchwork",
  "here",
  "rriskDistributions",
  "dtplyr",
  "rms",
  "DescTools",
  "MESS",
  "lubridate",
  "lemon",
  "boot",
  "furrr",
  "dtplyr",
  "data.table",
  "tidytable",
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
  "RGeode"
)

seed <- 1000

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

covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

`%!in%` = Negate(`%in%`)

plotting_theme <- theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
        panel.grid = element_blank(),
        legend.position = "bottom",
        strip.placement = "outside")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

time_periods <- tribble(~idx,~period,~date_start,~date_end,
                        -1, "POLYMOD",             as_date("01/01/2008",format="%d/%m/%Y"), as_date("01/01/2008",format="%d/%m/%Y"),
                        0, "BBC Pandemic",         as_date("01/09/2017",format="%d/%m/%Y"), as_date("01/12/2018",format="%d/%m/%Y"),
                        1, "Lockdown 1",           as_date("23/03/2020",format="%d/%m/%Y"), as_date("03/06/2020",format="%d/%m/%Y"),
                        2, "Lockdown 1 easing",    as_date("04/06/2020",format="%d/%m/%Y"), as_date("29/07/2020",format="%d/%m/%Y"),
                        3, "Relaxed restrictions", as_date("30/07/2020",format="%d/%m/%Y"), as_date("03/09/2020",format="%d/%m/%Y"),
                        4, "School reopening",     as_date("04/09/2020",format="%d/%m/%Y"), as_date("24/10/2020",format="%d/%m/%Y"),
                        5, "Lockdown 2",           as_date("05/11/2020",format="%d/%m/%Y"), as_date("02/12/2020",format="%d/%m/%Y"),
                        6, "Lockdown 2 easing",    as_date("03/12/2020",format="%d/%m/%Y"), as_date("19/12/2020",format="%d/%m/%Y"),
                        7, "Lockdown 3",           as_date("05/01/2021",format="%d/%m/%Y"), as_date("07/03/2021",format="%d/%m/%Y"),
                        8, "Lockdown 3 + schools", as_date("08/03/2021",format="%d/%m/%Y"), as_date("31/03/2021",format="%d/%m/%Y"),
                        9, "Step 2 + schools",     as_date("16/04/2021",format="%d/%m/%Y"), as_date("16/05/2021",format="%d/%m/%Y")) %>% 
  mutate(period=factor(period,levels = period))

#Load contact data
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

contacts_comix_o18 <-
  read.csv(here("data", "comix_contact_distributions_o18.csv")) %>%
  mutate(
    date = as_date(date),
    age="adults"
  )

contacts_comix_u18 <-
  read.csv(here("data", "comix_contact_distributions_u18.csv")) %>%
  mutate(
    date = as_date(date),
    age="children"
  )

contacts_comix <- bind_rows(contacts_comix_o18, contacts_comix_u18)

#summarise number of contacts
contact_data <- contacts_comix %>% 
  bind_rows(contacts_bbc)%>% 
  mutate(e_other=rowSums(across(c(e_work,e_school,e_other)),na.rm = T)) %>% 
  select(date,e_home,e_other) %>% 
  bind_rows(contacts_polymod) %>% 
  mutate(e_all = rowSums(across(c(e_home,e_other)),na.rm = T),
         date=as_date(date))%>% 
  fuzzyjoin::fuzzy_inner_join(time_periods,
                   by=c("date"="date_start","date"="date_end"),
                   match_fun=list(`>=`,`<=`))

# impute out of HH values > 250 for BBC pandemic by fitting distribution to values from non-lockdown periods
contact_data %>%summarise.(n=n(),over_250=sum(e_other>=250),.by=period) %>% mutate.(prop=over_250/n*100,.by=period)

dist_over_250 <- contact_data %>% 
  filter(period%in%c("Relaxed restrictions","School reopening","Step 2 + schools")) %>% 
  filter(e_other>=250)  %>% 
  pull(e_other) %>% 
  fitdistr(.,"exponential")

data.frame(x=seq(250,4000,by=1)) %>% 
  mutate(y=dexp(x,rate=dist_over_250$estimate[1])) %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  geom_histogram(data=contact_data %>% 
                   filter(period%in%c("Relaxed restrictions","School reopening")) %>% 
                   filter(e_other>=250),
                 aes(x=e_other,y=stat(density)),alpha=0.5)

data.frame(x=rexptr(n=0.002*nrow(contact_data %>% filter(period=="BBC Pandemic")),
                  lambda=dist_over_250$estimate[1],range = c(250,Inf))) %>% 
  ggplot()+
  geom_histogram(aes(x=x,stat(density)))+
  geom_point(data= data.frame(x=seq(250,4000,by=1)) %>% mutate(y=dexp(x,rate=dist_over_250$estimate[1])),aes(x=x,y=y))

dat_append <- data.frame(e_other=round(rexptr(n=0.0025*nrow(contact_data %>% filter(period=="BBC Pandemic")),
                                        lambda=dist_over_250$estimate[1],range = c(250,Inf))),
                         e_home=sample(size=0.0025*nrow(contact_data %>% filter(period=="BBC Pandemic")),
                                      x=contact_data %>% filter(period=="BBC Pandemic") %>% pull(e_home))) %>% 
  mutate(e_all=e_home+e_other,period="BBC Pandemic",idx=3)

contact_data <- contact_data %>% bind_rows(dat_append)

prop_n <- function(df, threshold=10, col=e_all, op=">="){
  df %>% 
    summarise(n=n(),
              s = sum(match.fun(op)({{col}},{{threshold}})),
              "prop_{{threshold}}":= s/n) %>% 
    rename("n_{{threshold}}":= s)
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

test_model_choice <- function(boolean){
  if(boolean){
    innova_higher_mod
  }else{
    innova_mod
  }
}

culture_mod <- glm(culture~vl,data=pickering,family="binomial") 

inf_model_choice <- function(boolean){
  #browser()
  if(boolean){
    innova_mod
  }else{
    culture_mod
  }
}

make_trajectories <- function(
    n_sims = 100,
    asymp_parms = asymp_fraction,
    variant_info, 
    browsing = FALSE
){
  
  if (browsing) browser()
  
  set.seed(seed)
  #simulate CT trajectories
  
  #inf <- tidytable(sim=1:n_sims)
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
    mutate.(prolif=case_when.(heterogen_vl~rnormTrunc(n=n(),mean=mean_prolif,sd=sd_prolif,min = 1,max=14),
                              TRUE~mean_prolif),
            clear=case_when.(heterogen_vl~rnormTrunc(n=n(), mean=mean_clear, sd=sd_clear, min = 1,max=30),
                             TRUE~mean_clear),
            # prolif=ifelse(asymptomatic,prolif*0.8,prolif),
            # clear=ifelse(asymptomatic,clear*0.8,clear),
            end=prolif+clear,
            onset_t=prolif+rnormTrunc(n=n(),mean = 2,sd=1.5,min=0,max=end)
    ) %>%
    select.(-c(mean_prolif, sd_prolif, mean_clear, sd_clear,clear)) %>%
    pivot_longer.(cols = -c(sim,variant,onset_t, asymptomatic, heterogen_vl,
                            mean_peakvl,sd_peakvl),
                  values_to = "x") %>%
    mutate.(y=case_when(name=="start" ~ 40,#convert_Ct_logGEML(40),
                        name=="end"   ~ 40,#convert_Ct_logGEML(40),
                        name=="prolif"~case_when.(heterogen_vl~rnormTrunc(n=n(),mean=mean_peakvl,sd=sd_peakvl,min=0,max=40),
                                                  TRUE~mean_peakvl))) %>% 
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

test_times <- function(type,onset_t,sampling_freq=3){
  #browser()
  
   if(asymptomatic){
     initial_t <- sample(size=1,x = c(0:29))
  }else{
    initial_t <- sample(size=1,x = c(0:onset_t))
  }
  
  if(!is.na(sampling_freq)){
  test_timings <- data.frame(test_t = seq(from=initial_t,to=30,by=sampling_freq)) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  } else {
    test_timings <- data.frame(test_t = Inf) %>% 
    mutate(test_no = paste0("test_", row_number())) 
  }

  
  return(test_timings)
}

earliest_pos <- function(df){
  #browser()
  
  x_q <- df[(test_label)]
  
  if (nrow(x_q) == 0L){
    return(tidytable(test_no="None",test_p=0,test_t=Inf))
  } else {
    return(x_q %>% select.(test_no,test_p,test_t) %>% slice_min.(test_t))
  }
}

detector <- function(test_p, u = NULL){
  
  if (is.null(u)){
    u <- runif(n = length(test_p))
  }
  
  # true positive if the test exceeds a random uniform
  # when uninfected, PCR will be 0
  TP <- test_p > u
  
}

inf_and_test <- function(traj,sampling_freq=c(NA,3)){
  #browser()
  
  message(sprintf("\n%s == SCENARIO %d ======", Sys.time(), traj$sim[1]))
  
  traj %>% as.data.frame() %>% 
    mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m,start=start,end=end)))  %>% 
    unnest_wider(infectiousness) %>% 
    ungroup() %>%
    mutate.(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf))) %>% 
    #testing
    crossing(sampling_freq = sampling_freq) %>% 
    mutate.(test_times = pmap(
      .f = test_times,
      list(
        sampling_freq = sampling_freq,
        onset_t = onset_t,
        type = type
      )
    )) %>%
    unnest.(test_times,.drop=F) %>%
    mutate.(
      ct = pmap_dbl(.f = calc_sensitivity, list(model = m, x = test_t)),
      test_p = stats::predict(innova_mod, type = "response", newdata = data.frame(ct = ct)),
      test_label = detector(test_p = test_p,  u = runif(n = n(), 0, 1))
    ) %>%
    nest(ct, test_t, test_no, test_p, test_label) %>%
    mutate.(earliest_positive = map(.f = earliest_pos, .x = data)) %>%
    unnest.(earliest_positive,.drop=F) %>%
    select.(-data)
} 

sample_contacts <- function(time_period){
  sample(contact_data %>% filter(period==time_period) %>% pull(e_home),size=1)
}

sec_case_gen <- function(df){
  
  message(sprintf("\n%s", Sys.time()))

  df1 <- df %>%
    mutate.(self_iso_symp=ifelse(type=="symptomatic",rbinom(n=n(),size=1,prob=prop_self_iso_symp),0),
            self_iso_test=rbinom(n=n(),size=1,prob=prop_self_iso_test),
            test_t = ifelse(self_iso_test==0,Inf,test_t)) %>% 
    #select.(-u) %>%
    mutate.(.by=time_period,
      contacts_repeated = sample(contact_data$e_home[contact_data$period==time_period],size=n(),replace=T),
      trunc_t=case_when.(
        # if symptomatic, adhering to self isolation, and either not tested or test neg,
        # truncate at onset
        type == "symptomatic" & is.infinite(test_t) & self_iso_symp != 0 ~ onset_t,
        # if symptomatic, adhering to self isolation, and have onset before test, truncate at onset
        type == "symptomatic" &
          is.finite(test_t) & onset_t < test_t & self_iso_symp != 0 ~ onset_t,
        # if symptomatic, adhering to self isolation, and have onset after pos test, truncate at test
        type == "symptomatic" &
          is.finite(test_t) & test_t < onset_t & self_iso_symp != 0 ~ test_t,
        # if symptomatic, not adhering to self isolation, and have a positive test, truncate at test
        TRUE ~ test_t)) %>% 
    mutate.(data=pmap(.l=list(x=infectiousness,
                              contacts_repeated=contacts_repeated),.f=rep_contacts_inf)) %>% 
    unnest.(data) %>% 
    mutate.(.by=time_period,
      contacts_casual = sample(contact_data$e_other[contact_data$period==time_period],size=n(),replace=T)) %>%
    mutate.(contacts_casual=ifelse(t>=trunc_t,0L, contacts_casual)) %>% 
    uncount.(contacts_casual) %>% 
    mutate.(nhh_duration=sample(contacts_nhh_duration,size=n(),replace=T),
      n_casual_infected = rbernoulli(n=n(),p = culture*nhh_duration)) %>%
    summarise.(.by=c(sim,
                       idx,
                       t,
                       ct,
                       onset_t,
                       type,
                       variant,
                       time_period,
                       prop_self_iso_symp,
                       prop_self_iso_test,
                       contacts_repeated,
                       n_repeated_infected,
                       sampling_freq),
                       contacts_casual=n(),
                       n_casual_infected=sum(n_casual_infected))

}

hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}
