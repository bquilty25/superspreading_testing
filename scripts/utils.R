# Load required packages scripts
pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","dtplyr","rms","DescTools","MESS","lubridate","lemon","boot","furrr","dtplyr","data.table","tidytable","ggtext","fst","extraDistr","emdbook","colorspace")

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

vl_params <- bind_rows(kissler_dat_est,hay_dat_est)

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
  out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
  return(out) 
}

#Load contact data
contacts_bbc_o18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_o18.csv")) 

contacts_bbc_u18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_u18.csv")) 

contacts_bbc <- bind_rows(contacts_bbc_o18, contacts_bbc_u18) %>%
  mutate(time_period = "pre",
         e_school=0)

contacts_comix_o18 <-
  read.csv(here("data", "comix_contact_distributions_o18.csv")) %>%
  mutate(
    date = as.Date(date),
    age="adults"
  )

contacts_comix_u18 <-
  read.csv(here("data", "comix_contact_distributions_u18.csv")) %>%
  mutate(
    date = as.Date(date),
    age="children"
  )

contacts_comix <- bind_rows(contacts_comix_o18, contacts_comix_u18)

comix_high_low <- contacts_comix %>%
  mutate(year = year(as.Date(date)),
         month = month(as.Date(date))) %>%
  filter(year == 2020 & month %in% c(8, 9) |
           year == 2021 & month %in% c(1, 2)) %>%
  mutate(
    time_period = case_when(
      year == 2020 & month %in% c(8, 9) ~ "aug_sept",
      year == 2021 & month %in% c(1, 2) ~ "jan_feb"
    )
  )

#summarise number of contacts
contact_data <- comix_high_low %>% 
  bind_rows(contacts_bbc)%>% 
  mutate(time_period=fct_relevel(time_period,"pre","aug_sept","jan_feb"),
         e_work_school=rowSums(across(c(e_work,e_school),na.rm = T)),
         e_all = rowSums(across(c(e_home,e_work_school,e_other)),na.rm = T)) %>% 
  select(-c(e_work,e_school))

##### KCL ANALYSIS ----
pickering <- readxl::read_xlsx(here("data","pickering_dat.xlsx")) %>% 
  select(-c(`Viral Growth`,...7,...8)) %>% 
  rename("culture"=...6) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)ifelse(x=="ND",NA,x)) %>% 
  mutate_at(.vars=vars(`SureScreen F`,Innova,Encode),
            .funs = function(x)case_when(x%in%c(0.5,1,2)~1,
                                         is.na(x)~NA_real_,
                                         TRUE~0)) %>% 
  mutate(id=row_number()) %>% 
  rename(ct=`Ct N1`)

innova_mod <- glm(Innova~ct,data=pickering,family="binomial") 
culture_mod <- glm(culture~ct,data=pickering,family="binomial") 

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
    n_cases = 100, n_sims = 100,
    asymp_parms = asymp_fraction,
    variant_info, browsing = FALSE
){
  
  if (browsing) browser()
  
  set.seed(seed)
  #simulate CT trajectories
  
  #inf <- tidytable(sim=1:n_sims)
  inf <- rbbinom(n = n_sims,
                 size=1,
                 alpha = asymp_parms$shape1,
                 beta = asymp_parms$shape2) %>%
    as_tidytable() %>%
    rename.("asymptomatic"=x) %>%
    mutate.(sim=row_number.(),
            asymptomatic=as.logical(asymptomatic))

  traj <- inf %>% 
    crossing.(start=0) %>% 
    crossing.(vl_params) %>% 
    mutate.(prolif=round(rnormTrunc(n=n(),mean=mean_prolif,sd=sd_prolif,min = 0.25,max=14)),
            clear=round(rnormTrunc(n=n(), mean=mean_clear, sd=sd_clear, min = 2,max=30)),
            # prolif=ifelse(asymptomatic,prolif*0.8,prolif),
            # clear=ifelse(asymptomatic,clear*0.8,clear),
            end=prolif+clear,
            onset_t=prolif+rnormTrunc(n=n(),mean = 2,sd=1.5,min=0,max=end)
    ) %>%
    select.(-c(mean_prolif, sd_prolif, mean_clear, sd_clear,clear)) %>%
    pivot_longer.(cols = -c(sim,variant,#,onset_t
                            mean_peakvl,sd_peakvl),
                  values_to = "x") %>%
    mutate.(y=case_when(name=="start" ~ 40,#convert_Ct_logGEML(40),
                        name=="end"   ~ 40,#convert_Ct_logGEML(40),
                        name=="prolif"~rnormTrunc(n=n(),mean=mean_peakvl,sd=sd_peakvl,min=0,max=40))) %>% 
    select.(-c(mean_peakvl,sd_peakvl))
  
  models <- traj %>%
    nest.(data = -c(variant,onset_t)) %>%  
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
    select.(c(sim, variant, onset_t, prolif, start, end, m)) %>% 
    arrange(sim)
  
}

inf_curve_func <- function(m,start=0,end=30,trunc_t){
  #browser()
  x <- tidytable(t=seq(start,end,by=1)) %>% 
    mutate.(vl=m(t))
  
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