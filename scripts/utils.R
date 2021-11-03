# Load required packages scripts
pacman::p_load("fitdistrplus","EnvStats","tidyverse","patchwork","here","rriskDistributions","DescTools","MESS","lubridate","lemon","boot","furrr","data.table","tidytable","ggtext","fst","qs","tictoc","scales","fuzzyjoin")

set.seed(2021)
plan(multisession,workers=8)

covid_pal <- c("#e66101", "#5e3c99", "#0571b0")
`%!in%` = Negate(`%in%`)

plotting_theme <- theme_minimal()+
  theme(axis.ticks = element_line(),
        panel.border = element_rect(fill=NA),
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
                        9, "Step 2 + schools",     as_date("16/04/2021",format="%d/%m/%Y"), as_date("16/05/2021",format="%d/%m/%Y"))

#Load contact data
contacts_polymod <- 
  read.csv(here("data","2008_Mossong_POLYMOD_contact_common.csv")) %>% 
  pivot_longer.(cols=c(cnt_home,cnt_school,cnt_work,cnt_transport,cnt_leisure,cnt_otherplace)) %>% 
  filter(value) %>% 
  select.(-value) %>% 
  count.(name,part_id) %>% 
  pivot_wider(names_from = name,values_from = N) %>% 
  mutate(e_home=cnt_home,
    e_other=rowSums(across(c(cnt_work,cnt_school,cnt_transport,cnt_leisure,cnt_otherplace)),na.rm = T)) %>% 
  select.(part_id,e_home,e_other) %>% 
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
  fuzzy_inner_join(time_periods,
                   by=c("date"="date_start","date"="date_end"),
                   match_fun=list(`>=`,`<=`))

prop_n <- function(df, threshold=10, col=e_all, op=">="){
  df %>% 
    summarise(n=n(),
              s = sum(match.fun(op)({{col}},{{threshold}})),
              "prop_{{threshold}}":= s/n) %>% 
    rename("n_{{threshold}}":= s)
}

# # McAloon et al. incubation period meta-analysis
#https://bmjopen.bmj.com/content/10/8/e039652
inc_parms <- list(mu_inc = 1.63,
                  sigma_inc = 0.5)

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

#https://science.sciencemag.org/content/sci/early/2021/05/24/science.abi5273.full.pdf
#assuming 1 and 3 days are 95% interval:
peak_to_onset <- rriskDistributions::get.norm.par(p=c(0.025,0.975),q=c(1,3),plot = F)

make_trajectories <- function(n_cases=100, n_sims=100, seed=1000,asymp_parms=asymp_fraction,variant=c("delta")){
  
  set.seed(seed)
  #simulate CT trajectories
  #browser()

  inf <- data.frame(sim=1:n_sims) %>% 
    mutate.(prop_asy = rbeta(n = n(),
                            shape1 = asymp_parms$shape1,
                            shape2 = asymp_parms$shape2)) 
  
  inf %<>%  
    mutate.(x = map.(.x = prop_asy,
                   .f = ~make_proportions(prop_asy     = .x,
                                          n_cases      = n_cases))) %>% 
    unnest.(x) 
  
  traj <- inf %>% 
    crossing.(start=0) %>% 
    crossing.(variant=variant) %>% 
    # duration from: https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(20)30172-5/fulltext
    # scaling of asymptomatics taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    # mutate(end=case_when(type == "symptomatic"  ~ qnormTrunc(p = u, mean=17, 
    #                                                          sd=approx_sd(15.5,18.6), min = 0),
    #                      type == "asymptomatic" ~ qnormTrunc(p = u, mean=17*0.6, 
    #                                                          sd=approx_sd(15.5,18.6), min = 0))) %>% 
    # incubation period from https://bmjopen.bmj.com/content/10/8/e039652.info
    mutate.(peak=rnormTrunc(n=n(),mean=3.2,sd=approx_sd(2.4, 4.2),min=0),
           clear=case_when.(type == "symptomatic"  ~ rnormTrunc(n=n(), mean=10.9, 
                                                               sd=approx_sd(7.8, 14.2), min = 0),
                           type == "asymptomatic" ~ rnormTrunc(n=n(), mean=7.8, 
                                                               sd=approx_sd(6.1, 9.7), min = 0)),
           clear=case_when.(variant=="delta"~1.4*clear,
                           #variant=="delta"&name=="end"~0.71*y,
                           TRUE~clear),
           end=peak+clear,
           onset_t=peak+rnormTrunc(n=n(),mean = 2,sd=0.5,min=0,max=end)
           
           #onset_t=qlnormTrunc(p = u,
           #                            meanlog=1.63,
           #                            sdlog=0.5,
           #                            min = 0,
           #                            max=end)
    ) %>% 
    select.(-clear) %>% 
    pivot_longer.(cols = -c(sim,prop_asy,idx,type,variant,onset_t),
                 values_to = "x") %>%
    # peak CT taken from https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v2
    mutate.(y=case_when(name=="start" ~ 40,
                       name=="end"  ~ 40,
                       name=="peak"&variant!="delta" ~ rnorm(n=n(),mean=22.4,sd=approx_sd(20.7, 24)),
                       name=="peak"&variant=="delta" ~ 0.7*rnorm(n=n(),mean=22.4,sd=approx_sd(20.7, 24))
    ))  #multiply Ct by 0.7 for Delta variant (24Ct /34Ct): https://virological.org/t/viral-infection-and-transmission-in-a-large-well-traced-outbreak-caused-by-the-delta-sars-cov-2-variant/724
  
  #browser()
  
  models <- traj %>%
    nest.(data = -c(idx,type,variant,onset_t)) %>%  
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
    left_join.(x_model)

}

inf_curve_func <- function(m,start=0,end=30,trunc_t){
  #browser()
  x <- data.frame(t=seq(start,end,by=1)) %>% 
    mutate(u=runif(n=n(),0,1))
  
  #predict CTs for each individual per day
  x$ct <- m(x$t)
  
  #predict culture probability given CTs
  x$culture <-  stats::predict(culture_mod, type = "response", newdata = x)

  sum_inf <- sum(x$culture)
  
  return(list(infectiousness=x,sum_inf=sum_inf))
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

## sample asymp proportions
make_proportions <- function(prop_asy, n_cases){
  
  props <- c("symptomatic"  = (1 - prop_asy),
             "asymptomatic" = prop_asy)
  
  x <- data.frame(type=rbinom(n=n_cases,size=1,prob = prop_asy)) %>% 
    mutate(type=ifelse(type==1,"asymptomatic","symptomatic"),
           idx=row_number())
  
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
  
   if(type=="asymptomatic"){
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

rep_contacts_inf <- function(x,contacts_repeated){
  #browser()
  
  #repeated contacts
  # initial values at t=0

  infected <- c()
  infected[1] <- 0

  #chance of infecting repeated contacts on a given day is binomial w/o replacement
  
  for(i in 2:nrow(x)){
    infected[i] <- infected[i-1] + rbinom(1,prob=x$culture[i],size=contacts_repeated-infected[i-1])
  }
  
  x$n_repeated_infected <- c(infected[1],diff(infected))
  
 return(x)
}

hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}