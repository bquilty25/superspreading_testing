
source("scripts/utils.R")
source("scripts/lee_infectivity.R")

res <- fst::read.fst(res %>% bind_rows.(),path = "results_ss.fst")

inf_curve1 <- make_trajectories(n_cases = 500, n_sims = 200) %>% 
  as_tibble() %>% 
  mutate(infectiousness = pmap(inf_curve_func, .l = list(m = m)))  %>% 
  unnest_wider(infectiousness) %>% 
  ungroup() %>%
  mutate.(norm_sum = (sum_inf - min(sum_inf)) / (max(sum_inf) - min(sum_inf))) %>% 
  mutate.(test_t=Inf,
         test_no="None",
         test_p=0,
         prop_self_iso_symp = c(0),
         prop_self_iso_test = c(1)) %>% 
  crossing(het_contacts=c(TRUE,FALSE),
           het_vl=c(TRUE,FALSE)) 

sec_case_hetero <- function(df){
  
  message(sprintf("\n%s", Sys.time()))
  
 #browser()
  df %>%
    mutate.(self_iso_symp=ifelse(type=="symptomatic",rbinom(n=n(),size=1,prob=prop_self_iso_symp),0),
            self_iso_test=rbinom(n=n(),size=1,prob=prop_self_iso_test),
            test_t = ifelse(self_iso_test==0,Inf,test_t)) %>% 
    select.(-u) %>%
    mutate.(
      contacts_hh_home = case_when.(
        het_contacts ~ sample(
          contact_data %>%
            filter(time_period == "pre") %>%
            pull(e_home),
          size = n(),
          replace = T
        ),
        !het_contacts ~ contact_data %>%
          filter(time_period == "pre") %>%
          summarise(e_home=ceiling(mean(e_home,na.rm=T))) %>% 
          pull(e_home)
        ),
      contacts_hh_work_school = case_when(
        het_contacts ~ sample(
          contact_data %>%
            filter(time_period == "pre") %>%
            pull(e_work_school),
          size = n(),
          replace = T
        ),
        !het_contacts ~ contact_data %>%
          filter(time_period == "pre") %>%
          summarise(e_work_school=ceiling(mean(e_work_school,na.rm=T))) %>% 
          pull(e_work_school)),
      trunc_t=case_when.(
        # if symptomatic, adhering to self isolation, and either not tested or test neg,
        # truncate at onset
        type=="symptomatic"&is.infinite(test_t)&self_iso_symp!=0~onset_t,
        # if symptomatic, adhering to self isolation, and have onset before test, truncate at onset
        type=="symptomatic"&is.finite(test_t)&onset_t<test_t&self_iso_symp!=0~onset_t,
        # if symptomatic, adhering to self isolation, and have onset after pos test, truncate at test
        type=="symptomatic"&is.finite(test_t)&test_t<onset_t&self_iso_symp!=0~test_t,
        # if symptomatic, not adhering to self isolation, and have a positive test, truncate at test
        TRUE ~ test_t), 
      #if symp onset or test positive, contacts outside of home should cease; for school/work, multiply 
      #number of contacts by proportion of time since infection which occurs pre-onset
      contacts_hh_work_school=case_when.(!is.infinite(trunc_t)~ ceiling(contacts_hh_work_school*(trunc_t/(end-start))),
                                               TRUE                 ~ ceiling(contacts_hh_work_school))
    ) %>% 
    mutate.(
      contacts_hh   = contacts_hh_home + contacts_hh_work_school,
      n_hh_infected = rbinom(n = n(), 
                                   size = contacts_hh, 
                                   prob = case_when.(het_vl~norm_sum,
                                                     !het_vl~mean(norm_sum,na.rm=T)))) %>%
    #Daily nhh contacts
    unnest.(infectiousness) %>%
    mutate.(norm_daily = case_when.(het_vl ~ culture / sum_inf * norm_sum,
                                    !het_vl ~ mean(culture / sum_inf * norm_sum,na.rm=T))) %>% 
    mutate.(
      contacts_nhh = case_when.(
        t>trunc_t  ~ 0L,
        het_contacts ~ sample(
          contact_data %>%
            filter(time_period == "pre") %>%
            pull(e_other),
          size = n(),
          replace = T),
        !het_contacts ~ contact_data %>%
          filter(time_period == "pre") %>%
          summarise(e_other=ceiling(mean(e_other,na.rm=T))) %>% 
          pull(e_other))
    ) %>%
    mutate.(n_nhh_infected = rbinom(n = n(), 
                                       size = contacts_nhh, 
                                       prob =
                                         case_when.(het_vl~norm_daily,
                                                    !het_vl~mean(norm_daily,na.rm=T)))) %>% 
    summarise.(.by=c(sim,
                     idx,
                     type,
                     norm_sum,
                     het_contacts,
                     het_vl,
                     prop_self_iso_symp,
                     prop_self_iso_test,
                     contacts_hh,
                     n_hh_infected),
               contacts_nhh = sum(contacts_nhh),
               n_nhh_infected = sum(n_nhh_infected)) %>%
    mutate.(n_total_infected = n_hh_infected + n_nhh_infected) 
}

res1 <- inf_curve1 %>% 
  group_split.(het_contacts,het_vl) %>% 
  map(~sec_case_hetero(.x))

dists <- res1 %>% 
  bind_rows.() %>%  
  nest_by.(sim,het_contacts,het_vl) %>% 
  mutate.(dists=map(.x=data,.f= .%>% 
                      pull(n_total_infected) %>% 
                      fitdist("nbinom")),
          params=map(dists,~c(R=.x$estimate[[2]],k=.x$estimate[[1]]))) %>% 
  unnest_longer(params,values_to = "value",indices_to = "param") %>% 
  nest_by.(het_contacts,
           het_vl,
           param) %>% 
  mutate.(q = map(.x = data, ~quantile(.x$value,probs=c(0.5,0.025,0.975)))) %>% 
 unnest_wider(q) 

ggplot(res1 %>% bind_rows())+geom_histogram(aes(x=n_total_infected))+facet_grid(het_vl~het_contacts,labeller = label_both)#+scale_y_log10()

res1 %>% bind_rows() %>% 
  group_by(het_contacts,het_vl) %>% 
  count(n_total_infected>10) %>% 
  mutate(prop=prop.table(n)) %>% 
  as.tibble() %>% base::print(n=Inf)

dists %>% select(het_contacts,het_vl,param,`50%`) %>% pivot_wider(,names_from=param,values_from = `50%`)
