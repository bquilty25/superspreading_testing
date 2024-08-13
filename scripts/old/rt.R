pacman::p_load(EpiNow2,covidregionaldata)

uk_dat <- get_national_data(countries = "UK")

uk_dat_2020 <- uk_dat %>% filter(date<"2021-01-01") %>% select(date,cases_new) %>% rename("confirm"=cases_new)

generation_time <- get_generation_time(disease = "SARS-CoV-2", source = "ganyani")
incubation_period <- get_incubation_period(disease = "SARS-CoV-2", source = "lauer")

estimates <- epinow(reported_cases = uk_dat_2020,
                    generation_time = generation_time,
                    rt = rt_opts(prior = list(mean = 2, sd = 0.2)),
                    stan = stan_opts(cores = 4))
