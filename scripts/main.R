# Estimate secondary case distribution pre-pandemic (R0, BBC) and with various levels of contact reduction from CoMix
source("scripts/utils.R")
source("scripts/duration.R")

N_sims <- 10000
# Make VL trajectories
traj <- vl_params %>%
  filter.(variant %in% c("wild")) %>%
  mutate.(variant = fct_drop(variant)) %>%
  crossing(heterogen_vl = c(TRUE, FALSE)) %>%
  group_split.(variant, heterogen_vl) %>%
  map.(~ make_trajectories(n_sims = N_sims, asymp_parms = asymp_fraction, variant_info = .x, browsing = F)) %>%
  bind_rows.()

infctsnss_params <- generate_params(culture_mod, N_sims) %>%
  as_tidytable(.) %>%
  rename.(beta0 = "(Intercept)", beta1 = vl) %>%
  mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

# Calculate daily infectiousness and test positivity, remove never-infectious
traj_ <- traj %>%
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end, interval = 1
  ))) %>%
  unnest.(infectiousness) %>%
  crossing.(
    lower_inf_thresh = c(FALSE)
  ) %>%
  mutate.(
    culture_p = culture_prob(vl, beta0, beta1),
    infectious = rbernoulli(
      n = n(),
      p = 1 - exp(-beta_inf * culture_p * median_contact_duration)
    ),
    test_p = stats::predict(
      object = innova_mod,
      type = "response",
      newdata = tidytable(vl = vl)
    ),
    test = rbernoulli(
      n = n(),
      p = test_p
    ),
    .by = c(lower_inf_thresh)
  ) %>%
  replace_na.(list(
    test = FALSE,
    infectious = FALSE
  )) %>%
  select.(-c(prolif, start, end))

# Calibrate beta_inf so pre-pandemic mean R matches target
message("Calibrating beta_inf (heterogen_vl=T, heterogen_contacts=T) ...")
beta_inf <- calibrate_beta(target_R0 = 2.5, n_calib = 10000)
message(sprintf("Calibrated beta_inf = %.4f", beta_inf))
qsave(beta_inf, "results/calibrated_beta.qs")

# Calibrate betas for heterogeneity-off combinations so that R0=2.5 in each case,
# isolating the effect of heterogeneity on k rather than on both k and R0.
message("Calibrating beta_inf (heterogen_vl=T, heterogen_contacts=F) ...")
beta_inf_vl_on_contacts_off <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = TRUE, heterogen_contacts_flag = FALSE)
message(sprintf("Calibrated beta_inf (vl=T, contacts=F) = %.4f", beta_inf_vl_on_contacts_off))

message("Calibrating beta_inf (heterogen_vl=F, heterogen_contacts=T) ...")
beta_inf_vl_off_contacts_on <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = FALSE, heterogen_contacts_flag = TRUE)
message(sprintf("Calibrated beta_inf (vl=F, contacts=T) = %.4f", beta_inf_vl_off_contacts_on))

message("Calibrating beta_inf (heterogen_vl=F, heterogen_contacts=F) ...")
beta_inf_both_off <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = FALSE, heterogen_contacts_flag = FALSE)
message(sprintf("Calibrated beta_inf (vl=F, contacts=F) = %.4f", beta_inf_both_off))

qsave(
  list(
    vl_on_contacts_on  = beta_inf,
    vl_on_contacts_off = beta_inf_vl_on_contacts_off,
    vl_off_contacts_on = beta_inf_vl_off_contacts_on,
    vl_off_contacts_off = beta_inf_both_off
  ),
  "results/calibrated_betas.qs"
)

# baseline
testing_scenarios <- traj %>%
  filter.(heterogen_vl == T) %>%
  select.(-m) %>%
  crossing.(
    prop_self_iso_test = c(0),
    sampling_freq = c(7),
    event_size = NA
  ) %>%
  mutate.(
    self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
    begin_testing = rdunif(n(), 0, sampling_freq)
  )

time_periods_of_interest <-
  crossing(time_periods) %>%
  filter(date_end < as.Date("2021-01-01"), period != "POLYMOD") %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end)) %>%
  crossing(heterogen_contacts = c(T))

processed_infections_baseline <- run_model(testing_scenarios = testing_scenarios, contact_dat = contact_data, scenarios = time_periods_of_interest, browsing = F)

rm(testing_scenarios)
rm(time_periods_of_interest)

print("baseline done")
qsave(processed_infections_baseline, "results/processed_infections_baseline.qs")
rm(processed_infections_baseline)
gc()

# heterogen onoff — beta recalibrated per combination so R0=2.5 in all cases

time_periods_base_heterogen <- crossing(time_periods) %>%
  filter(date_end < as.Date("2021-01-01"), period != "POLYMOD") %>%
  # filter(period%in%c("Pre-pandemic","1st lockdown","School reopening")) %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end))

heterogen_combos <- list(
  list(het_vl = TRUE,  het_contacts = TRUE,  beta = beta_inf),
  list(het_vl = TRUE,  het_contacts = FALSE, beta = beta_inf_vl_on_contacts_off),
  list(het_vl = FALSE, het_contacts = TRUE,  beta = beta_inf_vl_off_contacts_on),
  list(het_vl = FALSE, het_contacts = FALSE, beta = beta_inf_both_off)
)

beta_inf_saved <- beta_inf
processed_infections_heterogen_on_off <- map(heterogen_combos, function(combo) {
  beta_inf <<- combo$beta
  ts <- traj %>%
    filter.(heterogen_vl == combo$het_vl) %>%
    select.(-m) %>%
    crossing.(prop_self_iso_test = c(0), sampling_freq = NA, event_size = NA) %>%
    mutate.(
      self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
      begin_testing = rdunif(n(), 0, sampling_freq)
    )
  run_model(
    testing_scenarios = ts,
    contact_dat = contact_data,
    scenarios = time_periods_base_heterogen %>% crossing(heterogen_contacts = combo$het_contacts),
    browsing = F
  )
}) %>% bind_rows()
beta_inf <- beta_inf_saved
rm(beta_inf_saved, time_periods_base_heterogen, heterogen_combos)

print("heterogen on/off done")
qsave(processed_infections_heterogen_on_off, "results/processed_infections_heterogen_on_off.qs")
rm(processed_infections_heterogen_on_off)
gc()

# testing
testing_scenarios <- traj %>%
  filter.(heterogen_vl == T) %>%
  select.(-m) %>%
  crossing.(
    prop_self_iso_test = seq(0, 1, by = 0.1),
    sampling_freq = c(1, 3, 7),
    event_size = NA
  ) %>%
  mutate.(
    self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
    begin_testing = rdunif(n(), 0, sampling_freq)
  )

time_periods_of_interest <-
  crossing(time_periods) %>%
  filter(period %in% c("Pre-pandemic", "1st lockdown", "School reopening")) %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end)) %>%
  crossing(heterogen_contacts = c(T))

processed_infections_testing <- run_model(testing_scenarios = testing_scenarios, contact_dat = contact_data, scenarios = time_periods_of_interest, browsing = F)

rm(testing_scenarios)
rm(time_periods_of_interest)

print("testing done")
qsave(processed_infections_testing, "results/processed_infections_testing.qs")
rm(processed_infections_testing)
gc()

# event testing

testing_scenarios <- traj %>%
  filter.(heterogen_vl == T) %>%
  select.(-m) %>%
  crossing.(
    prop_self_iso_test = seq(0, 1, by = 0.1), # c(0,.25,.5,0.75,1),
    sampling_freq = c(NA),
    event_size = c(10, 20, 50, NA)
  ) %>%
  filter.(!(prop_self_iso_test > 0 & is.na(sampling_freq) & is.na(event_size))) %>%
  mutate.(
    self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
    begin_testing = rdunif(n(), 0, sampling_freq)
  )

time_periods_of_interest <-
  crossing(time_periods) %>%
  filter(period %in% c("Pre-pandemic", "1st lockdown", "School reopening")) %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end)) %>%
  crossing(heterogen_contacts = c(T))

processed_infections_events <- run_model(testing_scenarios = testing_scenarios, contact_dat = contact_data, scenarios = time_periods_of_interest, browsing = F)

rm(testing_scenarios)
rm(time_periods_of_interest)
# source("scripts/results.R")

print("event testing done")
qsave(processed_infections_events, "results/processed_infections_events.qs")
rm(processed_infections_events)
gc()

### Sensitivity analysis

# imputing upper tail of distribution for BBC Pandemic
# baseline
testing_scenarios <- traj %>%
  filter.(heterogen_vl == T) %>%
  select.(-m) %>%
  crossing.(
    prop_self_iso_test = c(0),
    sampling_freq = c(7),
    event_size = NA
  ) %>%
  mutate.(
    self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
    begin_testing = rdunif(n(), 0, sampling_freq)
  )

time_periods_of_interest <-
  crossing(time_periods) %>%
  filter(date_end < as.Date("2021-01-01"), period == "Pre-pandemic") %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end)) %>%
  crossing(heterogen_contacts = c(T))

processed_infections_sens <- run_model(
  testing_scenarios = testing_scenarios,
  contact_dat = contact_data_adjusted,
  scenarios = time_periods_of_interest, browsing = F
)

rm(testing_scenarios)
rm(time_periods_of_interest)

print("sens done")
qsave(processed_infections_sens, "results/processed_infections_sens.qs")
rm(processed_infections_sens)
gc()

### Sensitivity analysis: amplified VL heterogeneity ----
# Doubles SDs of peak VL, proliferation, and clearance to test whether
# inflating unmeasured between-person infectiousness heterogeneity can
# reverse the conclusion that contacts dominate overdispersion.

vl_params_amplified <- vl_params %>%
  mutate(
    sd_peakvl = sd_peakvl * 2,
    sd_prolif = sd_prolif * 2,
    sd_clear = sd_clear * 2
  )

traj_amplified <- vl_params_amplified %>%
  filter.(variant %in% c("wild")) %>%
  mutate.(variant = fct_drop(variant)) %>%
  crossing(heterogen_vl = c(TRUE, FALSE)) %>%
  group_split.(variant, heterogen_vl) %>%
  map.(~ make_trajectories(
    n_sims = N_sims, asymp_parms = asymp_fraction,
    variant_info = .x,
    max_prolif = 28, max_clear = 60, max_peakvl = 80,
    browsing = F
  )) %>%
  bind_rows.()

traj_amplified <- traj_amplified %>% left_join.(infctsnss_params, by = "sim")

traj_amplified_ <- traj_amplified %>%
  mutate.(infectiousness = pmap(inf_curve_func, .l = list(
    m = m, start = start, end = end, interval = 1
  ))) %>%
  unnest.(infectiousness) %>%
  crossing.(lower_inf_thresh = c(FALSE)) %>%
  mutate.(
    culture_p = culture_prob(vl, beta0, beta1),
    infectious = rbernoulli(n = n(), p = 1 - exp(-beta_inf * culture_p * median_contact_duration)),
    test_p = stats::predict(
      object = innova_mod, type = "response",
      newdata = tidytable(vl = vl)
    ),
    test = rbernoulli(n = n(), p = test_p),
    .by = c(lower_inf_thresh)
  ) %>%
  replace_na.(list(test = FALSE, infectious = FALSE)) %>%
  select.(-c(prolif, start, end))


# Calibrate betas for amplified VL analysis — one per heterogeneity combination
# so R0=2.5 in each case, isolating heterogeneity effects on k.
message("Calibrating beta_inf for amplified VL (het_vl=T, het_contacts=T) ...")
beta_inf_amp_vl_on_contacts_on <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = TRUE, heterogen_contacts_flag = TRUE, traj_data = traj_amplified)
message(sprintf("Calibrated beta_inf (amp, vl=T, contacts=T) = %.4f", beta_inf_amp_vl_on_contacts_on))

message("Calibrating beta_inf for amplified VL (het_vl=T, het_contacts=F) ...")
beta_inf_amp_vl_on_contacts_off <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = TRUE, heterogen_contacts_flag = FALSE, traj_data = traj_amplified)
message(sprintf("Calibrated beta_inf (amp, vl=T, contacts=F) = %.4f", beta_inf_amp_vl_on_contacts_off))

message("Calibrating beta_inf for amplified VL (het_vl=F, het_contacts=T) ...")
beta_inf_amp_vl_off_contacts_on <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = FALSE, heterogen_contacts_flag = TRUE, traj_data = traj_amplified)
message(sprintf("Calibrated beta_inf (amp, vl=F, contacts=T) = %.4f", beta_inf_amp_vl_off_contacts_on))

message("Calibrating beta_inf for amplified VL (het_vl=F, het_contacts=F) ...")
beta_inf_amp_vl_off_contacts_off <- calibrate_beta(target_R0 = 2.5, n_calib = 10000,
  heterogen_vl_flag = FALSE, heterogen_contacts_flag = FALSE, traj_data = traj_amplified)
message(sprintf("Calibrated beta_inf (amp, vl=F, contacts=F) = %.4f", beta_inf_amp_vl_off_contacts_off))

qsave(
  list(
    vl_on_contacts_on  = beta_inf_amp_vl_on_contacts_on,
    vl_on_contacts_off = beta_inf_amp_vl_on_contacts_off,
    vl_off_contacts_on = beta_inf_amp_vl_off_contacts_on,
    vl_off_contacts_off = beta_inf_amp_vl_off_contacts_off
  ),
  "results/calibrated_betas_amplified.qs"
)

time_periods_base_amp <- crossing(time_periods) %>%
  filter(date_end < as.Date("2021-01-01"), period != "POLYMOD") %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end))

heterogen_combos_amp <- list(
  list(het_vl = TRUE,  het_contacts = TRUE,  beta = beta_inf_amp_vl_on_contacts_on),
  list(het_vl = TRUE,  het_contacts = FALSE, beta = beta_inf_amp_vl_on_contacts_off),
  list(het_vl = FALSE, het_contacts = TRUE,  beta = beta_inf_amp_vl_off_contacts_on),
  list(het_vl = FALSE, het_contacts = FALSE, beta = beta_inf_amp_vl_off_contacts_off)
)

beta_inf_saved <- beta_inf
processed_infections_vl_sens <- map(heterogen_combos_amp, function(combo) {
  beta_inf <<- combo$beta
  ts <- traj_amplified %>%
    filter.(heterogen_vl == combo$het_vl) %>%
    select.(-m) %>%
    crossing.(prop_self_iso_test = c(0), sampling_freq = NA, event_size = NA) %>%
    mutate.(
      self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
      begin_testing = rdunif(n(), 0, sampling_freq)
    )
  run_model(
    testing_scenarios = ts,
    contact_dat       = contact_data,
    scenarios         = time_periods_base_amp %>% crossing(heterogen_contacts = combo$het_contacts),
    traj_full         = traj_amplified,
    traj_processed    = traj_amplified_,
    browsing          = F
  )
}) %>% bind_rows()
beta_inf <- beta_inf_saved
rm(beta_inf_saved, time_periods_base_amp, heterogen_combos_amp)

rm(traj_amplified, traj_amplified_, vl_params_amplified)

print("sens 2 done")
qsave(processed_infections_vl_sens, "results/processed_infections_vl_sens.qs")
rm(processed_infections_vl_sens)
gc()
#
### Additional analysis: testing effectiveness under heterogeneous vs homogeneous contacts ----
# Addresses reviewer request to show whether testing has a differential
# impact in a model with versus without contact heterogeneity.

time_periods_base_testing <- crossing(time_periods) %>%
  filter(period %in% c("Pre-pandemic", "1st lockdown", "School reopening")) %>%
  mutate(scenario_id = row_number()) %>%
  select(-c(date_start, date_end))

ts_testing_by_heterogen <- traj %>%
  filter.(heterogen_vl == T) %>%
  select.(-m) %>%
  crossing.(
    prop_self_iso_test = seq(0, 1, by = 0.1),
    sampling_freq = c(1, 3, 7),
    event_size = NA
  ) %>%
  mutate.(
    self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
    begin_testing = rdunif(n(), 0, sampling_freq)
  )

# Run separately for heterogen_contacts = T and F so each uses R0=2.5-calibrated beta
beta_inf_saved <- beta_inf
processed_infections_testing_by_heterogen <- map(
  list(
    list(het_contacts = TRUE,  beta = beta_inf),
    list(het_contacts = FALSE, beta = beta_inf_vl_on_contacts_off)
  ),
  function(combo) {
    beta_inf <<- combo$beta
    run_model(
      testing_scenarios = ts_testing_by_heterogen,
      contact_dat       = contact_data,
      scenarios         = time_periods_base_testing %>% crossing(heterogen_contacts = combo$het_contacts),
      browsing          = F
    )
  }
) %>% bind_rows()
beta_inf <- beta_inf_saved
rm(beta_inf_saved, time_periods_base_testing, ts_testing_by_heterogen)

print("testing by heterogen done")
qsave(processed_infections_testing_by_heterogen, "results/processed_infections_testing_by_heterogen.qs")
rm(processed_infections_testing_by_heterogen)
gc()

# all outputs saved above
