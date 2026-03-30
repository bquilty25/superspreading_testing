# run_new_analyses.R
# Runs only the two new analysis blocks added in revision:
#   1. Amplified VL heterogeneity sensitivity
#   2. Testing effectiveness by contact heterogeneity (heterogen vs homogeneous)
# Existing .qs outputs are loaded, not re-run.

source("scripts/utils.R")
source("scripts/duration.R")

N_sims <- 10000
set.seed(seed)

# ---- Rebuild base trajectories (needed for infctsnss_params) ----
traj <- vl_params %>%
    filter.(variant %in% c("wild")) %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = c(TRUE, FALSE)) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = N_sims, asymp_parms = asymp_fraction,
        variant_info = .x, browsing = F
    )) %>%
    bind_rows.()

infctsnss_params <- generate_params(culture_mod, N_sims) %>%
    as_tidytable(.) %>%
    rename.(beta0 = "(Intercept)", beta1 = vl) %>%
    mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

traj_ <- traj %>%
    mutate.(infectiousness = pmap(inf_curve_func, .l = list(
        m = m, start = start, end = end, interval = 1
    ))) %>%
    unnest.(infectiousness) %>%
    crossing.(lower_inf_thresh = c(FALSE)) %>%
    mutate.(
        culture_p = culture_prob(vl, beta0, beta1),
        infectious = rbernoulli(n = n(), p = pmin(culture_p * median_contact_duration, 1)),
        test_p = stats::predict(
            object = innova_mod, type = "response",
            newdata = tidytable(vl = vl)
        ),
        test = rbernoulli(n = n(), p = test_p),
        .by = c(lower_inf_thresh)
    ) %>%
    replace_na.(list(test = FALSE, infectious = FALSE)) %>%
    select.(-c(prolif, start, end))

# ---- 1. Amplified VL heterogeneity sensitivity ----
message("Running VL sensitivity (amplified SD) ...")

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
        variant_info = .x, browsing = F
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
        infectious = rbernoulli(n = n(), p = pmin(culture_p * median_contact_duration, 1)),
        test_p = stats::predict(
            object = innova_mod, type = "response",
            newdata = tidytable(vl = vl)
        ),
        test = rbernoulli(n = n(), p = test_p),
        .by = c(lower_inf_thresh)
    ) %>%
    replace_na.(list(test = FALSE, infectious = FALSE)) %>%
    select.(-c(prolif, start, end))

testing_scenarios <- traj_amplified %>%
    select.(-m) %>%
    crossing.(
        prop_self_iso_test = c(0),
        sampling_freq = NA,
        event_size = NA
    ) %>%
    mutate.(
        self_iso_test = rbernoulli(n = n(), prop_self_iso_test),
        begin_testing = rdunif(n(), 0, sampling_freq)
    )

time_periods_of_interest <- crossing(time_periods) %>%
    filter(date_end < as.Date("2021-01-01"), period != "POLYMOD") %>%
    mutate(scenario_id = row_number()) %>%
    select(-c(date_start, date_end)) %>%
    crossing(heterogen_contacts = c(T, F))

processed_infections_vl_sens <- run_model(
    testing_scenarios = testing_scenarios,
    contact_dat       = contact_data,
    scenarios         = time_periods_of_interest,
    traj_full         = traj_amplified,
    traj_processed    = traj_amplified_,
    browsing          = F
)
qsave(processed_infections_vl_sens, "results/processed_infections_vl_sens.qs")
message("VL sensitivity done.")

rm(testing_scenarios, time_periods_of_interest)

# ---- 2. Testing effectiveness by contact heterogeneity ----
message("Running testing × contact heterogeneity ...")

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

time_periods_of_interest <- crossing(time_periods) %>%
    filter(period %in% c("Pre-pandemic", "1st lockdown", "School reopening")) %>%
    mutate(scenario_id = row_number()) %>%
    select(-c(date_start, date_end)) %>%
    crossing(heterogen_contacts = c(T, F))

processed_infections_testing_by_heterogen <- run_model(
    testing_scenarios = testing_scenarios,
    contact_dat       = contact_data,
    scenarios         = time_periods_of_interest,
    browsing          = F
)
qsave(
    processed_infections_testing_by_heterogen,
    "results/processed_infections_testing_by_heterogen.qs"
)
message("Testing × heterogeneity done.")

message("All new analyses complete. Run run_new_figures.R to generate plots.")
