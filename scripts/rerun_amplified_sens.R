## Rerun amplified VL sensitivity analysis only (max_peakvl = 40)
source("scripts/utils.R")
source("scripts/duration.R")

N_sims <- 10000
beta_inf <- qread("results/calibrated_beta.qs")

# Rebuild baseline trajectories (needed for infctsnss_params)
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
    as_tidytable() %>%
    rename.(beta0 = "(Intercept)", beta1 = vl) %>%
    mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

### Sensitivity analysis: amplified VL heterogeneity (max_peakvl = 40) --------
vl_params_amplified <- vl_params %>%
    mutate(
        sd_peakvl = sd_peakvl * 2,
        sd_prolif = sd_prolif * 2,
        sd_clear  = sd_clear * 2
    )

traj_amplified <- vl_params_amplified %>%
    filter.(variant %in% c("wild")) %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = c(TRUE, FALSE)) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = N_sims, asymp_parms = asymp_fraction,
        variant_info = .x,
        max_prolif = 28, max_clear = 60, max_peakvl = 40,
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

# Calibrate betas
message("Calibrating beta_inf for amplified VL (het_vl=T, het_contacts=T) ...")
beta_inf_amp_vl_on_contacts_on <- calibrate_beta(
    target_R0 = 2.5, n_calib = 10000,
    heterogen_vl_flag = TRUE, heterogen_contacts_flag = TRUE, traj_data = traj_amplified
)
message(sprintf("  = %.4f", beta_inf_amp_vl_on_contacts_on))

message("Calibrating beta_inf for amplified VL (het_vl=T, het_contacts=F) ...")
beta_inf_amp_vl_on_contacts_off <- calibrate_beta(
    target_R0 = 2.5, n_calib = 10000,
    heterogen_vl_flag = TRUE, heterogen_contacts_flag = FALSE, traj_data = traj_amplified
)
message(sprintf("  = %.4f", beta_inf_amp_vl_on_contacts_off))

message("Calibrating beta_inf for amplified VL (het_vl=F, het_contacts=T) ...")
beta_inf_amp_vl_off_contacts_on <- calibrate_beta(
    target_R0 = 2.5, n_calib = 10000,
    heterogen_vl_flag = FALSE, heterogen_contacts_flag = TRUE, traj_data = traj_amplified
)
message(sprintf("  = %.4f", beta_inf_amp_vl_off_contacts_on))

message("Calibrating beta_inf for amplified VL (het_vl=F, het_contacts=F) ...")
beta_inf_amp_vl_off_contacts_off <- calibrate_beta(
    target_R0 = 2.5, n_calib = 10000,
    heterogen_vl_flag = FALSE, heterogen_contacts_flag = FALSE, traj_data = traj_amplified
)
message(sprintf("  = %.4f", beta_inf_amp_vl_off_contacts_off))

qsave(
    list(
        vl_on_contacts_on = beta_inf_amp_vl_on_contacts_on,
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
    list(het_vl = TRUE, het_contacts = TRUE, beta = beta_inf_amp_vl_on_contacts_on),
    list(het_vl = TRUE, het_contacts = FALSE, beta = beta_inf_amp_vl_on_contacts_off),
    list(het_vl = FALSE, het_contacts = TRUE, beta = beta_inf_amp_vl_off_contacts_on),
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

print("sens 2 done")
qsave(processed_infections_vl_sens, "results/processed_infections_vl_sens.qs")
message("Saved processed_infections_vl_sens.qs")
