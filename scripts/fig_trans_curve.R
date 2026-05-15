## Figure: Realised transmission probability over time (pre-pandemic)
## Two-panel: Household (left) | Non-household (right)
## Y-axis: Transmission probability
## X-axis: Time since symptom onset (days)

source("scripts/utils.R")
source("scripts/duration.R")

dir.create("results/manuscript_figures", recursive = TRUE, showWarnings = FALSE)

beta_inf <- qread("results/calibrated_beta.qs")

time_step <- 0.25
n_sims <- 5000

# --- 1. VL trajectories (wild-type, VL heterogeneity on) ----------------------
traj <- vl_params %>%
    filter.(variant == "wild") %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = TRUE) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = n_sims,
        asymp_parms = asymp_fraction,
        variant_info = .x,
        browsing = FALSE
    )) %>%
    bind_rows.()

infctsnss_params <- generate_params(culture_mod, n_sims) %>%
    as_tidytable() %>%
    rename.(beta0 = "(Intercept)", beta1 = vl) %>%
    mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

# --- 2. Per-simulation contact durations (sampled once per sim) ---------------
# Use pre-pandemic distributions for both HH and NHH
hh_durs <- contacts_hh_duration %>%
    filter.(period == "Pre-pandemic") %>%
    pull.(cnt_duration)
nhh_durs <- contacts_nhh_duration %>%
    filter.(period == "Pre-pandemic") %>%
    pull.(cnt_duration)

set.seed(seed)
sim_durs <- tidytable(
    sim          = 1:n_sims,
    hh_duration  = sample(hh_durs, n_sims, replace = TRUE),
    nhh_duration = sample(nhh_durs, n_sims, replace = TRUE)
)

# --- 3. Expand to time series -------------------------------------------------
plot_dat <- traj %>%
    mutate.(infectivity = pmap(inf_curve_func, .l = list(
        m        = m,
        start    = start,
        end      = end,
        interval = time_step
    ))) %>%
    unnest.(infectivity) %>%
    mutate.(
        culture_p  = culture_prob(vl, beta0, beta1),
        infectious = rbernoulli(n(), p = 1 - exp(-beta_inf * culture_p * median_contact_duration))
    ) %>%
    filter.(any(infectious), .by = sim) %>% # exclude never-infectious
    filter.(t >= 0, t <= 20) %>%
    left_join.(sim_durs, by = "sim") %>%
    mutate.(
        p_hh  = 1 - exp(-beta_inf * culture_p * hh_duration),
        p_nhh = 1 - exp(-beta_inf * culture_p * nhh_duration)
    ) %>%
    select.(-c(prolif, start, end))

# --- 4. Summarise over simulations -------------------------------------------
summary_dat <- plot_dat %>%
    mutate.(t_bin = round(t / time_step) * time_step) %>%
    summarise.(
        hh_mean = mean(p_hh),
        hh_lo = quantile(p_hh, 0.05),
        hh_hi = quantile(p_hh, 0.95),
        nhh_mean = mean(p_nhh),
        nhh_lo = quantile(p_nhh, 0.05),
        nhh_hi = quantile(p_nhh, 0.95),
        .by = t_bin
    )

# --- 5. Build panels ----------------------------------------------------------
make_panel <- function(dat, trans_col, trans_var, title_label) {
    lo_var <- paste0(trans_var, "_lo")
    hi_var <- paste0(trans_var, "_hi")
    mean_var <- paste0(trans_var, "_mean")

    ggplot(dat, aes(x = t_bin)) +
        geom_ribbon(
            aes(
                ymin = .data[[lo_var]],
                ymax = .data[[hi_var]],
                fill = "90% interval"
            ),
            alpha = 0.35,
            colour = NA
        ) +
        geom_line(
            aes(
                y        = .data[[mean_var]],
                linetype = "Mean"
            ),
            colour = trans_col,
            linewidth = 0.8
        ) +
        scale_fill_manual(
            name   = NULL,
            values = c("90% interval" = trans_col)
        ) +
        scale_linetype_manual(
            name   = NULL,
            values = c("Mean" = "dashed")
        ) +
        scale_y_continuous(
            name   = "Transmission probability",
            limits = c(0, 1),
            breaks = seq(0, 1, by = 0.25),
            labels = label_number(accuracy = 0.01)
        ) +
        scale_x_continuous(
            name   = "Days since first detectable by PCR (Ct<40)",
            breaks = breaks_width(5)
        ) +
        labs(title = title_label) +
        guides(
            fill     = guide_legend(order = 2, override.aes = list(alpha = 0.35)),
            linetype = guide_legend(order = 1, title = NULL)
        ) +
        plotting_theme +
        theme(
            plot.title       = element_text(hjust = 0.5, size = 10),
            legend.position  = "bottom",
            legend.key.width = unit(1, "cm")
        )
}

p_hh <- make_panel(summary_dat, bi_col_pal[1], "hh", "Household contacts")
p_nhh <- make_panel(summary_dat, bi_col_pal[2], "nhh", "Non-household contacts")

fig <- p_hh | p_nhh

ggsave(
    "results/manuscript_figures/fig3_trans_prob_curve.png",
    fig,
    dpi = 600, width = 210, height = 120, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig3_trans_prob_curve.pdf",
    fig,
    dpi = 600, width = 210, height = 120, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig3_trans_prob_curve.eps",
    fig,
    width = 210, height = 120, units = "mm", device = cairo_ps
)

message("Saved fig3_trans_prob_curve.{png,pdf,eps} to results/manuscript_figures/")

# --- Culture probability figure (Marc et al. style) --------------------------
culture_summary <- plot_dat %>%
    mutate.(t_bin = round(t / time_step) * time_step) %>%
    summarise.(
        cp_mean = mean(culture_p),
        cp_lo = quantile(culture_p, 0.05),
        cp_hi = quantile(culture_p, 0.95),
        .by = t_bin
    )

fig_culture <- ggplot(culture_summary, aes(x = t_bin)) +
    geom_ribbon(
        aes(ymin = cp_lo, ymax = cp_hi, fill = "90% interval"),
        alpha = 0.35, colour = NA
    ) +
    geom_line(
        aes(y = cp_mean, linetype = "Mean"),
        colour = bi_col_pal[1], linewidth = 0.8
    ) +
    scale_fill_manual(name = NULL, values = c("90% interval" = bi_col_pal[1])) +
    scale_linetype_manual(name = NULL, values = c("Mean" = "dashed")) +
    scale_y_continuous(
        name   = "Probability of culturing virus",
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25),
        labels = label_number(accuracy = 0.01)
    ) +
    scale_x_continuous(
        name   = "Days since first detectable by PCR (Ct<40)",
        breaks = breaks_width(5)
    ) +
    guides(
        fill     = guide_legend(order = 2, override.aes = list(alpha = 0.35)),
        linetype = guide_legend(order = 1)
    ) +
    plotting_theme +
    theme(legend.position = "bottom", legend.key.width = unit(1, "cm"))

ggsave(
    "results/manuscript_figures/fig_culture_prob_curve.png",
    fig_culture,
    dpi = 600, width = 120, height = 100, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_culture_prob_curve.pdf",
    fig_culture,
    dpi = 600, width = 120, height = 100, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_culture_prob_curve.eps",
    fig_culture,
    width = 120, height = 100, units = "mm", device = cairo_ps
)

message("Saved fig_culture_prob_curve.{png,pdf,eps} to results/manuscript_figures/")

# --- Amplified VL sensitivity analysis figure --------------------------------
# Doubles SDs of peak VL, proliferation, and clearance.
# max_peakvl capped at 40 (Ct detection limit) — same as baseline.

beta_inf_amp <- qread("results/calibrated_betas_amplified.qs")$vl_on_contacts_on

vl_params_amplified <- vl_params %>%
    mutate(
        sd_peakvl = sd_peakvl * 2,
        sd_prolif = sd_prolif * 2,
        sd_clear  = sd_clear * 2
    )

traj_amp <- vl_params_amplified %>%
    filter.(variant == "wild") %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = TRUE) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = n_sims,
        asymp_parms = asymp_fraction,
        variant_info = .x,
        max_prolif = 28,
        max_clear = 60,
        max_peakvl = 40, # reverted from 80
        browsing = FALSE
    )) %>%
    bind_rows.()

traj_amp <- traj_amp %>% left_join.(infctsnss_params, by = "sim")

set.seed(seed)
sim_durs_amp <- tidytable(
    sim          = 1:n_sims,
    hh_duration  = sample(hh_durs, n_sims, replace = TRUE),
    nhh_duration = sample(nhh_durs, n_sims, replace = TRUE)
)

plot_dat_amp <- traj_amp %>%
    mutate.(infectivity = pmap(inf_curve_func, .l = list(
        m        = m,
        start    = start,
        end      = end,
        interval = time_step
    ))) %>%
    unnest.(infectivity) %>%
    mutate.(
        culture_p  = culture_prob(vl, beta0, beta1),
        infectious = rbernoulli(n(), p = 1 - exp(-beta_inf_amp * culture_p * median_contact_duration))
    ) %>%
    filter.(any(infectious), .by = sim) %>%
    filter.(t >= 0, t <= 20) %>%
    left_join.(sim_durs_amp, by = "sim") %>%
    mutate.(
        p_hh  = 1 - exp(-beta_inf_amp * culture_p * hh_duration),
        p_nhh = 1 - exp(-beta_inf_amp * culture_p * nhh_duration)
    ) %>%
    select.(-c(prolif, start, end))

summary_dat_amp <- plot_dat_amp %>%
    mutate.(t_bin = round(t / time_step) * time_step) %>%
    summarise.(
        hh_mean = mean(p_hh),
        hh_lo = quantile(p_hh, 0.05),
        hh_hi = quantile(p_hh, 0.95),
        nhh_mean = mean(p_nhh),
        nhh_lo = quantile(p_nhh, 0.05),
        nhh_hi = quantile(p_nhh, 0.95),
        .by = t_bin
    )

p_hh_amp <- make_panel(summary_dat_amp, bi_col_pal[1], "hh", "Household contacts")
p_nhh_amp <- make_panel(summary_dat_amp, bi_col_pal[2], "nhh", "Non-household contacts")

fig_amp <- p_hh_amp | p_nhh_amp

ggsave(
    "results/manuscript_figures/fig_trans_prob_amplified.png",
    fig_amp,
    dpi = 600, width = 210, height = 120, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_trans_prob_amplified.pdf",
    fig_amp,
    dpi = 600, width = 210, height = 120, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_trans_prob_amplified.eps",
    fig_amp,
    width = 210, height = 120, units = "mm", device = cairo_ps
)

message("Saved fig_trans_prob_amplified.{png,pdf,eps} to results/manuscript_figures/")
