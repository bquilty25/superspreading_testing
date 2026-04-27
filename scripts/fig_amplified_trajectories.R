## Amplified VL heterogeneity: trajectory plots (Fig 2 style)
## Mirrors curve_plot.R but uses vl_params_amplified (doubled SDs, max_peakvl=40)

source("scripts/utils.R")
source("scripts/duration.R")

dir.create("results/manuscript_figures", recursive = TRUE, showWarnings = FALSE)

time_step <- 0.1
n_sims <- 5000

vl_params_amplified <- vl_params %>%
    mutate(
        sd_peakvl = sd_peakvl * 2,
        sd_prolif = sd_prolif * 2,
        sd_clear  = sd_clear * 2
    )

traj_normal <- vl_params %>%
    filter.(variant %in% c("wild")) %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = c(TRUE, FALSE)) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = n_sims,
        asymp_parms = asymp_fraction,
        variant_info = .x,
        browsing = F
    )) %>%
    bind_rows.() %>%
    mutate.(sd_type = "Normal SD")

traj_amplified <- vl_params_amplified %>%
    filter.(variant %in% c("wild")) %>%
    mutate.(variant = fct_drop(variant)) %>%
    crossing(heterogen_vl = c(TRUE, FALSE)) %>%
    group_split.(variant, heterogen_vl) %>%
    map.(~ make_trajectories(
        n_sims = n_sims,
        asymp_parms = asymp_fraction,
        variant_info = .x,
        max_prolif = 28, max_clear = 60, max_peakvl = 40,
        browsing = F
    )) %>%
    bind_rows.() %>%
    mutate.(sd_type = "Amplified SD")

traj <- bind_rows.(traj_normal, traj_amplified) %>%
    mutate.(sd_type = factor(sd_type, levels = c("Normal SD", "Amplified SD")))

infctsnss_params <- generate_params(culture_mod, n_sims) %>%
    as_tidytable() %>%
    rename.(beta0 = "(Intercept)", beta1 = vl) %>%
    mutate.(sim = row_number())

traj <- traj %>% left_join.(infctsnss_params, by = "sim")

# Integer-step data (for days_inf, auc)
plot_dat <- traj %>%
    mutate.(infectivity = pmap(inf_curve_func, .l = list(
        m = m, start = start, end = end
    ))) %>%
    unnest.(infectivity) %>%
    crossing.(lower_inf_thresh = c(FALSE)) %>%
    mutate.(
        culture_p = culture_prob(vl, beta0, beta1),
        infectious = rbernoulli(n = n(), p = culture_p),
        test_p = stats::predict(innova_mod,
            type = "response",
            newdata = tidytable(vl = vl)
        ),
        test = rbernoulli(n = n(), p = test_p),
        .by = c(lower_inf_thresh)
    ) %>%
    replace_na.(list(test = FALSE, infectious = FALSE)) %>%
    select.(-c(prolif, start, end))

# Fine-step data (for smooth curves)
plot_dat1 <- traj %>%
    mutate.(infectivity = pmap(inf_curve_func, .l = list(
        m = m, start = start, end = end, interval = time_step
    ))) %>%
    unnest.(infectivity) %>%
    crossing.(lower_inf_thresh = c(FALSE)) %>%
    mutate.(
        culture_p = culture_prob(vl, beta0, beta1),
        infectious = rbernoulli(n = n(), p = culture_p),
        test_p = stats::predict(innova_mod,
            type = "response",
            newdata = tidytable(vl = vl)
        ),
        test = rbernoulli(n = n(), p = test_p),
        .by = c(lower_inf_thresh)
    ) %>%
    replace_na.(list(test = FALSE, infectious = FALSE)) %>%
    select.(-c(prolif, start, end))

sd_colours <- c("Normal SD" = tri_col_pal[1], "Amplified SD" = bi_col_pal[2])

# Precompute median trajectories for panels A and C
log_median <- plot_dat1 %>%
    filter.(heterogen_vl == TRUE) %>%
    summarise.(median_vl = 10^median(vl), .by = c(t, sd_type))

culture_median <- plot_dat1 %>%
    filter.(heterogen_vl == TRUE) %>%
    summarise.(median_culture_p = median(culture_p), .by = c(t, sd_type))

# --- log_plot: spaghetti of VL trajectories -----------------------------------
log_plot <- plot_dat1 %>%
    filter.(sim <= 1000, heterogen_vl == TRUE) %>%
    mutate.(vl = 10^vl) %>%
    ggplot(aes(colour = sd_type)) +
    geom_line(
        aes(x = t, y = vl, group = interaction(sim, sd_type)),
        alpha = 0.04
    ) +
    geom_line(
        data = log_median,
        aes(x = t, y = median_vl, group = sd_type),
        linewidth = 1.2
    ) +
    geom_hline(yintercept = 1.6e7, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(
        name = "", values = sd_colours,
        guide = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2))
    ) +
    scale_x_continuous(
        name   = "Days since first detectable by PCR (Ct<40)",
        breaks = breaks_width(5)
    ) +
    scale_y_log10(name = "RNA copies/ml", labels = label_log()) +
    coord_cartesian() +
    plotting_theme

# --- culture_plot: spaghetti of culture probability curves -------------------
culture_plot <- plot_dat1 %>%
    filter.(sim <= 1000, heterogen_vl == TRUE) %>%
    ggplot(aes(colour = sd_type)) +
    geom_line(
        aes(x = t, y = culture_p, group = interaction(sim, sd_type)),
        alpha = 0.04
    ) +
    geom_line(
        data = culture_median,
        aes(x = t, y = median_culture_p, group = sd_type),
        linewidth = 1.2
    ) +
    scale_colour_manual(
        name = "", values = sd_colours,
        guide = "none"
    ) +
    ylab("Relative infectivity") +
    scale_x_continuous(
        name   = "Days since first detectable by PCR (Ct<40)",
        breaks = breaks_width(5)
    ) +
    plotting_theme

# --- inf_plot: culture probability vs VL -------------------------------------
prob_culture <- infctsnss_params %>%
    crossing(ct = seq(10, 40, by = 0.5)) %>%
    mutate.(
        vl        = convert_Ct_logGEML(ct),
        culture_p = culture_prob(vl, beta0, beta1)
    ) %>%
    arrange.(sim)

inf_plot <- prob_culture %>%
    mutate.(vl = 10^vl) %>%
    group_by(vl) %>%
    summarise(
        median_p = median(culture_p, na.rm = TRUE),
        lower    = quantile(culture_p, 0.025, na.rm = TRUE),
        upper    = quantile(culture_p, 0.975, na.rm = TRUE)
    ) %>%
    ggplot(aes(x = vl)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = bi_col_pal[1], alpha = 0.3) +
    geom_line(aes(y = median_p), colour = bi_col_pal[1], linewidth = 1) +
    geom_vline(xintercept = 1.6e7, linetype = "dashed") +
    scale_x_log10(name = "RNA copies/ml", labels = label_log(), limits = c(10^3.5, NA)) +
    labs(y = "Probability of\nculturing virus", caption = "Pickering et al. 2021") +
    plotting_theme

# --- days_inf_plot: distribution of infectious days --------------------------
days_inf <- plot_dat %>%
    summarise.(n_inf = sum(infectious == TRUE), .by = c(sim, heterogen_vl, sd_type))

days_inf_plot <- days_inf %>%
    filter.(heterogen_vl == TRUE) %>%
    mutate.(n_inf_days = floor(n_inf)) %>%
    filter.(n_inf_days >= 0) %>%
    count.(n_inf_days, sd_type) %>%
    mutate.(prop = n / sum(n), .by = sd_type) %>%
    ggplot(aes(x = n_inf_days, y = prop, fill = sd_type, colour = sd_type)) +
    geom_col(position = "identity", alpha = 0.5, width = 0.8) +
    scale_fill_manual(name = "", values = sd_colours, guide = "none") +
    scale_colour_manual(name = "", values = sd_colours, guide = "none") +
    ylab("Probability") +
    scale_x_continuous(
        name   = "Days infectious",
        breaks = breaks_width(1),
        limits = c(-0.5, NA),
        expand = expansion(add = c(0, 0.5))
    ) +
    plotting_theme

# --- auc_plot: distribution of area under infectivity curve ------------------
auc_dat <- plot_dat %>%
    filter.(heterogen_vl == TRUE) %>%
    summarise.(sum_inf = sum(culture_p), .by = c(sim, sd_type))

auc_plot <- auc_dat %>%
    ggplot(aes(x = sum_inf, fill = sd_type, colour = sd_type)) +
    geom_density(alpha = 0.25, adjust = 2) +
    scale_fill_manual(name = "", values = sd_colours, guide = "none") +
    scale_colour_manual(name = "", values = sd_colours, guide = "none") +
    lims(y = c(0, NA)) +
    labs(x = "Area under infectivity curve (AU)", y = "Probability density") +
    plotting_theme

# --- Combine and save --------------------------------------------------------
(log_plot | inf_plot | culture_plot) / (auc_plot | days_inf_plot) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(
    "results/manuscript_figures/fig2_vl_amplified.png",
    dpi = 600, width = 300, height = 150, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig2_vl_amplified.pdf",
    dpi = 600, width = 300, height = 150, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig2_vl_amplified.eps",
    width = 300, height = 150, units = "mm", device = cairo_ps
)

message("Saved fig2_vl_amplified.{png,pdf,eps} to results/manuscript_figures/")
