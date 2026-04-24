#### Results ####
source("scripts/utils.R")

dir.create("results/manuscript_figures", recursive = TRUE, showWarnings = FALSE)

# Safe fitdist wrapper: returns NA estimates instead of crashing on degenerate data
safe_nb_est <- function(x) {
  tryCatch(
    fitdist(x, "nbinom")$estimate,
    error = function(e) c(mu = NA_real_, size = NA_real_)
  )
}

# load simulation output
processed_infections_baseline <- qread("results/processed_infections_baseline.qs")
processed_infections_heterogen_on_off <- qread("results/processed_infections_heterogen_on_off.qs")
processed_infections_testing <- qread("results/processed_infections_testing.qs")
processed_infections_events <- qread("results/processed_infections_events.qs")
processed_infections_sens <- qread("results/processed_infections_sens.qs")
processed_infections_vl_sens <- qread("results/processed_infections_vl_sens.qs")
processed_infections_testing_by_heterogen <- qread("results/processed_infections_testing_by_heterogen.qs")

# inferred generation time
processed_infections_baseline %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(
      all_of(key_grouping_var[-1]),
      sampling_freq, prop_self_iso_test, t
    )
  ) %>%
  ggplot(aes(x = t, y = (sum_inf) / sum(sum_inf))) +
  geom_col(fill = bi_col_pal[2]) +
  stat_function(
    fun = dweibull,
    args = list(shape = 2.826, scale = 5.665),
    colour = bi_col_pal[1]
  ) +
  geom_text(
    aes(
      x = Inf, y = Inf,
      hjust = 1,
      vjust = 1,
      label = "Ferretti et al.\nWeibull shape: 2.826\nWeibull scale: 5.665"
    ),
    colour = bi_col_pal[1]
  ) +
  scale_x_continuous("Days since first detectable by PCR (Ct<40)",
    breaks = breaks_width(5),
    limits = c(0, 20)
  ) +
  scale_y_continuous("Proportion", labels = label_percent()) +
  plotting_theme

ggsave("results/gen_time.png", width = 210, height = 120, dpi = 600, units = "mm", bg = "white")

# R and K estimates over time
tic()
boot_res <- processed_infections_baseline %>%
  # filter(period=="Pre-pandemic") %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test)
  ) %>%
  summarise.(
    dists = list(bootdist(fitdist(sum_inf, "nbinom"),
      bootmethod = "nonparam",
      parallel = "multicore",
      ncpus = 8
    )$CI %>%
      as.data.frame() %>%
      rownames_to_column(var = "name") %>%
      rename(
        "lo" = `2.5%`,
        "hi" = `97.5%`
      )),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, -sim)
  )
toc()
qsave(boot_res, "results/R_and_k_bootstrap_ests.qs")

boot_res_sum <- boot_res %>%
  unnest.(dists) %>%
  filter.(variant == "wild") %>%
  mutate.(
    name = as.factor(name),
    name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0")
  ) %>%
  as_tibble()
write.csv(boot_res_sum, "results/R_and_k_bootstrap_ests.csv")

# # Derive 1/k rows (degree of overdispersion); CIs invert because 1/x is decreasing
# boot_res_sum <- boot_res_sum |>
#   bind_rows(
#     boot_res_sum |>
#       filter(name == "size") |>
#       mutate(
#         name = "inv_k",
#         tmp = lo,
#         lo = 1 / hi,
#         hi = 1 / tmp,
#         Median = 1 / Median
#       ) |>
#       select(-tmp)
#   ) |>
#   mutate(name = fct_relevel(name, "mu", "size", "inv_k", "prop_ss_10", "prop_ss_0"))

boot_res_sum %>%
  ggplot(aes(y = Median, ymin = lo, ymax = hi, x = period, colour = name, fill = name, group = name)) +
  geom_line() +
  geom_point() +
  geom_ribbon(alpha = 0.4, colour = NA) +
  geom_segment(
    data = rt_by_time_period %>% mutate(name = "mu"),
    aes(
      x = period,
      xend = period,
      y = lo,
      yend = hi
    ),
    alpha = 0.25,
    size = 10
  ) +
  facet_grid2(
    rows = vars(name), switch = "y", scales = "free",
    labeller = labeller(name = c(
      "mu" = "R", "size" = "k",
      # "inv_k" = "1/k (overdispersion)",
      "prop_ss_0" = "Proportion infecting\n 0 others (%)",
      "prop_ss_10" = "Proportion infecting\n over 10 others (%)"
    )),
    axes = "all",
    remove_labels = "x",
  ) +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, NA)),
    scale_y_log10(limits = c(NA, NA), breaks = log_breaks()),
    # scale_y_log10(limits = c(NA, NA), breaks = log_breaks()),
    scale_y_continuous(limits = c(0, NA)),
    scale_y_continuous(limits = c(0, NA))
  )) +
  geom_hline(aes(linetype = name, yintercept = 1), colour = quad_col_pal[1]) +
  scale_colour_manual(values = rep(bi_col_pal[1], 5), guide = "none") +
  scale_fill_manual(values = rep(bi_col_pal[1], 5), guide = "none") +
  scale_linetype_manual(values = c("dashed", NA, NA, NA, NA), guide = "none") +
  lims(y = c(0, NA)) +
  labs(
    y = "",
    x = "Time period"
  ) +
  plotting_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("results/manuscript_figures/fig3_rk.png", width = 150, height = 100, dpi = 600, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig3_rk.pdf", width = 150, height = 100, units = "mm", bg = "white")

# heterogen_on_off
tic()
boot_res_heterogen <- processed_infections_heterogen_on_off %>%
  # filter(period=="Pre-pandemic") %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test)
  ) %>%
  summarise.(
    dists = list(bootdist(fitdist(sum_inf, "nbinom"),
      bootmethod = "nonparam",
      parallel = "multicore",
      ncpus = 8
    )$CI %>%
      as.data.frame() %>%
      rownames_to_column(var = "name") %>%
      rename(
        "lo" = `2.5%`,
        "hi" = `97.5%`
      )),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, -sim)
  )
toc()
qsave(boot_res_heterogen, "results/R_and_k_bootstrap_ests_heterogen.qs")

# other_est <- tribble(~study,~xmin,~xmax,~ymin,~ymax, ~y,
#                      "Endo et al. 2020", -Inf, Inf, 0.05, 0.2, 0.1,
#                      "Riou & Althaus 2020", -Inf, Inf, 0.014, 6.95, 0.54,
#                      "Adam et al. 2020", -Inf, Inf, 0.45, 0.72, 0.58,
#                      "Laxminarayan et al. 2020", -Inf, Inf, 0.49, 0.52, 0.51)

add_heterogen_labels <- function(df) {
  df %>%
    mutate.(
      heterogen_label = case_when.(
        heterogen_vl & heterogen_contacts ~ "Variable viral load, overdispersed contacts",
        heterogen_vl & !heterogen_contacts ~ "Variable viral load, Poisson contacts",
        !heterogen_vl & heterogen_contacts ~ "Equal viral load, overdispersed contacts"
      ),
      heterogen_label = fct_relevel(
        heterogen_label,
        "Variable viral load, overdispersed contacts",
        "Variable viral load, Poisson contacts",
        "Equal viral load, overdispersed contacts"
      ),
      heterogen_vl = ifelse(heterogen_vl, "Heterogeneous viral load", "Homogeneous viral load"),
      heterogen_contacts = ifelse(heterogen_contacts, "Heterogeneous contacts", "Homogeneous contacts")
    )
}

boot_res_heterogen_sum <- boot_res_heterogen %>%
  unnest.(dists) %>%
  mutate.(name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0")) %>%
  filter.(
    variant == "wild",
    name %in% c("size"),
    !(heterogen_vl == FALSE & heterogen_contacts == FALSE)
  ) %>%
  add_heterogen_labels()

prop_dat_heterogen <- processed_infections_heterogen_on_off %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test)
  ) %>%
  summarise.(
    props = list({
      x <- sum_inf
      B <- 2000
      boot_ss10 <- replicate(B, mean(sample(x, length(x), replace = TRUE) > 10) * 100)
      boot_ss0  <- replicate(B, mean(sample(x, length(x), replace = TRUE) <= 0) * 100)
      tidytable(
        name    = c("prop_ss_10", "prop_ss_0"),
        Median  = c(mean(x > 10) * 100, mean(x <= 0) * 100),
        lo      = c(quantile(boot_ss10, 0.025), quantile(boot_ss0, 0.025)),
        hi      = c(quantile(boot_ss10, 0.975), quantile(boot_ss0, 0.975))
      )
    }),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, -sim)
  ) %>%
  unnest.(props) %>%
  filter.(
    variant == "wild",
    !(heterogen_vl == FALSE & heterogen_contacts == FALSE)
  ) %>%
  add_heterogen_labels()

boot_res_heterogen_sum <- bind_rows.(boot_res_heterogen_sum, prop_dat_heterogen) %>%
  mutate.(name = fct_relevel(name, "size", "prop_ss_10", "prop_ss_0"))

write.csv(boot_res_heterogen_sum, "results/R_and_k_bootstrap_ests_heterogen.csv")

(heterogen_plot <- (boot_res_heterogen_sum %>% ggplot(aes(y = Median, ymin = lo, ymax = hi, x = period, colour = name)) +
  geom_line(aes(colour = heterogen_label, fill = heterogen_label, group = heterogen_label, linetype = heterogen_label)) +
  geom_point(aes(colour = heterogen_label, fill = heterogen_label, group = heterogen_label, linetype = heterogen_label, shape = heterogen_label)) +
  geom_lineribbon(
    data = . %>% filter.(name == "size"),
    aes(colour = heterogen_label, fill = heterogen_label, group = heterogen_label, linetype = heterogen_label),
    alpha = 0.4
  ) +
  scale_colour_manual(values = c(bi_col_pal[1], bi_col_pal[2], bi_col_pal[1])) +
  scale_fill_manual(values = c(bi_col_pal[1], bi_col_pal[2], bi_col_pal[1])) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"), name = "") +
  scale_shape_manual(values = c(16, 17, 1), name = "") +
  labs(
    y = "",
    x = "Time period",
    linetype = "",
    colour = "",
    fill = ""
  ) +
  facet_grid2(name ~ .,
    switch = "y",
    independent = "y",
    scales = "free_y",
    axes = "all",
    remove_labels = "x",
    labeller = labeller(
      name = c(
        "mu" = "Mean R", "size" = "Overdispersion (k)",
        "prop_ss_0" = "Proportion infecting\n 0 others (%)",
        "prop_ss_10" = "Proportion infecting\n >10 others (%)"
      )
    )
  ) +
  ggh4x::facetted_pos_scales(y = list(
    scale_y_log10(expand = expansion(mult = 0.15)),
    scale_y_continuous(limits = c(0, NA)),
    scale_y_continuous(limits = c(0, NA))
  )) +
  plotting_theme +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.direction = "vertical",
    legend.position = "bottom"
  )
)
)


ggsave(heterogen_plot, file = "results/manuscript_figures/fig4_heterogen.png", width = 200, height = 220, dpi = 600, units = "mm", bg = "white")
ggsave(heterogen_plot, file = "results/manuscript_figures/fig4_heterogen.pdf", width = 200, height = 220, units = "mm", bg = "white")

#### Heterogeneity input distributions plot ----

het_vl_curves <- processed_infections_heterogen_on_off %>%
  filter.(period == "Pre-pandemic", heterogen_contacts == TRUE, sim <= 500) %>%
  mutate.(
    vl_label = ifelse(heterogen_vl, "Variable viral load", "Equal viral load (median trajectory)")
  ) %>%
  ggplot(aes(x = t, y = culture_p, group = interaction(sim, heterogen_vl))) +
  geom_line(
    data = . %>% filter.(heterogen_vl == TRUE),
    colour = bi_col_pal[1], alpha = 0.05
  ) +
  geom_line(
    data = . %>% filter.(heterogen_vl == FALSE, sim == 1),
    colour = bi_col_pal[2], linewidth = 1.2
  ) +
  scale_x_continuous(name = "Days since infection", breaks = breaks_width(5)) +
  scale_y_continuous(name = "Relative infectivity (culture probability)", limits = c(0, 1)) +
  annotate("text", x = Inf, y = Inf, label = "Variable VL (individual curves)",
    hjust = 1.1, vjust = 2, colour = bi_col_pal[1], size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = "Equal VL (median trajectory)",
    hjust = 1.1, vjust = 3.8, colour = bi_col_pal[2], size = 3.5) +
  plotting_theme

het_contacts_dat <- processed_infections_heterogen_on_off %>%
  filter.(period == "Pre-pandemic", heterogen_vl == TRUE) %>%
  summarise.(daily_contacts = mean(total_contacts), .by = c(sim, heterogen_contacts)) %>%
  mutate.(contacts_label = ifelse(heterogen_contacts, "Overdispersed contacts", "Poisson contacts (mean)"))

het_contacts_plot <- het_contacts_dat %>%
  ggplot(aes(x = daily_contacts, fill = contacts_label, colour = contacts_label)) +
  geom_density(alpha = 0.4, adjust = 1.5) +
  scale_x_continuous(name = "Mean daily contacts over infectious period", limits = c(0, NA)) +
  scale_y_continuous(name = "Density") +
  scale_fill_manual(values = c("Overdispersed contacts" = bi_col_pal[1],
                               "Poisson contacts (mean)" = bi_col_pal[2])) +
  scale_colour_manual(values = c("Overdispersed contacts" = bi_col_pal[1],
                                 "Poisson contacts (mean)" = bi_col_pal[2])) +
  labs(fill = "", colour = "") +
  plotting_theme +
  theme(legend.position = "bottom")

het_inputs_plot <- het_vl_curves / het_contacts_plot +
  plot_annotation(tag_levels = "A")

ggsave(het_inputs_plot, file = "results/manuscript_figures/fig_heterogen_inputs.png",
  width = 160, height = 200, dpi = 600, units = "mm", bg = "white")
ggsave(het_inputs_plot, file = "results/manuscript_figures/fig_heterogen_inputs.pdf",
  width = 160, height = 200, units = "mm", bg = "white")

#### Sensitivity analysis ----

res_sens <- processed_infections_baseline %>%
  filter.(period == "Pre-pandemic") %>%
  mutate.(contacts = "Unadjusted") %>%
  bind_rows.(processed_infections_sens %>%
    mutate.(contacts = "Adjusted")) %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, contacts)
  ) %>%
  summarise.(
    dists = list(fitdist(sum_inf, "nbinom")),
    dist_means = list(safe_nb_est(sum_inf) %>% enframe() %>% pivot_wider(names_from = name, values_from = value)),
    n = n(),
    ss_10 = sum(sum_inf > 10),
    ss_20 = sum(sum_inf > 20),
    ss_0 = sum(sum_inf <= 0),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, contacts, -sim)
  ) %>%
  mutate.(
    prop_ss_10 = ss_10 / n * 100,
    # prop_ss_20 =ss_20/n,
    prop_ss_0 = ss_0 / n * 100
  ) %>%
  unnest.(dist_means) %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>%
  mutate.(
    name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0"),
    contacts = fct_relevel(contacts, "Unadjusted")
  ) %>%
  filter.(variant == "wild") 

res_sens %>%
  ggplot(aes(y = value, x = contacts, colour = name, group = name, shape = contacts)) +
  geom_point(show.legend = F) +
  # geom_text_repel(data=. %>% filter.(period=="Pre-pandemic",name=="mu"),
  #                 aes(x=contacts,y=value,label=paste0("R0 = ",round(value,1))),family="Lato",
  #                 nudge_y = -0.2)+
  geom_hline(aes(linetype = name, yintercept = 1), colour = quad_col_pal[1]) +
  scale_colour_manual(values = quad_col_pal, guide = "none") +
  scale_linetype_manual(values = c("dashed", NA, NA, NA), guide = "none") +
  facet_grid2(~name,
    switch = "y", scales = "free_y", independent = "y",
    labeller = labeller(name = c(
      "mu" = "R", "size" = "k",
      "prop_ss_0" = "Proportion infecting\n 0 others (%)",
      "prop_ss_10" = "Proportion infecting\n over 10 others (%)"
    )),
    axes = "all",
    remove_labels = "x"
  ) +
  lims(y = c(0, NA)) +
  labs(
    y = "",
    x = "Contact data tail adjustment"
  ) +
  plotting_theme +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("results/manuscript_figures/fig_sensitivity.png", width = 200, height = 100, dpi = 600, units = "mm", bg = "white")

processed_infections_baseline %>%
  # filter.(prop_self_iso_test==0,sampling_freq==3) %>%
  ggplot(aes(x = vl, y = total_contacts, colour = total_infections)) +
  geom_jitter(alpha = 0.5) +
  scale_y_log10("Daily contacts") +
  labs(x = "Viral load (RNA copies/ml)") +
  scale_colour_viridis_c("Infected contacts", trans = "log10", na.value = NA, option = "turbo") +
  facet_grid(prop_self_iso_test ~ period + sampling_freq) +
  plotting_theme

# ggsave("contacts_infections.png")

#### LFT testing ----

## regular testing ----

res_testing <- processed_infections_testing %>%
  summarise.(sum_inf = sum(total_infections), .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, event_size)) %>%
  summarise.(
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, event_size, -sim),
    dist_means = list(safe_nb_est(sum_inf) %>% enframe() %>% pivot_wider(names_from = name, values_from = value)),
    n = n(),
    ss_10 = sum(sum_inf > 10),
    ss_0 = sum(sum_inf <= 0)
  ) %>%
  mutate.(
    prop_ss_10 = ss_10 / n * 100,
    prop_ss_0 = ss_0 / n * 100
  ) %>%
  unnest.(dist_means) %>%
  filter.(variant == "wild") %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>%
  mutate.(name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0")) %>%
  ggplot(aes(y = value, x = prop_self_iso_test * 100, colour = factor(sampling_freq), group = sampling_freq, linetype = factor(sampling_freq))) +
  # geom_point()+
  geom_line() +
  geom_hline(data = ~ .x[.x$name == "mu", ], aes(yintercept = 1), colour = bi_col_pal[1], linetype = "dashed") +
  scale_colour_manual(values = tri_col_pal) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Testing frequency (days between tests)") +
  facet_grid2(name ~ period,
    # scales="free",
    scales = "free_y",
    remove_labels = "x",
    axes = "all",
    # independent = "y",
    labeller = labeller(name = c(
      "mu" = "R", "size" = "k",
      "prop_ss_0" = "Proportion infecting\n 0 others (%)",
      "prop_ss_10" = "Proportion infecting\n over 10 others (%)"
    )),
    switch = "y"
  ) +
  ggh4x::facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 3)),
    scale_y_log10(),
    scale_y_continuous(limits = c(0, NA)),
    scale_y_continuous(limits = c(0, NA))
  )) +
  lims(y = c(0, NA)) +
  labs(
    y = "",
    x = "Uptake of/adherence to lateral flow testing (%)",
    colour = "Testing frequency (days between tests)"
  ) +
  plotting_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

ggsave("results/lft_impact_testing.png", width = 210, height = 150, dpi = 600, units = "mm", bg = "white")
ggsave("results/lft_impact_testing.pdf", width = 210, height = 150, dpi = 600, units = "mm", bg = "white")

## events ----

res_events <- processed_infections_events %>%
  summarise.(sum_inf = sum(total_infections), .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, event_size)) %>%
  summarise.(
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, event_size, -sim),
    dist_means = list(safe_nb_est(sum_inf) %>%
      enframe() %>%
      pivot_wider(names_from = name, values_from = value)),
    n = n(),
    ss_10 = sum(sum_inf > 10),
    ss_0 = sum(sum_inf <= 0)
  ) %>%
  mutate.(
    prop_ss_10 = ss_10 / n * 100,
    prop_ss_0 = ss_0 / n * 100
  ) %>%
  unnest.(dist_means) %>%
  drop_na.(event_size) %>%
  filter.(variant == "wild") %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>%
  mutate.(name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0")) %>%
  ggplot(aes(y = value, x = prop_self_iso_test * 100, colour = factor(event_size), group = event_size, linetype = factor(event_size))) +
  # geom_point()+
  geom_line() +
  geom_hline(data = ~ .x[.x$name == "mu", ], aes(yintercept = 1), colour = bi_col_pal[1], linetype = "dashed") +
  scale_colour_manual(values = tri_col_pal) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Minimum event size\nto prompt testing") +
  facet_grid2(name ~ period,
    # scales="free",
    scales = "free_y",
    remove_labels = "x",
    # independent = "y",
    axes = "all",
    labeller = labeller(name = c(
      "mu" = "R", "size" = "k",
      "prop_ss_0" = "Proportion infecting\n 0 others (%)",
      "prop_ss_10" = "Proportion infecting\n over 10 others (%)"
    )),
    switch = "y"
  ) +
  ggh4x::facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 3)),
    scale_y_log10(),
    scale_y_continuous(limits = c(0, NA)),
    scale_y_continuous(limits = c(0, NA))
  )) +
  # lims(y=c(0,NA))+
  labs(
    y = "Mean parameter value",
    x = "Uptake of/adherence to pre-event lateral flow testing (%)",
    colour = "Minimum event size\nto prompt testing"
  ) +
  plotting_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

ggsave("results/lft_impact_events.png", width = 210, height = 150, dpi = 600, units = "mm", bg = "white")
ggsave("results/lft_impact_events.pdf", width = 210, height = 150, units = "mm", bg = "white")


testing_plot / events_plot + plot_annotation(tag_levels = "A")
ggsave("results/manuscript_figures/fig5_testing.png", dpi = 600, width = 210, height = 325, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig5_testing.pdf", width = 210, height = 300, units = "mm", bg = "white")

#### Sensitivity analysis: amplified VL heterogeneity ----

boot_res_vl_sens <- bind_rows(
  processed_infections_heterogen_on_off %>% mutate.(vl_sd_multiplier = "Standard (1\u00d7 SD)"),
  processed_infections_vl_sens %>% mutate.(vl_sd_multiplier = "Amplified (2\u00d7 SD)")
) %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, vl_sd_multiplier)
  ) %>%
  summarise.(
    dists = list(bootdist(fitdist(sum_inf, "nbinom"),
      bootmethod = "nonparam", parallel = "multicore", ncpus = 8
    )$CI %>%
      as.data.frame() %>%
      rownames_to_column(var = "name") %>%
      rename("lo" = `2.5%`, "hi" = `97.5%`)),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, vl_sd_multiplier, -sim)
  ) %>%
  unnest.(dists) %>%
  filter.(
    variant == "wild", name == "size",
    !(heterogen_vl == FALSE & heterogen_contacts == FALSE)
  ) %>%
  mutate.(
    heterogen_label = case_when.(
      heterogen_vl & heterogen_contacts ~ "Variable VL, overdispersed contacts",
      heterogen_vl & !heterogen_contacts ~ "Variable VL, Poisson contacts",
      !heterogen_vl & heterogen_contacts ~ "Equal VL, overdispersed contacts",
      !heterogen_vl & !heterogen_contacts ~ "Equal VL, Poisson contacts"
    ),
    heterogen_label = fct_relevel(
      heterogen_label,
      "Variable VL, overdispersed contacts",
      "Variable VL, Poisson contacts",
      "Equal VL, overdispersed contacts"
    ),
    vl_sd_multiplier = fct_relevel(vl_sd_multiplier, "Standard (1\u00d7 SD)")
  )
qsave(boot_res_vl_sens, "results/k_bootstrap_ests_vl_sens.qs")
write.csv(boot_res_vl_sens, "results/k_bootstrap_ests_vl_sens.csv")

vl_sens_plot <- boot_res_vl_sens %>%
  ggplot(aes(
    y = Median, ymin = lo, ymax = hi, x = period,
    colour = heterogen_label, fill = heterogen_label,
    group = heterogen_label, linetype = heterogen_label
  )) +
  geom_line() +
  geom_point(aes(shape = heterogen_label)) +
  geom_lineribbon(alpha = 0.25) +
  facet_wrap(~vl_sd_multiplier, ncol = 2) +
  scale_colour_manual(values = c(bi_col_pal[1], bi_col_pal[2], bi_col_pal[1])) +
  scale_fill_manual(values = c(bi_col_pal[1], bi_col_pal[2], bi_col_pal[1])) +
  scale_linetype_manual(values = c("solid", "solid", "dashed"), name = "") +
  scale_shape_manual(values = c(16, 17, 1), name = "") +
  scale_y_log10() +
  labs(
    y = "Overdispersion (k)",
    x = "Time period",
    colour = "", fill = "", linetype = "",
    caption = "Left: VL parameters as estimated from Kissler et al.\nRight: SDs of peak Ct, proliferation, and clearance all doubled."
  ) +
  plotting_theme +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.direction = "vertical", legend.position = "bottom"
  )

ggsave(vl_sens_plot,
  file = "results/manuscript_figures/fig_heterogen_vl_sens.png",
  width = 280, height = 150, dpi = 600, units = "mm", bg = "white"
)
ggsave(vl_sens_plot,
  file = "results/manuscript_figures/fig_heterogen_vl_sens.pdf",
  width = 280, height = 150, units = "mm", bg = "white"
)
#
#### Additional figure: testing effectiveness by contact heterogeneity ----

testing_heterogen_plot <- processed_infections_testing_by_heterogen %>%
  summarise.(
    sum_inf = sum(total_infections),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test)
  ) %>%
  summarise.(
    dist_means = list(safe_nb_est(sum_inf) %>%
      enframe() %>% pivot_wider(names_from = name, values_from = value)),
    n = n(),
    ss_10 = sum(sum_inf > 10),
    ss_0 = sum(sum_inf <= 0),
    .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test, -sim)
  ) %>%
  mutate.(prop_ss_10 = ss_10 / n * 100, prop_ss_0 = ss_0 / n * 100) %>%
  unnest.(dist_means) %>%
  pivot_longer.(c(prop_ss_10, prop_ss_0, size, mu)) %>%
  mutate.(
    name = fct_relevel(name, "mu", "size", "prop_ss_10", "prop_ss_0"),
    contact_label = ifelse(heterogen_contacts, "Overdispersed contacts", "Poisson contacts")
  ) %>%
  filter.(variant == "wild") %>%
  drop_na.(sampling_freq) %>%
  ggplot(aes(
    y = value, x = prop_self_iso_test * 100,
    colour = factor(sampling_freq), linetype = contact_label,
    group = interaction(sampling_freq, contact_label)
  )) +
  geom_line() +
  geom_hline(
    data = . %>% filter.(name == "mu"),
    aes(yintercept = 1), colour = quad_col_pal[1], linetype = "dashed"
  ) +
  scale_colour_manual(values = tri_col_pal) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  facet_grid2(name ~ period,
    scales = "free_y", remove_labels = "x", axes = "all",
    labeller = labeller(name = c(
      "mu"          = "R",
      "size"        = "k",
      "prop_ss_0"   = "Proportion infecting\n0 others (%)",
      "prop_ss_10"  = "Proportion infecting\n>10 others (%)"
    )),
    switch = "y"
  ) +
  ggh4x::facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 3)),
    scale_y_log10(),
    scale_y_continuous(limits = c(0, NA)),
    scale_y_continuous(limits = c(0, NA))
  )) +
  labs(
    y = "",
    x = "Uptake / adherence (%)",
    colour = "Testing frequency\n(days between tests)",
    linetype = "Contact distribution"
  ) +
  plotting_theme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

ggsave(testing_heterogen_plot,
  file = "results/manuscript_figures/fig_testing_heterogen.png",
  width = 280, height = 200, dpi = 600, units = "mm", bg = "white"
)
ggsave(testing_heterogen_plot,
  file = "results/manuscript_figures/fig_testing_heterogen.pdf",
  width = 280, height = 200, units = "mm", bg = "white"
)
