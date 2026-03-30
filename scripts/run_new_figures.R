# run_new_figures.R
# Generates figures from the two new analysis outputs and places them in figures/

source("scripts/utils.R")

processed_infections_heterogen_on_off <- qread("results/processed_infections_heterogen_on_off.qs")
processed_infections_vl_sens <- qread("results/processed_infections_vl_sens.qs")
processed_infections_testing_by_heterogen <- qread("results/processed_infections_testing_by_heterogen.qs")

#### Figure: VL sensitivity (standard vs amplified SD) ----

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

vl_sens_plot <- boot_res_vl_sens %>%
    ggplot(aes(
        y = Median, ymin = lo, ymax = hi, x = period,
        colour = heterogen_label, fill = heterogen_label,
        group = heterogen_label, linetype = heterogen_label
    )) +
    geom_line() +
    geom_point() +
    geom_lineribbon(alpha = 0.25) +
    facet_wrap(~vl_sd_multiplier, ncol = 2) +
    scale_colour_manual(values = quad_col_pal[1:3]) +
    scale_fill_manual(values = quad_col_pal[1:3]) +
    scale_y_log10() +
    coord_cartesian(ylim = c(0.1, 10)) +
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
    file = "figures/fig_heterogen_vl_sens.png",
    width = 280, height = 150, dpi = 600, units = "mm", bg = "white"
)
ggsave(vl_sens_plot,
    file = "figures/fig_heterogen_vl_sens.pdf",
    width = 280, height = 150, units = "mm", bg = "white"
)
message("VL sensitivity figure saved.")

#### Figure: testing effectiveness by contact heterogeneity ----

testing_heterogen_plot <- processed_infections_testing_by_heterogen %>%
    summarise.(
        sum_inf = sum(total_infections),
        .by = c(all_of(key_grouping_var), sampling_freq, prop_self_iso_test)
    ) %>%
    summarise.(
        dist_means = list(fitdist(sum_inf, "nbinom")$estimate %>%
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
    file = "figures/fig_testing_heterogen.png",
    width = 280, height = 200, dpi = 600, units = "mm", bg = "white"
)
ggsave(testing_heterogen_plot,
    file = "figures/fig_testing_heterogen.pdf",
    width = 280, height = 200, units = "mm", bg = "white"
)
message("Testing × heterogeneity figure saved.")
