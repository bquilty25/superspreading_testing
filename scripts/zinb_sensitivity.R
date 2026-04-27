# Sensitivity analysis: Zero-Inflated Negative Binomial vs Negative Binomial
# Fits NB and ZINB to the per-individual offspring distribution (secondary
# infections) from the baseline model across 3 time periods.

source("scripts/utils.R")

dir.create("results/manuscript_figures", recursive = TRUE, showWarnings = FALSE)

colour_pal <- c("#E69F00", "#0072B2", "#009E73") # pre-pandemic=orange, lockdown=blue, school=green

periods_keep <- c("Pre-pandemic", "1st lockdown", "School reopening")

# ---- Load baseline simulation output ---------------------------------------
processed_infections_baseline <- qread("results/processed_infections_baseline.qs")

# ---- Derive per-individual offspring counts --------------------------------
# Each row in processed_infections_baseline is one infectious individual on one
# day; sum total_infections across days to get the full offspring count per
# individual. Use the no-testing baseline scenario only
# (sampling_freq=Inf, prop_self_iso_test=0).
offspring <- processed_infections_baseline %>%
    filter.(period %in% periods_keep) %>%
    summarise.(
        offspring = sum(total_infections),
        .by = c(sim, idx_id, period)
    ) %>%
    mutate.(period = factor(period, levels = periods_keep))

cat("\nOffspring counts per period:\n")
offspring %>%
    summarise.(n = n(), mean_off = mean(offspring), zero_frac = mean(offspring == 0), .by = period) %>%
    as.data.frame() %>%
    print()

# ---- ZINB log-likelihood (per-observation, unconstrained params) ----------
ll_zinb <- function(params, x) {
    mu <- exp(params[1])
    k <- exp(params[2])
    pi <- plogis(params[3])
    lp <- ifelse(
        x == 0,
        log(pi + (1 - pi) * dnbinom(0L, size = k, mu = mu)),
        log(1 - pi) + dnbinom(x, size = k, mu = mu, log = TRUE)
    )
    -sum(lp[is.finite(lp)])
}

# ---- Fit both models -------------------------------------------------------
fit_both <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)

    # NB via fitdist (reliable profile-likelihood MLE)
    nb_fit <- tryCatch(
        suppressWarnings(fitdist(x, "nbinom")),
        error = function(e) NULL
    )
    if (is.null(nb_fit)) {
        return(NULL)
    }

    nb_mu <- unname(nb_fit$estimate["mu"])
    nb_k <- unname(nb_fit$estimate["size"])
    nb_ll <- nb_fit$loglik
    nb_aic <- nb_fit$aic
    nb_bic <- nb_fit$bic

    # ZINB warm-started from NB
    zinb_opt <- tryCatch(
        optim(
            c(log(nb_mu), log(nb_k), qlogis(0.05)),
            ll_zinb,
            x = x,
            method = "L-BFGS-B",
            control = list(maxit = 2000)
        ),
        error = function(e) NULL
    )
    if (is.null(zinb_opt) || zinb_opt$convergence != 0) {
        return(NULL)
    }

    zinb_mu <- exp(zinb_opt$par[1])
    zinb_k <- exp(zinb_opt$par[2])
    zinb_pi <- plogis(zinb_opt$par[3])
    zinb_ll <- -zinb_opt$value
    zinb_np <- 3L
    zinb_aic <- -2 * zinb_ll + 2 * zinb_np
    zinb_bic <- -2 * zinb_ll + zinb_np * log(n)

    lrt_stat <- 2 * (zinb_ll - nb_ll)
    lrt_p <- pchisq(lrt_stat, df = 1, lower.tail = FALSE)

    nb_pred_zero <- dnbinom(0L, size = nb_k, mu = nb_mu)
    zinb_pred_zero <- zinb_pi + (1 - zinb_pi) * dnbinom(0L, size = zinb_k, mu = zinb_mu)
    obs_zero <- mean(x == 0)

    tibble(
        n = n,
        obs_zero_frac = obs_zero,
        nb_mu = nb_mu, nb_k = nb_k,
        nb_ll = nb_ll, nb_aic = nb_aic, nb_bic = nb_bic,
        nb_pred_zero = nb_pred_zero,
        zinb_mu = zinb_mu, zinb_k = zinb_k, zinb_pi = zinb_pi,
        zinb_ll = zinb_ll, zinb_aic = zinb_aic, zinb_bic = zinb_bic,
        zinb_pred_zero = zinb_pred_zero,
        delta_aic = zinb_aic - nb_aic,
        delta_bic = zinb_bic - nb_bic,
        lrt_stat = lrt_stat,
        lrt_p = lrt_p
    )
}

# ---- Run fits per period ---------------------------------------------------
fit_results <- offspring %>%
    summarise.(
        fits = list(fit_both(offspring)),
        .by = period
    ) %>%
    filter.(!sapply(fits, is.null)) %>%
    unnest.(fits) %>%
    mutate.(
        sig = case_when(
            lrt_p < 0.001 ~ "***",
            lrt_p < 0.01 ~ "**",
            lrt_p < 0.05 ~ "*",
            TRUE ~ "n.s."
        )
    )

# ---- Print summary table ---------------------------------------------------
cat("\n=== ZINB vs NB: secondary infections per individual ===\n\n")
fit_results %>%
    select.(
        period, n, obs_zero_frac,
        nb_mu, nb_k, zinb_mu, zinb_k, zinb_pi,
        delta_aic, delta_bic, lrt_p, sig
    ) %>%
    mutate.(across.(where(is.numeric), \(x) round(x, 4))) %>%
    as.data.frame() %>%
    print()

# ---- Figure A: ΔAIC --------------------------------------------------------
p_delta_aic <- fit_results %>%
    ggplot(aes(x = period, y = delta_aic, fill = period)) +
    geom_col(width = 0.5, colour = "grey30") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -2, linetype = "dotted", colour = "grey60") +
    geom_text(
        aes(
            label = sig,
            y = delta_aic + sign(delta_aic) * abs(delta_aic) * 0.05 + sign(delta_aic) * 0.3
        ),
        size = 3.5
    ) +
    scale_fill_manual(values = colour_pal, guide = "none") +
    scale_x_discrete(limits = periods_keep) +
    labs(
        x = "",
        y = "\u0394AIC (ZINB \u2212 NB)"
    ) +
    plotting_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

# ---- Figure B: Observed vs predicted zero fraction -------------------------
zero_df <- fit_results %>%
    select.(period, obs_zero_frac, nb_pred_zero, zinb_pred_zero) %>%
    pivot_longer.(
        cols = c(obs_zero_frac, nb_pred_zero, zinb_pred_zero),
        names_to = "model", values_to = "zero_frac"
    ) %>%
    mutate.(
        model = factor(model,
            levels = c("obs_zero_frac", "nb_pred_zero", "zinb_pred_zero"),
            labels = c("Observed", "NB predicted", "ZINB predicted")
        ),
        period = factor(period, levels = periods_keep)
    )

p_zeros <- zero_df %>%
    ggplot(aes(
        x = period, y = zero_frac * 100,
        colour = model, shape = model, linetype = model, group = model
    )) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 3) +
    scale_colour_manual(
        name = "",
        values = c("Observed" = "#333333", "NB predicted" = "#0072B2", "ZINB predicted" = "#D55E00")
    ) +
    scale_shape_manual(
        name = "",
        values = c("Observed" = 16, "NB predicted" = 17, "ZINB predicted" = 15)
    ) +
    scale_linetype_manual(
        name = "",
        values = c("Observed" = "solid", "NB predicted" = "dashed", "ZINB predicted" = "dotted")
    ) +
    scale_y_continuous(labels = label_percent(scale = 1)) +
    scale_x_discrete(limits = periods_keep) +
    labs(x = "", y = "Individuals causing zero\nsecondary infections (%)") +
    plotting_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

# ---- Figure C: NB k and mu by period --------------------------------------
params_df <- fit_results %>%
    select.(period, nb_mu, nb_k, zinb_mu, zinb_k) %>%
    pivot_longer.(cols = -period) %>%
    mutate.(
        model = if_else(grepl("^nb", name), "NB", "ZINB"),
        param = if_else(grepl("mu", name), "mu (mean)", "k (dispersion)"),
        period = factor(period, levels = periods_keep)
    )

p_params <- params_df %>%
    ggplot(aes(x = period, y = value, colour = model, shape = model, linetype = model, group = model)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 3) +
    facet_wrap2(~param, scales = "free_y", axes = "all") +
    scale_colour_manual(name = "Model", values = c("NB" = "#0072B2", "ZINB" = "#D55E00")) +
    scale_shape_manual(name = "Model", values = c("NB" = 16, "ZINB" = 17)) +
    scale_linetype_manual(name = "Model", values = c("NB" = "solid", "ZINB" = "dashed")) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = "", y = "") +
    plotting_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

# ---- Combine and save ------------------------------------------------------
fig_zinb <- (p_delta_aic | p_zeros) / p_params +
    plot_annotation(
        tag_levels = "A",
        caption = "Negative \u0394AIC => ZINB preferred; dotted line = \u22122\n* p<0.05  ** p<0.01  *** p<0.001 (LRT, df=1)"
    ) +
    plot_layout(heights = c(1, 1))

ggsave(
    "results/manuscript_figures/fig_zinb_sensitivity.png",
    fig_zinb,
    dpi = 600, width = 220, height = 210, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_zinb_sensitivity.pdf",
    fig_zinb,
    width = 220, height = 210, units = "mm", bg = "white"
)
ggsave(
    "results/manuscript_figures/fig_zinb_sensitivity.eps",
    fig_zinb,
    width = 220, height = 210, units = "mm", device = cairo_ps
)

fit_results %>%
    mutate.(across.(where(is.numeric), \(x) round(x, 4))) %>%
    write.csv("results/zinb_sensitivity_results.csv", row.names = FALSE)

message("Saved fig_zinb_sensitivity.{png,pdf,eps} and zinb_sensitivity_results.csv")
