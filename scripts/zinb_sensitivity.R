# Sensitivity analysis: Zero-Inflated Negative Binomial vs Negative Binomial
# Fits NB and ZINB to the per-individual total daily contact distribution
# across 3 time periods.

source("scripts/utils.R")

dir.create("results/manuscript_figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/manuscript_tables", recursive = TRUE, showWarnings = FALSE)

colour_pal <- c("#E69F00", "#0072B2", "#009E73") # pre-pandemic=orange, lockdown=blue, school=green

periods_keep <- c("Pre-pandemic", "1st lockdown", "School reopening")

# ---- Load contact data -----------------------------------------------------
contacts <- contact_data %>%
    filter.(period %in% periods_keep) %>%
    mutate.(
        contacts = as.integer(round(e_all)),
        period = factor(period, levels = periods_keep)
    ) %>%
    select.(period, contacts)

cat("\nContact counts per period:\n")
contacts %>%
    summarise.(n = n(), mean_c = mean(contacts), zero_frac = mean(contacts == 0), .by = period) %>%
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
fit_results <- contacts %>%
    summarise.(
        fits = list(fit_both(contacts)),
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
cat("\n=== ZINB vs NB: daily contacts per individual ===\n\n")
fit_results %>%
    select.(
        period, n, obs_zero_frac,
        nb_mu, nb_k, zinb_mu, zinb_k, zinb_pi,
        delta_aic, delta_bic, lrt_p, sig
    ) %>%
    mutate.(across.(where(is.numeric), \(x) round(x, 4))) %>%
    as.data.frame() %>%
    print()

# ---- Save results table ----------------------------------------------------
table_zinb <- fit_results %>%
    select.(
        period, n,
        nb_mu, nb_k,
        zinb_mu, zinb_k, zinb_pi,
        delta_aic, lrt_p
    ) %>%
    mutate.(across.(where(is.numeric) & !all_of("zinb_pi"), \(x) round(x, 3)),
            zinb_pi = formatC(zinb_pi, format = "e", digits = 0))

table_zinb %>%
    as.data.frame() %>%
    print()

write.csv(table_zinb, "results/manuscript_tables/table_s2_zinb_sensitivity.csv", row.names = FALSE)
message("Saved table_s2_zinb_sensitivity.csv")
