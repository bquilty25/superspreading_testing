source("scripts/utils.R")

cat("\n=== Proportion e_other >= 250 by period ===\n")
contact_data %>%
    filter(period %in% c("POLYMOD", "Relaxed restrictions", "School reopening", "Step 2 + schools", "Pre-pandemic")) %>%
    group_by(period) %>%
    summarise(n = n(), over_250 = sum(e_other >= 250), prop = over_250 / n) %>%
    print()

cat("\n=== Summary stats for e_other by period ===\n")
contact_data %>%
    filter(period %in% c("POLYMOD", "Relaxed restrictions", "School reopening", "Step 2 + schools", "Pre-pandemic")) %>%
    group_by(period) %>%
    summarise(
        n = n(), mean = mean(e_other), median = median(e_other),
        p95 = quantile(e_other, 0.95), p99 = quantile(e_other, 0.99),
        max = max(e_other)
    ) %>%
    print()

cat("\n=== GPD fit: CoMix non-lockdown (current approach) ===\n")
excess_comix <- contact_data %>%
    filter(period %in% c("Relaxed restrictions", "School reopening", "Step 2 + schools")) %>%
    filter(e_other >= 250) %>%
    pull(e_other) - 250
cat("n exceedances:", length(excess_comix), "\n")
fit_comix <- fitdistrplus::fitdist(
    excess_comix, "gpd",
    fix.arg = list(mu = 0),
    start = list(sigma = mean(excess_comix), xi = 0.5),
    lower = c(1e-6, -Inf)
)
cat("sigma =", fit_comix$estimate[["sigma"]], " xi =", fit_comix$estimate[["xi"]], "\n")
cat("AIC =", fit_comix$aic, "\n")

cat("\n=== GPD fit: POLYMOD ===\n")
excess_polymod <- contact_data %>%
    filter(period == "POLYMOD") %>%
    filter(e_other >= 250) %>%
    pull(e_other) - 250
cat("n exceedances:", length(excess_polymod), "\n")
if (length(excess_polymod) >= 5) {
    fit_polymod <- fitdistrplus::fitdist(
        excess_polymod, "gpd",
        fix.arg = list(mu = 0),
        start = list(sigma = mean(excess_polymod), xi = 0.5),
        lower = c(1e-6, -Inf)
    )
    cat("sigma =", fit_polymod$estimate[["sigma"]], " xi =", fit_polymod$estimate[["xi"]], "\n")
    cat("AIC =", fit_polymod$aic, "\n")
} else {
    cat("Too few exceedances to fit GPD\n")
}

cat("\n=== GPD fit: POLYMOD + CoMix combined ===\n")
excess_combined <- contact_data %>%
    filter(period %in% c("POLYMOD", "Relaxed restrictions", "School reopening", "Step 2 + schools")) %>%
    filter(e_other >= 250) %>%
    pull(e_other) - 250
cat("n exceedances:", length(excess_combined), "\n")
fit_combined <- fitdistrplus::fitdist(
    excess_combined, "gpd",
    fix.arg = list(mu = 0),
    start = list(sigma = mean(excess_combined), xi = 0.5),
    lower = c(1e-6, -Inf)
)
cat("sigma =", fit_combined$estimate[["sigma"]], " xi =", fit_combined$estimate[["xi"]], "\n")
cat("AIC =", fit_combined$aic, "\n")
