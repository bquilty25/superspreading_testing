# Rendered Number Check

This note summarises numeric differences between the legacy manuscript source in [docs/manuscript.md](docs/manuscript.md) and the latest rendered markdown produced from [docs/manuscript.qmd](docs/manuscript.qmd).

## Summary

Most differences appear to be intentional consequences of replacing approximate or stale manuscript values with inline R expressions backed by saved result objects. The main exceptions are a wording/statistic-label change in the contact-duration paragraph and a few testing-threshold statements that now differ materially from the legacy markdown.

## Notable Differences

### Contact summaries

- Legacy markdown: pandemic mean daily contacts `~6`; rendered markdown: `4.9`.
- Legacy markdown: pre-pandemic mean daily contacts `~12`; rendered markdown: `11.5`.
- Legacy markdown: pandemic average `k` `~0.6`; rendered markdown: `0.74`.
- Legacy markdown: mean daily contacts ranged from `~3` to `~7`; rendered markdown: `2.6` to `6.9`.
- Legacy markdown: post-first-lockdown contact `k` `~0.6`; rendered markdown: `0.48` to `0.68`.

Interpretation: these changes look intentional and data-driven.

### Contact-duration paragraph

- Legacy markdown: household contact duration reported as `480 minutes (8 hours) (95% CI: 180, 1080 minutes)`.
- Rendered markdown: household contact duration reported as `480 minutes (8 hours) (IQR: 180, 1080 minutes)`.

Interpretation: this is not just a numeric rounding change. The statistic label changed from `95% CI` to `IQR` and should be confirmed editorially.

### Viral-load / infectivity summaries

- Legacy markdown: `>70-fold`; rendered markdown: `73-fold`.
- Legacy markdown: peak culture probability `>0.8 for ~52%`; rendered markdown: `53%`.
- Legacy markdown: peak culture probability `>0.6 for ~68%`; rendered markdown: `68%`.

Interpretation: these changes look like expected rounding/precision updates from inline R.

### Secondary-infection estimates over time

- Legacy markdown: `R = 0.7` during first-lockdown easing; rendered markdown: `0.66`.
- Legacy markdown: `R = 1.0` during relaxed restrictions; rendered markdown: `1.02`.
- Legacy markdown: `R = 1.5` at school reopening; rendered markdown: `1.46`.
- Legacy markdown: `R = 1.0` during second lockdown; rendered markdown: `1.03`.
- Legacy markdown: heavier-tail sensitivity `R0 = 2.7 (95% CI 2.52–2.80)`; rendered markdown: `2.73 (95% CI 2.57–2.93)`.
- Legacy markdown: heavier-tail sensitivity `k = 0.46 (95% CI 0.43–0.49)`; rendered markdown: `0.46 (95% CI 0.42–0.5)`.

Interpretation: these changes look intentional and reflect more precise stored estimates.

### Heterogeneity paragraph

- Legacy markdown: Poisson contacts with variable viral load `k = 1.5-1.9`; rendered markdown: `1.46 to 1.88`.
- Legacy markdown: overdispersed contacts with variable viral load `k = 0.12-0.48`; rendered markdown: `0.11 to 0.51`.
- Legacy markdown: pre-pandemic overdispersed-contact comparison `0.74 vs. 0.48`; rendered markdown: `0.74 versus 0.51`.
- Legacy markdown: amplified viral-load sensitivity Poisson-contact range `0.75-0.89`; rendered markdown: `0.76 to 0.92`.
- Legacy markdown: amplified viral-load sensitivity overdispersed-contact range `0.08-0.32`; rendered markdown: `0.08 to 0.33`.

Interpretation: these changes look intentional and are consistent with exported result summaries.

### Rapid-testing thresholds

- Legacy markdown: every-3-day regular testing required `above 75%` uptake pre-pandemic; rendered markdown: `above 70%`.
- Legacy markdown: pre-event testing required `>85%` uptake for events with more than 10 others pre-pandemic; rendered markdown: `>80%`.
- Legacy markdown: school-reopening pre-event testing required `50%` adherence for events with more than 20 others; rendered markdown: `60%`.

Interpretation: these are the clearest substantive numeric changes. They appear to reflect updated values from the stored simulation outputs rather than the legacy prose.

### Repeated Discussion / conclusion summaries

- Legacy markdown: contact rates changed from `~12` to `~6` per day; rendered markdown: `11.5` to `4.9`.
- Legacy markdown: peak culture probability summaries `~68%` and `~52%`; rendered markdown: `68%` and `53%`.
- Legacy markdown: infectious period `approximately 2 days`; rendered markdown: still `2 days`.

Interpretation: repeated summary text now mirrors the inline-R-backed Results values.

## Recommended Follow-up

- Confirm whether the household-duration statement should report `IQR` or `95% CI`.
- If desired, align the legacy manuscript wording with the new data-backed testing thresholds to avoid future confusion.

## Audit Table

This table compares changed numeric statements between the legacy manuscript text in [docs/manuscript.md](docs/manuscript.md) and the latest rendered markdown from [docs/manuscript.qmd](docs/manuscript.qmd). The `Source` column refers to the file or computation now backing the rendered value.

| Section | Statement | Old value in `manuscript.md` | New rendered value | Source |
|---|---|---:|---:|---|
| Contacts | Pandemic mean daily contacts | `~6` | `4.9` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Pre-pandemic mean daily contacts | `~12` | `11.5` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Pandemic average `k` for contacts | `~0.6` | `0.74` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Range of daily contacts across 2020 periods | `~3` to `~7` | `2.6` to `6.9` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Post-first-lockdown contact `k` | `~0.6` | `0.48` to `0.68` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contact duration | Household duration interval label | `95% CI: 180, 1080 minutes` | `IQR: 180, 1080 minutes` | `data/contacts_duration.qs`, summarised in `scripts/duration.R` |
| Viral load | AUC fold difference | `>70-fold` | `73-fold` | `scripts/curve_plot.R` (`auc_dat`, `q95`) |
| Viral load | AUC 2.5th and 97.5th percentiles | `0.08` and `5.74` | `0.08` and `5.74` | `scripts/curve_plot.R` (`auc_dat`, `q95`) |
| Viral load | Gamma shape parameter | `1.42` | `1.42` | `scripts/curve_plot.R` (`auc_gamma_params`) |
| Viral load | Proportion infectious for zero days | `19.1%` | `19.1%` | `scripts/curve_plot.R` (`days_inf`, floored day count) |
| Viral load | Peak culture probability `>0.6` | `~68%` | `68%` | `scripts/curve_plot.R` threshold summary |
| Viral load | Peak culture probability `>0.8` | `~52%` | `53%` | `scripts/curve_plot.R` threshold summary |
| Secondary infections | Pre-pandemic `k` summary sentence | `0.48 (95% CI 0.47-0.51)` | `0.47 to 0.51 around a median of 0.48` | `results/R_and_k_bootstrap_ests.csv` |
| Secondary infections | First-lockdown-easing `R` | `0.7` | `0.66` | `results/R_and_k_bootstrap_ests.csv` |
| Secondary infections | Relaxed-restrictions `R` | `1.0` | `1.02` | `results/R_and_k_bootstrap_ests.csv` |
| Secondary infections | School-reopening `R` | `1.5` | `1.46` | `results/R_and_k_bootstrap_ests.csv` |
| Secondary infections | Second-lockdown `R` | `1.0` | `1.03` | `results/R_and_k_bootstrap_ests.csv` |
| Secondary infections | Heavier-tail sensitivity `R0` | `2.7 (95% CI 2.52–2.80)` | `2.73 (95% CI 2.57–2.93)` | `results/manuscript_tables/table_s1_gpd_sensitivity.csv` |
| Secondary infections | Heavier-tail sensitivity `k` CI | `0.46 (95% CI 0.43–0.49)` | `0.46 (95% CI 0.42–0.5)` | `results/manuscript_tables/table_s1_gpd_sensitivity.csv` |
| Heterogeneity | Poisson contacts, variable VL `k` | `1.5-1.9` | `1.46 to 1.88` | `results/R_and_k_bootstrap_ests_heterogen.csv` |
| Heterogeneity | Overdispersed contacts, equal VL `k` | `0.11-0.74` | `0.11 to 0.74` | `results/R_and_k_bootstrap_ests_heterogen.csv` |
| Heterogeneity | Overdispersed contacts, variable VL `k` | `0.12-0.48` | `0.11 to 0.51` | `results/R_and_k_bootstrap_ests_heterogen.csv` |
| Heterogeneity | Pre-pandemic equal vs variable VL comparison | `0.74 vs. 0.48` | `0.74 versus 0.51` | `results/R_and_k_bootstrap_ests_heterogen.csv` |
| Heterogeneity | Amplified VL sensitivity, Poisson contacts | `0.75-0.89` | `0.76 to 0.92` | `results/k_bootstrap_ests_vl_sens.csv` |
| Heterogeneity | Amplified VL sensitivity, overdispersed contacts | `0.08-0.32` | `0.08 to 0.33` | `results/k_bootstrap_ests_vl_sens.csv` |
| Testing | Pre-pandemic, daily regular testing threshold | `>60%` | `>60%` | `results/processed_infections_testing.qs` via `scripts/results.R` |
| Testing | Pre-pandemic, every-3-day regular testing threshold | `>75%` | `>70%` | `results/processed_infections_testing.qs` via `scripts/results.R` |
| Testing | School-reopening, every-3-day regular testing threshold | `>60%` | `>60%` | `results/processed_infections_testing.qs` via `scripts/results.R` |
| Testing | Pre-pandemic, pre-event testing \>10 threshold | `>85%` | `>80%` | `results/processed_infections_events.qs` via `scripts/results.R` |
| Testing | School-reopening, pre-event testing \>20 threshold | `50%` | `60%` | `results/processed_infections_events.qs` via `scripts/results.R` |
| Discussion repeat | Contact-rate summary | `~12` to `~6` | `11.5` to `4.9` | Same as contact summary above |
| Discussion repeat | High-infectiousness proportions | `~68%`, `~52%` | `68%`, `53%` | Same as viral-load summary above |
| Conclusion repeat | Infectious period | `approximately 2 days` | `2 days` | `scripts/curve_plot.R` (`days_inf`) |

## Audit Notes

- Some rows show no substantive numeric change, only a precision or presentation change. They are included because the statement was re-rendered from a different source path.
- The most substantive numeric differences are in the contact-summary paragraph, the heterogeneity paragraph, and the rapid-testing thresholds.
- The contact-duration row is included because the statistic label changed, even though the numbers themselves did not.