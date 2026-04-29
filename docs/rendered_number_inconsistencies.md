# Rendered Number Check

This note summarises numeric differences between the legacy manuscript source in [docs/manuscript.md](docs/manuscript.md) and the latest rendered markdown produced from [docs/manuscript.qmd](docs/manuscript.qmd).

## Summary

Most differences now appear to be direct consequences of deriving manuscript numbers from the current contact-data objects and saved model outputs rather than keeping the older hand-entered prose values. Compared with the previous audit pass, the biggest newly visible changes are in the contact table and contact-summary paragraph, plus the rapid-testing thresholds. There is also still a statistic-label change in the contact-duration paragraph.

## Assessment

### Likely Improvements

- The secondary-infection `R` and `k` values now match the saved bootstrap summaries and should be treated as the better source of truth than the older rounded markdown prose.
- The heterogeneity paragraph now matches the exported summary files for the main and amplified viral-load sensitivity analyses.
- The viral-load summary is mostly a precision update, with the `73-fold` AUC range and `68%`/`53%` peak-culture proportions consistent with the current model-backed setup.
- The repeated Discussion and conclusion values now track the same inline-R-backed quantities used in Results, which is structurally better than repeating stale prose numbers.

### Needs Confirmation

- Table 1 is now generated from the live `contact_data` object and differs materially from the legacy markdown table, especially in response counts; this needs a decision on whether the old table or the current derived table is canonical.
- The contact-summary paragraph now uses means and dispersion derived from the same live contact object that drives Table 1, so if Table 1 is disputed, these downstream prose values are disputed too.
- The contact-duration paragraph now reports `IQR` instead of `95% CI`; the rendered value is internally coherent, but the manuscript should use one statistic label deliberately.
- The rapid-testing thresholds now differ enough from the legacy prose that they should be explicitly accepted as the current model outputs before finalising the narrative.

## Notable Differences

### Contact summaries

- Legacy markdown: pandemic mean daily contacts `~6`; rendered markdown: `4.7`.
- Legacy markdown: pre-pandemic mean daily contacts `~12`; rendered markdown: `11.5`.
- Legacy markdown: pandemic average `k` `~0.6`; rendered markdown: `0.69`.
- Legacy markdown: mean daily contacts ranged from `~3` to `~7`; rendered markdown: `2.6` to `6.9`.
- Legacy markdown: post-first-lockdown contact `k` `~0.6`; rendered markdown: `0.48` to `0.68`.

Interpretation: these changes look intentional and data-driven.

### Contact table values

- The rendered Table 1 now differs from the legacy markdown table in both `N` and several percentage columns.
- Example rows:
	- 1st lockdown `N`: legacy `15906`; rendered `15971`.
	- School reopening `N`: legacy `20759`; rendered `21031`.
	- Step 2 + schools `N`: legacy `1771`; rendered `11607`.
	- Relaxed restrictions `Over 20 (%)`: legacy `3.2`; rendered `3.1`.
	- School reopening `Over 10 (%)`: legacy `9.0`; rendered `8.7`.

Interpretation: the table is now being generated from the current `contact_data` object in the manuscript setup, not preserved from the older static markdown table. This is a substantive source-path change, not just rounding.

### Table 1 row checks

- Pre-pandemic row is effectively unchanged.
- 1st lockdown is close to the legacy table, but not identical: `N` increased from `15906` to `15971`, and `Over 20 (%)` fell from `0.5` to `0.4`.
- School reopening is close but systematically lower in some columns: `N` increased from `20759` to `21031`, `Over 10 (%)` fell from `9.0` to `8.7`, and `Over 50 (%)` fell from `2.1` to `2.0`.
- Step 2 + schools is the clearest mismatch: `N` increased from `1771` to `11607`, with several percentages shifting accordingly.

Interpretation: this looks like a genuine dataset-definition or filtering difference, not a formatting issue.

### Contact-duration paragraph

- Legacy markdown: household contact duration reported as `480 minutes (8 hours) (95% CI: 180, 1080 minutes)`.
- Rendered markdown: household contact duration reported as `480 minutes (8 hours) (IQR: 180, 1080 minutes)`.

Interpretation: this is not just a numeric rounding change. The statistic label changed from `95% CI` to `IQR` and should be confirmed editorially.

### Viral-load / infectivity summaries

- Legacy markdown: `>70-fold`; rendered markdown: `73-fold`.
- Legacy markdown: zero infectious days `19.1%`; rendered markdown: `19.4%`.
- Legacy markdown: peak culture probability `>0.8 for ~52%`; rendered markdown: `53%`.
- Legacy markdown: peak culture probability `>0.6 for ~68%`; rendered markdown: `68%`.

Interpretation: most of these changes look like expected rounding/precision updates from model-backed inline R. The zero-days value is a small but real numerical shift relative to the old markdown.

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

- Legacy markdown: daily regular testing required `exceed 60%` uptake pre-pandemic; rendered markdown: `exceed 70%`.
- Legacy markdown: every-3-day regular testing required `above 75%` uptake pre-pandemic; rendered markdown: `above 80%`.
- Legacy markdown: school-reopening every-3-day regular testing required `exceeds 60%`; rendered markdown: `exceeds 70%`.
- Legacy markdown: pre-event testing required `>85%` uptake for events with more than 10 others pre-pandemic; rendered markdown: `>90%`.
- Legacy markdown: school-reopening pre-event testing required `50%` adherence for events with more than 20 others; rendered markdown: `60%`.

Interpretation: these are now the clearest substantive prose changes. They reflect the current saved simulation outputs rather than the legacy rounded narrative values.

### Repeated Discussion / conclusion summaries

- Legacy markdown: contact rates changed from `~12` to `~6` per day; rendered markdown: `11.5` to `4.7`.
- Legacy markdown: peak culture probability summaries `~68%` and `~52%`; rendered markdown: `68%` and `53%`.
- Legacy markdown: infectious period `approximately 2 days`; rendered markdown: still `2 days`.

Interpretation: repeated summary text now mirrors the inline-R-backed Results values.

## Recommended Follow-up

- Confirm whether Table 1 should remain dynamically generated from `contact_data`, because it now differs materially from the older markdown table, especially in response counts.
- Confirm whether the household-duration statement should report `IQR` or `95% CI`.
- Align the narrative testing-threshold wording with the current model-backed thresholds if the rendered manuscript is now the canonical source.

## Manuscript Pass

If the current model outputs are intended to be canonical, the following prose now looks defensible as written:

- The Results paragraph on secondary-infection estimates over time.
- The heterogeneity paragraph in Results.
- The repeated Discussion/conclusion summaries for `R`, contact means, and viral-load-derived infectiousness proportions.

The following items should still be reviewed before treating the manuscript as final:

- Table 1 and all contact-summary prose tied to it.
- The contact-duration statistic label.
- The rapid-testing threshold sentences in Results and Discussion.

## Audit Table

This table compares changed numeric statements between the legacy manuscript text in [docs/manuscript.md](docs/manuscript.md) and the latest rendered markdown from [docs/manuscript.qmd](docs/manuscript.qmd). The `Source` column refers to the file or computation now backing the rendered value.

| Section | Statement | Old value in `manuscript.md` | New rendered value | Source |
|---|---|---:|---:|---|
| Contacts | Pandemic mean daily contacts | `~6` | `4.7` | Dynamic `contact_summary` in [docs/manuscript.qmd](docs/manuscript.qmd) from `contact_data` |
| Contacts | Pre-pandemic mean daily contacts | `~12` | `11.5` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Pandemic average `k` for contacts | `~0.6` | `0.69` | Dynamic `contact_summary` in [docs/manuscript.qmd](docs/manuscript.qmd) from `contact_data` |
| Contacts | Range of daily contacts across 2020 periods | `~3` to `~7` | `2.6` to `6.9` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Contacts | Post-first-lockdown contact `k` | `~0.6` | `0.48` to `0.68` | Negative binomial fits in `scripts/contact_plots.R` from `contact_data` |
| Table 1 | 1st lockdown `N` | `15906` | `15971` | Dynamic `contact_table` in [docs/manuscript.qmd](docs/manuscript.qmd) from `contact_data` |
| Table 1 | School reopening `N` | `20759` | `21031` | Dynamic `contact_table` in [docs/manuscript.qmd](docs/manuscript.qmd) from `contact_data` |
| Table 1 | Step 2 + schools `N` | `1771` | `11607` | Dynamic `contact_table` in [docs/manuscript.qmd](docs/manuscript.qmd) from `contact_data` |
| Contact duration | Household duration interval label | `95% CI: 180, 1080 minutes` | `IQR: 180, 1080 minutes` | `data/contacts_duration.qs`, summarised in `scripts/duration.R` |
| Viral load | AUC fold difference | `>70-fold` | `73-fold` | `scripts/curve_plot.R` (`auc_dat`, `q95`) |
| Viral load | AUC 2.5th and 97.5th percentiles | `0.08` and `5.74` | `0.08` and `5.74` | `scripts/curve_plot.R` (`auc_dat`, `q95`) |
| Viral load | Gamma shape parameter | `1.42` | `1.42` | `scripts/curve_plot.R` (`auc_gamma_params`) |
| Viral load | Proportion infectious for zero days | `19.1%` | `19.4%` | Dynamic Figure 3 model summary in [docs/manuscript.qmd](docs/manuscript.qmd) |
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
| Testing | Pre-pandemic, daily regular testing threshold | `>60%` | `>70%` | `results/processed_infections_testing.qs` via manuscript setup summary |
| Testing | Pre-pandemic, every-3-day regular testing threshold | `>75%` | `>80%` | `results/processed_infections_testing.qs` via manuscript setup summary |
| Testing | School-reopening, every-3-day regular testing threshold | `>60%` | `>70%` | `results/processed_infections_testing.qs` via manuscript setup summary |
| Testing | Pre-pandemic, pre-event testing \>10 threshold | `>85%` | `>90%` | `results/processed_infections_events.qs` via manuscript setup summary |
| Testing | School-reopening, pre-event testing \>20 threshold | `50%` | `60%` | `results/processed_infections_events.qs` via `scripts/results.R` |
| Discussion repeat | Contact-rate summary | `~12` to `~6` | `11.5` to `4.7` | Same as contact summary above |
| Discussion repeat | High-infectiousness proportions | `~68%`, `~52%` | `68%`, `53%` | Same as viral-load summary above |
| Conclusion repeat | Infectious period | `approximately 2 days` | `2 days` | `scripts/curve_plot.R` (`days_inf`) |

## Audit Notes

- Some rows show no substantive numeric change, only a precision or presentation change. They are included because the statement was re-rendered from a different source path.
- The most substantive numeric differences are now in the contact-summary paragraph, Table 1 counts/percentages, and the rapid-testing thresholds.
- The contact-duration row is included because the statistic label changed, even though the numbers themselves did not.