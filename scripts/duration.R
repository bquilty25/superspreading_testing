source("scripts/utils.r")

contacts_duration <- qs::qread(file = "data/contacts_duration.qs")

contacts_hh_duration <- contacts_duration %>%
  drop_na.(cnt_minutes_max) %>%
  filter.(cnt_household == 1) %>%
  mutate.(
    cnt_duration = cnt_minutes_max / 1440,
    cnt_duration = ifelse(cnt_duration >= 1, 1, cnt_duration),
    period = ifelse(period %in% c("POLYMOD", "Pre-pandemic"), "Pre-pandemic", "Pandemic")
  ) %>%
  select.(part_id, cnt_duration, period)

contacts_nhh_duration <- contacts_duration %>%
  drop_na.(cnt_minutes_max) %>%
  filter.(cnt_household == 0) %>%
  mutate.(
    cnt_duration = cnt_minutes_max / 1440,
    cnt_duration = ifelse(cnt_duration >= 1, 1, cnt_duration),
    period = ifelse(period %in% c("POLYMOD", "Pre-pandemic"), "Pre-pandemic", "Pandemic")
  ) %>%
  select.(part_id, cnt_duration, period)

median_contact_duration <- median(c(contacts_hh_duration$cnt_duration, contacts_nhh_duration$cnt_duration), na.rm = TRUE)

contacts_duration %>%
  group_by(cnt_household) %>%
  summarise(
    q25 = quantile(cnt_minutes_max, 0.25, na.rm = T),
    q50 = quantile(cnt_minutes_max, 0.50, na.rm = T),
    q75 = quantile(cnt_minutes_max, 0.75, na.rm = T)
  )

contacts_duration %>%
  drop_na(cnt_minutes_max, cnt_household) %>%
  mutate(cnt_household = factor(cnt_household, levels = c("1", "0"))) %>%
  mutate(cnt_hours = pmin(cnt_minutes_max / 60, 24)) %>%
  ggplot() +
  geom_histogram(aes(x = cnt_hours, y = ..density.., fill = factor(cnt_household)), binwidth = 1) +
  scale_x_continuous("Per-contact time (hours)", breaks = scales::breaks_width(2)) +
  scale_y_continuous("Density") +
  MetBrewer::scale_fill_met_d(name = "Signac", override.order = FALSE, direction = -1, guide = F) +
  facet_wrap(~cnt_household, labeller = labeller(cnt_household = c(`1` = "Household", `0` = "Out of household"))) +
  plotting_theme +
  theme(legend.position = "none")

ggsave("results/manuscript_figures/fig_duration_hist.png", width = 210, height = 100, dpi = 600, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_duration_hist.pdf", width = 210, height = 100, dpi = 600, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_duration_hist.eps", width = 210, height = 100, units = "mm", device = cairo_ps)
