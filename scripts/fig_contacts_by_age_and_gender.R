# Plot UK ONS mid-2020, BBC Pandemic and CoMix age distributions
source("scripts/utils.R")

uk_age_dist <- read_csv("data/ons_uk_age_distn_2020.csv")

age_grps <- c("0-4", "5-11", "12-17", "18-29", "30-39", "40-49", "50-59", 
              "60-69", "70-120")
uk_age <- uk_age_dist %>%
  mutate(part_age = cut(age, 
                         breaks = c(as.numeric(sub("-.*","",age_grps)), Inf), 
                         labels = age_grps,
                         right = F)) %>%
  group_by(part_age) %>%
  summarise(pop = sum(population), .groups = "drop_last") %>%
  mutate(
    period = "UK mid-2020",
    proportion = pop/sum(pop)
  ) %>%
  select(period, part_age, proportion)

bbc_age_dist <- read_csv("data/bbc_age_distn.csv")

bbc_part_age <- bbc_age_dist %>%
  mutate(
    period = "Pre-pandemic",
    part_age = case_when(
      age_group %in% c("0-4") ~ "0-4",
      age_group %in% c("4-9", "10-12") ~ "5-11",
      age_group %in% c("12-17", "13-14", "15-17") ~ "12-17",
      age_group %in% c("18-19", "20-21", "22-24", "25-29") ~ "18-29",
      age_group %in% c("30-34", "35-39") ~ "30-39",
      age_group %in% c("40-44", "45-49") ~ "40-49",
      age_group %in% c("50-54", "55-59") ~ "50-59",
      age_group %in% c("60-64", "65-69") ~ "60-69",
      age_group %in% c("70-74", "75+") ~ "70-120"
    )
  ) %>%
  group_by(period, part_age) %>%
  summarise(
    proportion = sum(proportion),
    .groups = "drop"
  ) %>%
  mutate(
    part_age = factor(
      part_age,
      levels = age_grps
    )
  ) %>%
  arrange(part_age)

comix_part_age <- contact_data %>%
  drop_na(part_age) %>%
  group_by(period, part_age) %>%
  summarise(count = n(), .groups = "drop_last") %>% 
  mutate(proportion = count / sum(count))

part_age <- bind_rows(uk_age, bbc_part_age, comix_part_age) %>%
  mutate(period = factor(
    period,
    levels = c("UK mid-2020", "Pre-pandemic", "1st lockdown", "1st lockdown easing",
      "Relaxed restrictions", "School reopening", "2nd lockdown", 
      "2nd lockdown easing", "3rd lockdown", "3rd lockdown + schools", 
      "Step 2 + schools")
    )
  )

ggplot(part_age) + 
  geom_col(aes(x = period, y = proportion, fill = part_age), 
           position = position_stack(reverse = T)) +
  labs(x = "Time period", y = "Proportion") + 
  scale_fill_discrete(name = "Age group") + 
  plotting_theme +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("results/manuscript_figures/fig_age_distribution.png", dpi = 600, width = 200, height = 133, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_age_distribution.pdf", dpi = 600, width = 200, height = 133, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_age_distribution.eps", width = 200, height = 133, units = "mm", device = cairo_ps)

# Plot CoMix contact distributions by age and gender and time period
cnt_by_age <- contact_data %>%
  drop_na(part_age) %>%
  pivot_longer.(cols = c(e_home, e_other, e_all)) %>%
  ggplot() +
  geom_boxplot(aes(x = part_age, y = value, fill = period), 
               position = "dodge") +
  facet_wrap2( ~ name,
              labeller = labeller(name = c(
                "e_all" = "All contacts",
                "e_home" = "Household contacts",
                "e_other" = "Out-of-household contacts"
              )),
              axes = "all"
  ) +
  scale_y_continuous(trans = "pseudo_log", 
                     breaks = c(0, 1, 10, 100, 1000), 
                     expand = expansion(0, 0)) + 
  plotting_theme +
  theme(axis.title.y = element_blank()) + 
  labs(x = "Age group")

cnt_by_gender <- contact_data %>%
  filter(part_gender != "other") %>%
  drop_na(part_gender) %>%
  pivot_longer.(cols = c(e_home, e_other, e_all)) %>%
  ggplot() +
  geom_boxplot(aes(x = part_gender, y = value, fill = period), 
               position = "dodge") +
  facet_wrap2( ~ name, 
              axes = "all"
  ) +
  scale_y_continuous(trans = "pseudo_log", 
                     breaks = c(0, 1, 10, 100, 1000), 
                     expand = expansion(0, 0)) + 
  plotting_theme + 
  theme(axis.title.y = element_blank(), strip.text.x = element_blank()) + 
  scale_fill_discrete(name = "Period") +
  labs(x = "Gender")

(cnt_by_age + theme(legend.position = "none")) /
  cnt_by_gender +
  labs(tag = "Reported daily contacts") + 
  theme(plot.tag = element_text(angle = 90, vjust = 3, hjust = 1.5),
        plot.tag.position = "left")

ggsave("results/manuscript_figures/fig_contacts_by_age_and_gender.png", dpi = 600, width = 350, height = 200, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_contacts_by_age_and_gender.pdf", dpi = 600, width = 350, height = 200, units = "mm", bg = "white")
ggsave("results/manuscript_figures/fig_contacts_by_age_and_gender.eps", width = 350, height = 200, units = "mm", device = cairo_ps)
