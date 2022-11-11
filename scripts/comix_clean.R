#Load data
comix_contact_common <-
  read_csv("data/comix/CoMix_uk_contact_common.csv")

comix_contact_extra <-
  read_csv("data/comix/CoMix_uk_contact_extra.csv")

comix_hh_common <- 
  read_csv("data/comix/CoMix_uk_hh_common.csv")

comix_participant_common <-
  read_csv("data/comix/CoMix_uk_participant_common.csv")

comix_participant_extra <-
  read_csv("data/comix/CoMix_uk_participant_extra.csv")

comix_sday <- 
  read_csv("data/comix/CoMix_uk_sday.csv")

comix <- list(
  participants = comix_participant_common %>%
    left_join(comix_participant_extra) %>%
    left_join(comix_sday),
  contacts = comix_contact_common %>%
    left_join(comix_contact_extra) %>%
    left_join(comix_sday)
)

#Multiple contacts
nhh_cols_to_sum <- comix$contacts %>%
  select(contains("cnt_")) %>%
  select(cnt_work:cnt_shop) %>%
  colnames()

contacts_mult <- comix$contacts %>%
  mutate(e_other = cnt_home == F,
         e_home = cnt_home == T) %>%
  pivot_longer(c(e_home, e_other)) %>%
  filter(value == T) %>%
  count(part_id, day, month, year, name, .drop = F) %>%
  pivot_wider(values_from = n, names_from = name) %>%
  replace_na(list(e_home = 0, e_other = 0))

#Zero contacts
contacts_zero <- comix$participants %>%
  anti_join(comix$contacts) %>%
  select(part_id, day, month, year) %>%
  mutate(e_home = 0,
         e_other = 0)

contacts_comix <- bind_rows(contacts_mult,
                          contacts_zero) %>%
  arrange(part_id, day, month, year) %>%
  unite("date", c(day, month, year), sep = "-") %>%
  mutate(date = as.Date(date, format = "%d-%m-%y")) %>%
  mutate(e_all = e_home + e_other)
