source("scripts/utils.R")

contacts <- qread("data/contacts.qs")

dates_of_time_periods_of_interest <- 
  crossing(time_periods) %>% 
  filter(date_end<as.Date("2021-01-01"),period!="POLYMOD")

contacts_duration <- contacts[
  country=="uk" & 
    date >= min(dates_of_time_periods_of_interest$date_start) &
    date <= max(dates_of_time_periods_of_interest$date_end),
  .(cnt_minutes_max, cnt_household)
]
qsave(contacts_duration, "data/contacts_duration.qs")
