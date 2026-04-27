source("scripts/utils.R")

contacts <- qread("data/contacts.qs")

dates_of_time_periods_of_interest <-
  crossing(time_periods) %>%
  filter(period!="POLYMOD")
#   filter(date_end<as.Date("2021-01-01"),period!="POLYMOD")

contacts_duration <- contacts[
  country=="uk" & 
    date >= min(dates_of_time_periods_of_interest$date_start) &
    date <= max(dates_of_time_periods_of_interest$date_end),
  .(part_id = gsub("uk_", "", part_wave_uid), date, cnt_minutes_max, cnt_household, cnt_total_time)
]

contacts_polymod1 <- fread("data/POLYMOD/2008_Mossong_POLYMOD_contact_common.csv")
survey_dates <- fread("data/POLYMOD/2008_Mossong_POLYMOD_sday.csv")
survey_dates[, date := as.Date(as.character(sday_id), "%Y%m%d")]
contacts_polymod1 <- merge(contacts_polymod1, survey_dates, by = "part_id")
parts_polymod <- fread("data/POLYMOD/2008_Mossong_POLYMOD_participant_common.csv")
hh_polymod <- fread("data/POLYMOD/2008_Mossong_POLYMOD_hh_common.csv")
parts_polymod <- merge(parts_polymod, hh_polymod, by = "hh_id")
contacts_polymod1 <- merge(contacts_polymod1, parts_polymod, by = "part_id")

contacts_duration_polymod <- contacts_polymod1[
  country == "GB",
  .(part_id,
     date,
     cnt_household = as.numeric(cnt_home),
     cnt_total_time = fcase(
       duration_multi == 1, "<5m",
       duration_multi == 2, "5m-14m",
       duration_multi == 3, "15m-59m",
       duration_multi == 4, "60m-4h",
       duration_multi == 5, "4h+"
     )
  )
]

contacts_duration_polymod[, date := nafill(date, type = "locf")]
  
contacts_duration <- rbind(contacts_duration_polymod, contacts_duration, fill = T)

# ggplot(contacts_duration[!is.na(cnt_total_time)]) + geom_boxplot(aes(x=cnt_total_time,y=cnt_minutes_max))

# Impute missing exact contact durations for contacts with only range for 
# duration from observed data
set.seed(1)

contacts_duration[,
  cnt_minutes_max := if (all(is.na(cnt_minutes_max))) {
    cnt_minutes_max
  } else {
    na_idx <- is.na(cnt_minutes_max)
    vals <- cnt_minutes_max[!na_idx]
    cnt_minutes_max[na_idx] <- sample(vals, sum(na_idx), replace = TRUE)
    cnt_minutes_max
  },
  by = cnt_total_time
]

contacts_duration <- contacts_duration %>%
  fuzzyjoin::fuzzy_left_join(time_periods,
                              by = c("date" = "date_start", "date" = "date_end"),
                              match_fun = list(`>=`, `<=`)
  )

setDT(contacts_duration)

# ggplot(contacts_duration[!is.na(period)]) + geom_boxplot(aes(x = period, y = cnt_minutes_max))

qsave(contacts_duration, "data/contacts_duration.qs")

# for (prd in unique(contacts_duration[!is.na(period),period])){
#   print(ggplot(contacts_duration[!is.na(cnt_total_time) & period == prd], aes(x = cnt_total_time)) + geom_bar())
# }
