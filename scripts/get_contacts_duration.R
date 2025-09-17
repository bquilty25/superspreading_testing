source("scripts/utils.R")

contacts <- qread("data/contacts.qs")

contacts_duration <- contacts[,.(cnt_minutes_max, cnt_household)]
qsave(contacts_duration, "data/contacts_duration.qs")
