library(qs)
library(data.table)

# Participant info
pt <- qread("~/Filr/Net Folders/EPH Shared/comix_shared/private_do_not_share/part_min.qs")

# Contact count data
cnt <- qread("../comix_analysis/data/part_cnts.qs")

# Merge with participant info
pt_cnt <- merge(pt,cnt,by="part_wave_uid",all=T)

# Columns to include
# cols <- c("panel","wave","date",grep("n_cnt",names(pt_cnt),value = T))
cols <- c("panel","wave","date","part_age_group","n_cnt_home","n_cnt_school","n_cnt_work","n_cnt_other")
pt_cnt <- pt_cnt[country=="uk",..cols]
cnt_cols <- grep("n_cnt",names(pt_cnt),value = T)
setnames(pt_cnt,cnt_cols,sub("n_cnt","e",cnt_cols))

# Split data into under-18s and over-18s
# Under-18s
pt_cnt_u18 <- pt_cnt[part_age_group %in% c("0-4","5-11","12-17")]
pt_cnt_u18[,part_age_group:=NULL]
# Over-18s
pt_cnt_o18 <- pt_cnt[!(part_age_group %in% c(NA,"0-4","5-11","12-17")),]
pt_cnt_o18[,part_age_group:=NULL]

# Save
write.csv(pt_cnt_u18,"data/comix_contact_distributions_u18.csv",row.names = F)
write.csv(pt_cnt_o18,"data/comix_contact_distributions_o18.csv",row.names = F)
