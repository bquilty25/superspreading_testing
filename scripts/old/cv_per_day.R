


traj_ %>% filter.(variant=="wild",lower_inf_thresh==F,heterogen_vl==T) %>% summarise.(sum_inf=sum(culture_p),.by=sim) %>% summarise.(cv=cv(sum_inf))
