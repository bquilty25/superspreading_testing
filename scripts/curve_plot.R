## curves for paper

# Load required packages and utility scripts
source("packages.R")
source("utils.R")
source("plot_functions.R")
source("parameters.R")
source("lfa_test_data.R")


set.seed(2020)

traj <- readRDS(file.info(list.files("results/", full.names = T)) %>% 
                  as.data.frame() %>% 
                  rownames_to_column()%>% 
                  filter(str_detect(rowname,"traj_")) %>% 
                  slice_max(mtime) %>% 
                  pull(rowname))

trajectories_to_plot <- traj %>%
  #filter(type=="symptomatic") %>% 
  group_by(type) %>% 
  sample_n(1000) %>% 
  #ungroup() %>% 
  mutate(pred = map(.x = m, 
                    ~data.frame(x = seq(0, 25, length.out = 26)) %>%
                      mutate(ct = .x(x)))) %>%
  unnest(pred) %>%
  mutate(LFT = stats::predict(innova_mod,
                              newdata =
                                data.frame(ct = ct),
                              type = "response"
  ),
  culture=stats::predict(culture_mod,
                         newdata =
                           data.frame(ct = ct),
                         type = "response"
  ),
  PCR=ifelse(ct<35,1,0)) %>% 
  arrange(type)

trajectories_to_plot %>% 
  #pivot_longer(cols=c(ct,culture,LFT)) %>%  
  ggplot()+
  ggdist::stat_lineribbon(aes(x=factor(x),y=10^(12 - 0.328*ct)))+
  scale_fill_brewer()+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  #scale_color_manual(values="black",name="",guide=F)+
  #scale_x_continuous(name="Days since exposure",breaks = breaks_width(7))
  # scale_y_reverse(# Features of the first axis
  #   limits=c(39.99,NA),
  #   name = "Cycle\nthreshold")+
  plotting_theme+
  facet_wrap(~fct_rev(type))
