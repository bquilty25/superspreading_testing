# Load required packages scripts
library(qs)
library(ggdist)
library(fitdistrplus)
library(EnvStats)
library(tidyverse)
library(patchwork)
library(here)
library(rriskDistributions)
library(rms)
library(DescTools)
library(MESS)
library(lubridate)
library(lemon)
library(boot)
library(furrr)
library(data.table)
library(ggtext)
library(fst)
library(extraDistr)
library(emdbook)
library(colorspace)
library(fuzzyjoin)
library(ggpubr)
library(bench)
library(tictoc)
library(naniar)
library(scales)
library(ggforce)
library(RGeode)
library(tsibble)
library(MetBrewer)
library(ggrepel)
library(ggh4x)
library(geomtextpath)
library(ggnewscale)

if (packageVersion("tidytable") != "0.8.0") {
  remotes::install_version("tidytable", version = "0.8")
} else {
  library(tidytable)
}

seed <- 1000
set.seed(seed)

# Transmission scaling factor: calibrated in main.R so pre-pandemic mean R = target_R0
beta_inf <- 1

# plotting options
covid_pal <- c("#e66101", "#5e3c99", "#0571b0")

`%!in%` <- Negate(`%in%`)

plotting_theme <- theme_minimal(
  # base_family = "Lato"
) +
  theme(
    axis.ticks = element_line(colour = "#2E4C6D"),
    axis.title = element_text(colour = "#2E4C6D"),
    axis.text = element_text(colour = "#2E4C6D"),
    strip.text = element_text(colour = "#2E4C6D"),
    axis.line.x = element_line(colour = "#2E4C6D"),
    axis.line.y = element_line(colour = "#2E4C6D"),
    # panel.border = element_rect(fill=NA,colour="#2E4C6D"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    strip.placement = "outside",
    axis.line = element_line(colour = "#2E4C6D"),
    line = element_line(colour = "#2E4C6D"),
    text = element_text(
      colour = "#2E4C6D"
      # ,family = "Lato"
    )
  )

bi_col_pal <- c("#396EB0", "#FC997C")
tri_col_pal <- c("#396EB0", "#DADDFC", "#FC997C")
quad_col_pal <- c("#2E4C6D", "#396EB0", "#DADDFC", "#FC997C")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

### Load viral load data ----
kissler_dat <- read_csv("CtTrajectories_B117/output/shared_params_df.csv") %>%
  select(
    "alpha_peakvl" = dpmeanB,
    "wild_peakvl" = dpmeanW,
    "alpha_prolif" = wpmeanB,
    "wild_prolif" = wpmeanW,
    "alpha_clear" = wrmeanB,
    "wild_clear" = wrmeanW,
    "peakvl_sd" = dpsd,
    "prolif_sd" = wpsd,
    "clear_sd" = wrsd
  )

kissler_dat_means <- kissler_dat %>%
  select(-c(peakvl_sd:clear_sd)) %>%
  mutate(across(.cols = contains("vl"), ~ 40 - .x)) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean = median(value)) %>%
  separate(name, sep = "_", into = c("variant", "param"))

kissler_dat_sd <- kissler_dat %>%
  select(c(peakvl_sd:clear_sd)) %>%
  pivot_longer(everything(), names_to = "param") %>%
  group_by(param) %>%
  summarise(sd = median(value)) %>%
  mutate(param = str_extract(param, "[^_]+"))

kissler_dat_est <- kissler_dat_means %>%
  left_join(kissler_dat_sd) %>%
  pivot_wider(names_from = param, values_from = c(mean, sd))

hay_dat <- read_csv(file = "data/shared_params_df.csv") %>%
  select(
    "omicron_peakvl" = dpmeanB_trans,
    "delta_peakvl" = dpmeanW_trans,
    "omicron_prolif" = wpmeanB_trans,
    "delta_prolif" = wpmeanW_trans,
    "omicron_clear" = wrmeanB_trans,
    "delta_clear" = wrmeanW_trans,
    "peakvl_sd" = dpsd,
    "prolif_sd" = wpsd,
    "clear_sd" = wrsd
  )

hay_dat_means <- hay_dat %>%
  select(-c(peakvl_sd:clear_sd)) %>%
  mutate(across(.cols = contains("vl"), ~ 40 - .x)) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean = median(value)) %>%
  separate(name, sep = "_", into = c("variant", "param"))

hay_dat_sd <- hay_dat %>%
  select(c(peakvl_sd:clear_sd)) %>%
  pivot_longer(everything(), names_to = "param") %>%
  group_by(param) %>%
  summarise(sd = median(value)) %>%
  mutate(param = str_extract(param, "[^_]+"))

hay_dat_est <- hay_dat_means %>%
  left_join(hay_dat_sd) %>%
  pivot_wider(names_from = param, values_from = c(mean, sd))

vl_params <- bind_rows(kissler_dat_est, hay_dat_est) %>%
  mutate(variant = as_factor(variant))

# From Kissler et al. https://github.com/gradlab/CtTrajectories/blob/main/code/utilities/utils_analysis.R
convert_Ct_logGEML <- function(Ct, m_conv = -3.609714286, b_conv = 40.93733333) {
  out <- (Ct - b_conv) / m_conv * log10(10) + log10(250)
  return(out)
}

# Define time periods of interest
time_periods <- tribble(
  ~idx, ~period, ~date_start, ~date_end,
  -1, "POLYMOD", as_date("12/05/2005", format = "%d/%m/%Y"), as_date("05/09/2006", format = "%d/%m/%Y"),
  0, "Pre-pandemic", as_date("01/09/2017", format = "%d/%m/%Y"), as_date("01/12/2018", format = "%d/%m/%Y"),
  1, "1st lockdown", as_date("23/03/2020", format = "%d/%m/%Y"), as_date("03/06/2020", format = "%d/%m/%Y"),
  2, "1st lockdown easing", as_date("04/06/2020", format = "%d/%m/%Y"), as_date("29/07/2020", format = "%d/%m/%Y"),
  3, "Relaxed restrictions", as_date("30/07/2020", format = "%d/%m/%Y"), as_date("03/09/2020", format = "%d/%m/%Y"),
  4, "School reopening", as_date("04/09/2020", format = "%d/%m/%Y"), as_date("24/10/2020", format = "%d/%m/%Y"),
  5, "2nd lockdown", as_date("05/11/2020", format = "%d/%m/%Y"), as_date("02/12/2020", format = "%d/%m/%Y"),
  6, "2nd lockdown easing", as_date("03/12/2020", format = "%d/%m/%Y"), as_date("19/12/2020", format = "%d/%m/%Y"),
  7, "3rd lockdown", as_date("05/01/2021", format = "%d/%m/%Y"), as_date("07/03/2021", format = "%d/%m/%Y"),
  8, "3rd lockdown + schools", as_date("08/03/2021", format = "%d/%m/%Y"), as_date("31/03/2021", format = "%d/%m/%Y"),
  9, "Step 2 + schools", as_date("16/04/2021", format = "%d/%m/%Y"), as_date("16/05/2021", format = "%d/%m/%Y")
) %>%
  mutate(period = factor(period, levels = period))

# Scenarios to investigate
key_grouping_var <- c("sim", "variant", "period", "lower_inf_thresh", "heterogen_vl", "heterogen_contacts")

#### Load contact data ----
contacts_polymod <-
  read.csv(here("data", "POLYMOD/2008_Mossong_POLYMOD_contact_common.csv")) %>%
  pivot_longer(cols = c(cnt_home, cnt_school, cnt_work, cnt_transport, cnt_leisure, cnt_otherplace)) %>%
  filter(value) %>%
  select(-value) %>%
  count(name, part_id) %>%
  pivot_wider(names_from = name, values_from = n) %>%
  mutate(
    e_home = cnt_home,
    e_other = rowSums(across(c(cnt_work, cnt_school, cnt_transport, cnt_leisure, cnt_otherplace)), na.rm = T)
  ) %>%
  select(part_id, e_home, e_other) %>%
  complete(part_id = full_seq(part_id, 1), fill = list(e_home = 0, e_other = 0)) %>%
  mutate(date = as_date("01/01/2008", format = "%d/%m/%Y"))

contacts_bbc_o18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_o18.csv"))

contacts_bbc_u18 <-
  read.csv(here("2020-cov-tracing", "data", "contact_distributions_u18.csv"))

contacts_bbc <- bind_rows(contacts_bbc_o18, contacts_bbc_u18) %>%
  mutate(
    date = as_date("01/09/2017", format = "%d/%m/%Y"),
    e_school = 0
  )

# Clean and load Comix data
source("scripts/comix_clean.R")

# summarise number of contacts
contact_data <- contacts_bbc %>%
  mutate(e_other = rowSums(across(c(e_work, e_school, e_other)), na.rm = T)) %>%
  select(date, e_home, e_other) %>%
  bind_rows(contacts_polymod) %>%
  mutate(
    e_all = rowSums(across(c(e_home, e_other)), na.rm = T),
    date = as_date(date),
    part_id = as.character(part_id)
  ) %>%
  bind_rows(contacts_comix) %>%
  fuzzyjoin::fuzzy_inner_join(time_periods,
    by = c("date" = "date_start", "date" = "date_end"),
    match_fun = list(`>=`, `<=`)
  ) %>%
  mutate(
    id = gsub("[A-Z][0-9]{1,2}_", "", part_id),
    idx_id = paste(idx, id, sep = "_")
  )

#### impute out of HH values > 250 for Pre-pandemic by fitting distribution to values from non-lockdown periods ----

# Calculate proportion over 250 by time period
contact_data %>%
  filter(period %in% c("Relaxed restrictions", "School reopening", "Step 2 + schools")) %>%
  summarise.(n = n(), over_250 = sum(e_other >= 250)) %>%
  mutate.(prop = over_250 / n)

# 0.00160 or 0.16% over 250

# assume distribution of high contacts is exponential and fit distribution
dist_over_250 <- contact_data %>%
  filter(period %in% c("Relaxed restrictions", "School reopening", "Step 2 + schools")) %>%
  filter(e_other >= 250) %>%
  pull(e_other) %>%
  fitdistr(., "exponential")

# simulate individuals with high numbers of contacts for Pre-pandemic
dat_append <- data.frame(
  e_other = round(rexptr(
    n = 0.0016 * 1.0016 * nrow(contact_data %>%
      filter(period == "Pre-pandemic")),
    lambda = dist_over_250$estimate[1],
    range = c(250, Inf)
  )),
  e_home = sample(
    size = 0.0016 * 1.0016 * nrow(contact_data %>%
      filter(period == "Pre-pandemic")),
    x = contact_data %>%
      filter(period == "Pre-pandemic") %>%
      pull(e_home)
  )
) %>%
  mutate(e_all = e_home + e_other, period = "Pre-pandemic", 
         idx = 0,
         id = NA,
         idx_id = paste(idx, id, sep = "_")
  )

# append to data
contact_data_adjusted <- contact_data %>% bind_rows(dat_append)

#### Estimate between-person NHH contact variance from CoMix panel ----
# Uses the panel structure of CoMix (repeated observations per person across fortnightly waves)
# to decompose variance in log(e_other + 0.5) into between- and within-person components.
# nhh_re_sigma: between-person log-SD, used as the SD of the individual log-normal multiplier
#               in run_model() (mean-corrected so E[multiplier] = 1, preserving period means).
# nhh_re_theta: NB dispersion for day-to-day within-person variation, estimated from CoMix.

nhh_panel <- contact_data %>%
  filter(!period %in% c("Pre-pandemic", "POLYMOD")) %>%
  mutate(person_id = str_extract(part_id, "\\d+$")) %>%
  filter(!is.na(person_id))

nhh_panel_log <- nhh_panel %>%
  group_by(period) %>%
  mutate(
    log_nh    = log(e_other + 0.5),
    log_nh_dm = log_nh - mean(log_nh)
  ) %>%
  ungroup()

# Between-person variance: variance of person-level means on the period-demeaned log scale.
# Restrict to persons with >=2 CoMix observations (required to separate between- from within-person).
person_means_nhh <- nhh_panel_log %>%
  group_by(person_id) %>%
  filter(n() >= 2) %>%
  summarise(person_mean_dm = mean(log_nh_dm), .groups = "drop")

nhh_re_sigma <- sqrt(var(person_means_nhh$person_mean_dm))

# NB dispersion for day-to-day within-person variation.
# Include person-level offset (log of person mean) to absorb between-person variance;
# theta from this model reflects residual within-person overdispersion only.
nhh_person_means_offset <- nhh_panel %>%
  group_by(person_id) %>%
  summarise(person_mean_contact = mean(e_other + 0.5), .groups = "drop")

nhh_nb_fit <- MASS::glm.nb(
  e_other ~ 0 + period + offset(log(person_mean_contact)),
  data = nhh_panel %>% left_join(nhh_person_means_offset, by = "person_id")
)
nhh_re_theta <- nhh_nb_fit$theta


##### KCL ANALYSIS ----
pickering <- readxl::read_xlsx(here::here("data", "pickering_dat.xlsx")) %>%
  select(-c(`Viral Growth`, ...7, ...8)) %>%
  rename("culture" = ...6) %>%
  mutate_at(
    .vars = vars(`SureScreen F`, Innova, Encode),
    .funs = function(x) ifelse(x == "ND", NA, x)
  ) %>%
  mutate_at(
    .vars = vars(`SureScreen F`, Innova, Encode),
    .funs = function(x) {
      case_when(
        x %in% c(0.5, 1, 2) ~ 1,
        is.na(x) ~ NA_real_,
        TRUE ~ 0
      )
    }
  ) %>%
  mutate(id = row_number()) %>%
  rename(ct = `Ct N1`) %>%
  mutate(vl = (-(ct - 44.34) / 3.134))

innova_mod <- glm(Innova ~ vl,
  data = pickering,
  family = "binomial"
)

innova_higher_mod <- glm(Innova ~ vl,
  data = pickering %>%
    mutate(vl = vl + 2.5), family = "binomial"
)

# sensitivity analysis on test probability
test_model_choice <- function(boolean) {
  if (boolean) {
    innova_higher_mod
  } else {
    innova_mod
  }
}

culture_mod <- glm(culture ~ vl, data = pickering, family = "binomial")

generate_params <- function(mod, n) {
  mu <- coef(mod)
  Sigma <- vcov(mod)
  res <- mvrnorm(n, mu, Sigma)
  return(res)
}

culture_prob <- function(vl, beta0, beta1) {
  1 / (1 + exp(-(beta0 + beta1 * vl)))
}

# sensitivity analysis for infectiousness as culture prob or lft prob
inf_model_choice <- function(boolean) {
  # browser()
  if (boolean) {
    innova_mod
  } else {
    culture_mod
  }
}

# Historical Rt estimates
rt <- readxl::read_excel("data/221123_R_and_growth_rate_time_series_for_publication_v1.0.xlsx",
  range = "Table1_-_R!B10:D52",
  col_names = c("date", "lower", "upper")
)

rt_by_time_period <- rt %>%
  filter(date < as.Date("2021-01-01")) %>%
  fuzzyjoin::fuzzy_inner_join(time_periods,
    by = c("date" = "date_start", "date" = "date_end"),
    match_fun = list(`>=`, `<=`)
  ) %>%
  group_by(across(-c(date, lower, upper))) %>%
  summarise(lo = mean(lower), hi = mean(upper)) %>%
  ungroup() %>%
  filter(period != "Lockdown 1") # very minimal overlap with lockdown 1 period (Rt starts 29/5/2020)

# Create viral load trajectories for a given number of sims
make_trajectories <- function(
    n_sims = 100,
    asymp_parms = asymp_fraction,
    variant_info,
    max_prolif = 14,
    max_clear = 30,
    max_peakvl = 40,
    browsing = FALSE) {
  if (browsing) browser()

  set.seed(seed)
  # simulate CT trajectories

  inf <- rbbinom(
    n = n_sims,
    size = 1,
    alpha = asymp_parms$shape1,
    beta = asymp_parms$shape2
  ) %>%
    as_tidytable() %>%
    rename.("asymptomatic" = x) %>%
    mutate.(
      sim = row_number.(),
      asymptomatic = as.logical(asymptomatic)
    )

  traj <- inf %>%
    crossing.(start = 0) %>%
    crossing(variant_info) %>%
    mutate.(
      prolif = case_when.(
        heterogen_vl ~ rnormTrunc(
          n = n(), mean = mean_prolif,
          sd = sd_prolif, min = 1, max = max_prolif
        ),
        TRUE ~ median(rnormTrunc(
          n = n(), mean = mean_prolif,
          sd = sd_prolif, min = 1, max = max_prolif
        ))
      ),
      clear = case_when.(
        heterogen_vl ~ rnormTrunc(
          n = n(), mean = mean_clear,
          sd = sd_clear, min = 1, max = max_clear
        ),
        TRUE ~ median(rnormTrunc(
          n = n(), mean = mean_clear,
          sd = sd_clear, min = 1, max = max_clear
        ))
      ),
      end = prolif + clear,
      onset_t = prolif + rnorm(n = n(), mean = 2, sd = 1.5)
    ) %>%
    select.(-c(mean_prolif, sd_prolif, mean_clear, sd_clear, clear)) %>%
    pivot_longer.(
      cols = -c(
        sim, variant, onset_t, asymptomatic, heterogen_vl,
        mean_peakvl, sd_peakvl
      ),
      values_to = "x"
    ) %>%
    mutate.(y = case_when.(
      name == "start" ~ 40,
      name == "end" ~ 40,
      name == "prolif" ~ case_when.(
        heterogen_vl ~ rnormTrunc(
          n = n(),
          mean = mean_peakvl,
          sd = sd_peakvl, min = 0, max = max_peakvl
        ),
        TRUE ~ median(rnormTrunc(
          n = n(),
          mean = mean_peakvl,
          sd = sd_peakvl, min = 0, max = max_peakvl
        ))
      )
    )) %>%
    select.(-c(mean_peakvl, sd_peakvl))


  models <- traj %>%
    nest.(data = -c(sim, variant, onset_t, asymptomatic, heterogen_vl)) %>%
    mutate.(
      # Perform approxfun on each set of points
      m = map.(data, ~ approxfun(x = .x$x, y = .x$y))
    )

  # cannot pivot wider with "m" column - extract and rejoin
  x_model <- models %>%
    select.(-data)

  models <- models %>%
    select.(-m) %>%
    unnest.(data, .drop = F) %>%
    select.(-c(y)) %>%
    pivot_wider.(names_from = name, values_from = x) %>%
    left_join.(x_model) %>%
    select.(c(sim, variant, heterogen_vl, onset_t, prolif, start, end, m)) %>%
    arrange.(sim)
}

inf_curve_func <- function(m, start = 0, end = 30, interval = 1) {
  # browser()
  x <- tidytable(t = seq(start, end, by = interval)) %>%
    mutate.(
      ct = m(t),
      vl = convert_Ct_logGEML(ct)
    )

  return(infectiousness = x)
}

calc_sensitivity <- function(model, x) {
  # browser()
  if (!is.na(x)) {
    s <- model(x)
  } else {
    s <- NA_real_
  }

  return(s)
}

propresponsible <- function(R0, k, prop) {
  qm1 <- qnbinom(1 - prop, k + 1, mu = R0 * (k + 1) / k)
  remq <- 1 - prop - pnbinom(qm1 - 1, k + 1, mu = R0 * (k + 1) / k)
  remx <- remq / dnbinom(qm1, k + 1, mu = R0 * (k + 1) / k)
  q <- qm1 + 1
  1 - pnbinom(q - 1, k, mu = R0) - dnbinom(q, k, mu = R0) * remx
}

sample_filter <- function(condition, df, col, n) {
  sample(df %>% filter.(period == condition) %>% pull.(col), size = n, replace = T)
}

mean_filter <- function(condition, df, col) {
  mean(df %>% filter.(period == condition) %>% pull.(col))
}

#### Main Model ----
run_model <- function(testing_scenarios, scenarios, contact_dat = contact_data,
                      traj_full = traj, traj_processed = traj_,
                      within_person_re = TRUE, browsing = F) {
  if (browsing) {
    browser()
  }

  #### Generate infections of hh (household) contacts ####
  indiv_params <- traj_full %>%
    select.(-m) %>%
    crossing.(
      heterogen_contacts = unique(scenarios$heterogen_contacts),
      period = unique(scenarios %>%
        mutate.(period = fct_drop(period)) %>%
        pull.(period))
    ) %>%
    mutate.(
      # Sample row indices jointly so each individual gets HH and NHH contacts
      # from the same survey respondent, preserving the empirical joint distribution.
      .row_idx = if (heterogen_contacts[1]) {
        rows <- which(as.character(contact_dat$period) == as.character(period[1]))
        sample(rows, size = n(), replace = TRUE)
      } else {
        rep(NA_integer_, n())
      },
      part_id = contact_dat$part_id[.row_idx],
      idx_id = contact_dat$idx_id[.row_idx],
      hh_contacts = if (heterogen_contacts[1]) {
        contact_dat$e_home[.row_idx]
      } else {
        rpois(n(), mean_filter(period[1], contact_dat, "e_home"))
      },
      # nhh_contacts = if (heterogen_contacts[1]) {
      #   contact_dat$e_other[.row_idx]
      # } else {
      #   rpois(n(), mean_filter(period[1], contact_dat, "e_other"))
      # },
      .by = c(period, heterogen_contacts)
    ) #%>%
    # select.(-.row_idx)

  indiv_params_long <- indiv_params %>%
    left_join.(traj_processed)

  # simulate infections (and keep first instance)
  lookup <- split(
    contacts_hh_duration$cnt_duration,
    contacts_hh_duration$part_id
  )

  lookup_period <- split(
    contacts_hh_duration$cnt_duration,
    contacts_hh_duration$period
  )

  indiv_expanded <- indiv_params_long %>%
    uncount(hh_contacts, .id = "id", .remove = FALSE)

  idx_part <- match(indiv_expanded$part_id, names(lookup))

  period_vec <- as.character(indiv_expanded$period)
  period_vec[period_vec != "Pre-pandemic"] <- "Pandemic"

  idx_period <- match(period_vec, names(lookup_period))

  hh_duration_vec <- numeric(length(idx_part))

  for (i in seq_along(idx_part)) {
    vals <- lookup[[idx_part[i]]]
    if (is.null(vals)) {
      vals <- lookup_period[[idx_period[i]]]
    }
    hh_duration_vec[i] <- vals[sample.int(length(vals), 1)]
  }

  hh_infections <- indiv_expanded %>%
    mutate(
      hh_duration = ifelse(
        heterogen_contacts,
        hh_duration_vec,
        median(contacts_hh_duration$cnt_duration, na.rm = TRUE)
      ),
      infected = rbernoulli(n(), p = 1 - exp(-beta_inf * culture_p * hh_duration))
    ) %>%
    filter.(infected == T) %>%
    slice.(min(t), .by = c(all_of(key_grouping_var), hh_contacts, id)) %>%
    count.(t, all_of(key_grouping_var), hh_contacts, name = "hh_infected") %>%
    arrange.(sim)

  rm(indiv_expanded)
  gc()

  #### Calculate nhh infections ####
  lookup <- split(
    contacts_nhh_duration$cnt_duration,
    contacts_nhh_duration$part_id
  )

  lookup_period <- split(
    contacts_nhh_duration$cnt_duration,
    contacts_nhh_duration$period
  )

  lookup_nhh <- split(contact_dat$e_other, contact_dat$idx_id)

  mean_nhh_contacts <- contact_dat %>%
    summarise(mean_e_other = mean(e_other, na.rm = TRUE), .by = period)
    
  indiv_expanded <- indiv_params_long %>%
    left_join.(mean_nhh_contacts, by = "period") %>%
    mutate.(idx = match(idx_id, names(lookup_nhh))) %>%
    mutate.(
      nhh_contacts = if (!within_person_re && heterogen_contacts[1]) {
        # within_person_re=FALSE: draw independently each day (old behaviour, no correlation)
        sample_filter(condition = period[1], df = contact_dat, col = "e_other", n = n())
      } else {
        # within_person_re=TRUE: 
        # if heterogeneous contacts, sample NHH contacts from individual who HH 
        # contacts were sampled from, to preserve the empirical HH/NHH joint
        # distribution
        # otherwise sample from Poisson distribution with mean number of NHH 
        # contacts for given time period
        if (heterogen_contacts[1]) {
          if (period[1] == "Pre-pandemic"){
            contact_dat$e_other[.row_idx]
          } else {
            vapply(idx, function(iid) {
              if (is.na(iid) || length(lookup_nhh[[iid]]) == 0) return(NA_real_)
              vals <- lookup_nhh[[iid]]
              vals[sample.int(length(vals), 1)]
            }, numeric(1))            
          }
        } else {
          rpois(n(), mean_e_other)
        }
      },
      .by = all_of(key_grouping_var)
    ) %>%
    uncount.(nhh_contacts, .remove = F)

  idx_part <- match(indiv_expanded$part_id, names(lookup))

  period_vec <- as.character(indiv_expanded$period)
  period_vec[period_vec != "Pre-pandemic"] <- "Pandemic"

  idx_period <- match(period_vec, names(lookup_period))

  nhh_duration_vec <- numeric(length(idx_part))

  for (i in seq_along(idx_part)) {
    vals <- lookup[[idx_part[i]]]
    if (is.null(vals)) {
      vals <- lookup_period[[idx_period[i]]]
    }
    nhh_duration_vec[i] <- vals[sample.int(length(vals), 1)]
  }

  nhh_infections <- indiv_expanded %>%
    # Simulate infections
    mutate.(
      nhh_duration = ifelse(
        heterogen_contacts,
        nhh_duration_vec,
        median(contacts_nhh_duration$cnt_duration)
      ),
      nhh_infected = rbernoulli(n = n(), p = 1 - exp(-beta_inf * culture_p * nhh_duration))
    ) %>%
    summarise.(nhh_infected = sum(nhh_infected), .by = c(t, all_of(key_grouping_var), nhh_contacts, test)) %>%
    # Testing: determine if and when testing + isolating by specified sampling frequency, adherence
    right_join.(testing_scenarios) %>%
    mutate.(
      test_day = case_when.(
        (t - begin_testing) %% sampling_freq == 0 ~ TRUE,
        nhh_contacts > event_size ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    mutate.(
      earliest_pos = min(t[test & test_day]),
      test_iso = t >= earliest_pos & self_iso_test,
      .by = c(all_of(key_grouping_var), prop_self_iso_test, self_iso_test, begin_testing, sampling_freq, event_size)
    ) %>%
    filter.(test_iso == F) %>%
    select.(everything(), -test_iso, -test, -earliest_pos, -test_day)

  rm(indiv_expanded)
  gc()

  # Join nhh and hh contacts and summarise
  processed_infections <- indiv_params_long %>%
    right_join.(testing_scenarios) %>%
    left_join.(hh_infections) %>%
    left_join.(nhh_infections) %>%
    replace_na.(list(hh_infected = 0, nhh_infected = 0, nhh_contacts = 0)) %>%
    arrange.(period, lower_inf_thresh) %>%
    mutate.(
      total_contacts = nhh_contacts + hh_contacts,
      total_infections = nhh_infected + hh_infected
    )
}

# Calibrate beta_inf so that mean pre-pandemic secondary cases equals target_R0.
# Uses a reduced number of simulations (n_calib) for speed during optimisation;
# the calibrated value is then saved and re-used by all downstream scripts.
calibrate_beta <- function(target_R0 = 2.5, n_calib = 2000, lower = 0.01, upper = 100) {
  calib_scenarios <- crossing(time_periods) %>%
    filter(period == "Pre-pandemic") %>%
    mutate(scenario_id = row_number()) %>%
    select(-c(date_start, date_end)) %>%
    crossing(heterogen_contacts = TRUE)

  traj_full_calib <- traj %>%
    filter.(heterogen_vl == TRUE, sim <= n_calib)

  # Build traj_processed internally — culture_p doesn't depend on beta_inf
  traj_calib_processed <- traj_full_calib %>%
    mutate.(infectiousness = pmap(inf_curve_func, .l = list(
      m = m, start = start, end = end, interval = 1
    ))) %>%
    unnest.(infectiousness) %>%
    crossing.(lower_inf_thresh = c(FALSE)) %>%
    mutate.(
      culture_p = culture_prob(vl, beta0, beta1),
      infectious = rbernoulli(n = n(), p = 1 - exp(-beta_inf * culture_p * median_contact_duration)),
      test_p = stats::predict(
        object = innova_mod, type = "response",
        newdata = tidytable(vl = vl)
      ),
      test = rbernoulli(n = n(), p = test_p),
      .by = c(lower_inf_thresh)
    ) %>%
    replace_na.(list(test = FALSE, infectious = FALSE)) %>%
    select.(-c(prolif, start, end))

  calib_testing <- traj_full_calib %>%
    select.(-m) %>%
    crossing.(
      prop_self_iso_test = 0,
      sampling_freq      = 7L,
      event_size         = NA_real_
    ) %>%
    mutate.(
      self_iso_test = FALSE,
      begin_testing = 0L
    )

  obj_fn <- function(beta) {
    beta_inf <<- beta
    set.seed(12345) # fixed seed so contact draws are identical every evaluation → smooth monotone objective
    res <- run_model(
      testing_scenarios = calib_testing,
      scenarios         = calib_scenarios,
      contact_dat       = contact_data,
      traj_full         = traj_full_calib,
      traj_processed    = traj_calib_processed,
      within_person_re  = TRUE,
      browsing          = FALSE
    )
    mean_R <- res %>%
      group_by(sim) %>%
      summarise(tot = sum(total_infections), .groups = "drop") %>%
      pull(tot) %>%
      mean()
    message(sprintf("  beta_inf = %.4f  =>  mean R = %.4f  (target %.2f)", beta, mean_R, target_R0))
    mean_R - target_R0
  }

  uniroot(obj_fn, interval = c(lower, upper), tol = 0.001)$root
}

# function to calculate the proportion above or below a defined threshold
prop_n <- function(df, threshold = 10, col = e_all, op = ">=") {
  browser()
  df %>%
    summarise.(
      n = n(),
      s = sum(match.fun(op)({{ col }}, {{ threshold }})),
      "prop_{{threshold}}" := s / n
    ) %>%
    rename.("n_{{threshold}}" := s)
}

# https://www.medrxiv.org/content/10.1101/2020.04.25.20079103v3
asymp_fraction <- rriskDistributions::get.beta.par(
  q = c(0.24, 0.38),
  p = c(0.025, 0.975),
  show.output = F, plot = F
) %>%
  as.list()

approx_sd <- function(x1, x2) {
  (x2 - x1) / (qnorm(0.95) - qnorm(0.05))
}

# bootstrap confidence interval function
boot_ci <- function(x, nrep = 100) {
  trueval <- tibble(
    param = c("mu", "size"),
    mean = c(
      x$estimate[[2]],
      x$estimate[[1]]
    )
  )

  ci <- bootdist(f = x, niter = nrep)$CI %>%
    as.data.frame() %>%
    select(-Median) %>%
    rownames_to_column("param")

  left_join(trueval, ci)
}

# run code quietly
hush <- function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp <- code
  sink()
  return(tmp)
}

# logarithmic spaced sequence
# from library("emdbook"), because need only this
lseq <- function(from = 1, to = 100000, length.out = 6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

# quantile function
quibble2 <- function(x, q = c(0.25, 0.5, 0.75)) {
  tibble("{{ x }}" := quantile(x, q), "{{ x }}_q" := q)
}
