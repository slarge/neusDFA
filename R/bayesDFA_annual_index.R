# seasonal_dfa
library(MARSS)
library(tidyr)
# library(bayesdfa)
library(ggplot2)
library(dplyr)
library(future)
# library(rstan)
source("R/helper_functions.R")

index_ctl_annual_list <- data.frame(path = list.files(path = here::here("analysis/vast_seasonal_index/"),
                                                      pattern = ".*index_ctl_annual.rds$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(spp = gsub(pattern = ".*vast_seasonal_index/(.*)_seasonal(.*)",  "\\1", path),
         dat = purrr::map(path, ~readRDS(file = .x))) %>%
  tidyr::unnest(cols = c(dat))

dat <- index_ctl_annual_list %>%
  # filter(epu %in% c("Georges Bank", "Gulf of Maine")[1]) %>%
  mutate(epu_spp = paste0(spp, "_", epu)) %>%
  group_by(epu_spp) %>%
  mutate(est = scale(est, center = TRUE, scale = TRUE)) %>%
  select(epu_spp,
         # epu,
         # spp,
         year,
         est)

ggplot(dat, aes(x = year, y = est, group = epu_spp, color = epu_spp)) +
  geom_line(show.legend = FALSE) +
  theme_minimal()


temp_wide <- dat %>%
  tidyr::pivot_wider(names_from = year, values_from = est) %>%
  ungroup() %>%
  select(-epu_spp) %>%
  as.matrix()

dat_time <- as.integer(unique(dat$year))
dat_name <- as.character(unique(dat$epu_spp))

row.names(temp_wide) <- dat_name


dfa_mat <- expand_grid(m = 1:4,#1:nrow(dat),
                       R = c("diagonal and unequal", "diagonal and equal", "unconstrained" ),
                       covariate = c("none", "season_f")[1])


possibly_dfa_mod <- purrr::possibly(dfa_mod, otherwise = NA_character_)

plan(multisession, workers = parallelly::availableCores() -1)
dfa_out <- dfa_mat %>%
  # head(2) %>%
  mutate(mod = furrr::future_pmap(.l = list(m, R, covariate), .f = possibly_dfa_mod),
         AICc = purrr::map(mod, "AICc")) %>%
  arrange(AICc)

plan(sequential)


best_mod <- dfa_out %>%
  tidyr::unnest(AICc) %>%
  arrange(AICc) %>%
  head(1) %>%
  pull(mod)


table_mod <- dfa_out %>%
  select(-mod) %>%
  tidyr::unnest(AICc) %>%
  arrange(AICc)

best_mod <- best_mod[[1]]

###  make plot of states and CIs
# “model.ytT”, “xtT”, “model.resids”, “state.resids”, “qqplot.model.resids”, “qqplot.state.resids”, “ytT”, “acf.model.resids”
autoplot(best_mod,
         plot.type = "xtT",
         form = "dfa",
         rotate = TRUE)

### Plot of observations
autoplot(best_mod,
         plot.type = "model.ytT", #"ytT",
         form = "dfa",
         rotate = TRUE)

### Plot of residuals
autoplot(best_mod,
         plot.type = "acf.model.resids",
         form = "dfa",
         rotate = TRUE)

best_mod <- MARSSparamCIs(best_mod)
# the rotation matrix for the Z
z <- coef(best_mod, type = "Z")
H.inv <- varimax(z)$rotmat

# Get the Z, upZ, lowZ
# Z.up <- coef
z.low <- coef(best_mod, type = "Z", what="par.lowCI")
z.up <- coef(best_mod, type = "Z", what="par.upCI")

z.rot <- z %*% H.inv
z.rot.up <- z.up %*% H.inv
z.rot.low <- z.low %*% H.inv

df <- data.frame(name_code = rep(row.names(temp_wide), table_mod$m[1]),
                 trend = rep(1:table_mod$m[1], each = length(dat_name)),
                 Z = as.vector(z.rot),
                 Zup = as.vector(z.rot.up),
                 Zlow = as.vector(z.rot.low)) %>%
  mutate(epu = gsub("^.*_(.*)",  "\\1", name_code),
         spp = gsub("^(.*)_.*$",  "\\1", name_code))



ggplot(data = df,
       aes(x = spp,
           ymin = Zlow,
           ymax = Zup,
           y = Z,
           color = epu)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.5)) +
  facet_wrap(~trend) +
  coord_flip() +
  theme_bw()



ggplot(data = df,
       aes(x = spp,
           ymin = Zlow,
           ymax = Zup,
           y = Z,
           color = as.factor(trend))) +
  geom_hline(yintercept = 0, color = "black") +
  geom_pointrange(show.legend = FALSE, position = position_dodge(width = 0.5)) +
  facet_wrap(~epu) +
  coord_flip() +
  theme_bw()



