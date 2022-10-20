# seasonal_dfa
library(MARSS)
library(tidyr)
# library(bayesdfa)
library(ggplot2)
library(dplyr)
library(future)
library(furrr)
# library(rstan)
source("R/helper_functions.R")

possibly_dfa_mod <- purrr::possibly(dfa_mod, otherwise = NA_character_)


index_ctl_seasonal <- data.frame(path = list.files(path = here::here("analysis/vast_seasonal_index"),
                                                      pattern = ".*ln_index_ctl.rds$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(spp = gsub(pattern = ".*vast_seasonal_index/(.*)_seasonal(.*)",  "\\1", path),
         dat = purrr::map(path, ~ readRDS(file = .x))) %>%
  tidyr::unnest(cols = c(dat)) %>%
  rename(epu = EPU,
         est = Estimate)

dats <- index_ctl_seasonal %>%
  group_by(epu, spp) %>%
  # filter(epu %in% c("Georges Bank", "Gulf of Maine")[1]) %>%
  # mutate(epu_spp = paste0(spp, "_", epu)) %>%
  # group_by(epu_spp) %>%
  mutate(est = scale(est, center = TRUE, scale = TRUE)[,1]) %>%
  select(#epu_spp,
         year,
         season,
         year_season,
         epu,
         spp,
         est)

t1 <- ggplot(dats, aes(x = year_season, y = est, group = spp, color = spp)) +
  geom_line(show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Estimated seasonal density",
       subtitle = "n = 14 zooplankton taxa",
       x = "Year + season",
       y = expression(z-scored~ln(abundance%.%km^{-2}))) +
  facet_wrap(~epu, ncol = 1)
#
#
# temp_wide <- dat %>%
#   ungroup() %>%
#   select(-epu,
#          -spp,
#          -year,
#          -season) %>%
#   tidyr::pivot_wider(names_from = year_season, values_from = est) %>%
#   select(-epu_spp) %>%
#   ungroup() %>%
#   as.matrix()
#
# dat_time <- as.double(unique(dat$year_season))
# dat_name <- as.character(unique(dat$epu_spp))
#
# row.names(temp_wide) <- dat_name
#
# temp_wide <- read.csv(here::here("analysis/data/derived_data/zoop_temp_wide.csv"))
# dat_time <- as.numeric(gsub("X", "", colnames(temp_wide[-1])))
# dat_name <- as.character(rownames(temp_wide))
#
#
# cos_t <- cos(2 * pi * seq(dat_time)/4)
# sin_t <- sin(2 * pi * seq(dat_time)/4)
# season_f <- rbind(cos_t, sin_t)
#
#
# row.names(temp_wide) <- temp_wide$X
# temp_wide <- as.matrix(temp_wide[,-1])
#
# dfa_mat <- expand_grid(m = 1:4,#1:nrow(dat),
#                        R = c("diagonal and unequal", "diagonal and equal", "unconstrained" )[1:2],
#                        covariate = c("none", "season_f")[2])
#
#
# possibly_dfa_mod <- purrr::possibly(dfa_mod, otherwise = NA_character_)
#
# # plan(multisession, workers = parallelly::availableCores() -8)
# dfa_out <- dfa_mat %>%
#   # head(1) %>%
#   mutate(
#     # tt = purrr::pmap(.l = list(m, R, covariate), .f = function(m, R, covariate) sprintf("%s and %s and %s", m, R, covariate))#,
#     # mod = furrr::future_pmap(.l = list(m, R, covariate), .f = function(m, R, covariate) possibly_dfa_mod(m = m, R = R, covariate = covariate, just_testing = TRUE)),
#     mod = purrr::pmap(.l = list(m, R, covariate), .f = function(m, R, covariate) possibly_dfa_mod(m = m, R = R, covariate = covariate, just_testing = TRUE)),
#     AICc = purrr::map(mod, "AICc")
#   ) %>%
#   arrange(AICc)

# plan(sequential)


dat_time <- as.numeric(unique(dats$year_season))
# dat_name <- as.character(rownames(temp_wide))


cos_t <- cos(2 * pi * seq(dat_time)/4)
sin_t <- sin(2 * pi * seq(dat_time)/4)
season_f <- rbind(cos_t, sin_t)
cov_v = season_f
## First check out the base models without covariates
epu_long <- dats %>%
  select(-year,
         -season) %>%
  # group_by(EPU, common_name) %>%
  # mutate(count = sum(!is.na(cond))) %>%
  # filter(!all(is.na(cond)),
  #        count >= 20) %>%
  # ungroup() %>%
  nest(data = -epu) %>%
  left_join(expand_grid(epu = c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"),
                        m = c(1:4, 8, 10)[1:2],
                        R = c("diagonal and unequal", "diagonal and equal", "unconstrained")[2],
                        covariate = c("none", "season_f")[2]),
            by = "epu") %>%
  ungroup() %>%
  mutate(cov_v = ifelse(covariate == "none", NA, list(season_f)))

plan(multisession, workers = parallelly::availableCores() - 6)
dfa_out <- epu_long %>%
  # head(6) %>%
  mutate(
    mod = furrr::future_pmap(.l = list(data, m, R, cov_v),
                      .f = function(data, m, R, cov_v) dfa_mod(dat = data, m = m, R = R, cov_v = cov_v,
                                                               just_testing = TRUE,
                                                               data_wide = FALSE)),
    AICc = purrr::map(mod, "AICc")) %>%
  arrange(AICc)

saveRDS(dfa_out, file = "analysis/zoop_seasonal_dfa_epu.rds")
plan(sequential)


table_mod <- dfa_out %>%
  select(-mod,
         -cov_v,
         -data) %>%
  tidyr::unnest(AICc) %>%
  group_by(epu) %>%
  arrange(AICc, .by_group = TRUE) %>%
  slice_min(order_by = AICc)


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
         rotate = TRUE) +
  theme_minimal() -> s_plot

ggsave(filename = here::here("dfa_zoop.png"), s_plot)

### Plot of observations
autoplot(best_mod,
         plot.type = "fitted.ytT", #"ytT",
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
  geom_pointrange(show.legend = TRUE, position = position_dodge(width = 0.5)) +
  facet_wrap(~epu) +
  labs(color = "State",
       x = "",
       y = "factor loadings") +
  coord_flip() +
  theme_bw() -> l_plot

ggsave(filename = here::here("dfa_zoop_loadings.png"), l_plot)


