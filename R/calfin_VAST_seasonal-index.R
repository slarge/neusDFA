# VAST attempt 1 univariate model with seasonal effects
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")
library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)

## For some reason I need to make sure Rtools has path properly set, else TMB won't compile
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Load Data -----

data("ecomon_epu")

zoop_dat <- ecomon_epu %>%
  dplyr::filter(grepl("^calfin", spp),
                # EPU %in% c("GB", "GOM", "MAB"),
                as.numeric(year) >= 1994) %>%
  dplyr::mutate(areaswept_km2 = 1,
                catch_ab = abundance) %>%
  group_by(year, season) %>%
  slice_sample(prop = .5) %>%
  droplevels() %>%
  data.frame()

ecomon_epu <- NULL

# ggplot(zoop_dat, aes(x = lon, y = lat)) +
#   geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
#   geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = EPU, size = abundance), alpha = 0.5) +
#   facet_wrap(~year) +
#   NULL

# covariate_data <- with(zoop_dat, data.frame(Year = NA,
#                                             Lat = lat,
#                                             Lon = lon,
#                                             bottom_depth = bottom_depth/100))


# Load data and quick exploration of structure
# Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
year_set = sort(unique(zoop_dat[,'year']))
season_set = c("winter", "spring", "summer", "fall")

# Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
yearseason_grid = expand.grid("season" = season_set, "year" = year_set)
yearseason_levels = apply(yearseason_grid[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_labels = round(yearseason_grid[,'year'] + (as.numeric(factor(yearseason_grid[,'season'], levels = season_set))-1)/length(season_set), digits=1)

# Similar process, but for the observations
yearseason_i = apply(zoop_dat[,c("year","season")], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_i = factor(yearseason_i, levels = yearseason_levels)

# Add the year_season factor column to our sampling_data data set
zoop_dat$year_season = yearseason_i
zoop_dat$season = factor(zoop_dat$season, levels = season_set)

# Some last processing steps
zoop_dat = zoop_dat[, c("year", "season", "year_season", "lat", "lon", "areaswept_km2", "abundance")]

# Make dummy observation for each season-year combination
dummy_data = data.frame(
  year = yearseason_grid[,'year'],
  season = yearseason_grid[,'season'],
  year_season = yearseason_levels,
  lat = mean(zoop_dat[,'lat']),
  lon = mean(zoop_dat[,'lon']),
  areaswept_km2 = mean(zoop_dat[,'areaswept_km2']),
  abundance = 0,
  dummy = TRUE)

# Combine with sampling data
full_data = rbind(cbind(zoop_dat, dummy = FALSE), dummy_data)

# Create sample data
samp_dat = data.frame(
  "year_season" = as.numeric(full_data$year_season)-1,
  "Lat" = full_data$lat,
  "Lon" = full_data$lon,
  "abundance" = full_data$abundance,
  "areaswept_km2" = full_data$areaswept_km2,
  "Dummy" = full_data$dummy )

# Covariate data. Note here, case sensitive!
cov_dat = data.frame(
  "Year" = as.numeric(full_data$year_season)-1,
  "Year_Cov" = factor(full_data$year, levels = year_set),
  "Season" = full_data$season,
  "Lat" = full_data$lat,
  "Lon" = full_data$lon)

# Inspect
table("year_season"=cov_dat$Year, "Actual_year"=cov_dat$Year_Cov)
table("year_season"=cov_dat$Year, "Actual_season"=cov_dat$Season)

#####
## Model settings
#####


##  Random Fields -----

## Control the random fields part of the model.
## Omega = X is the number of random spatial fields to apply
## and Epsilon = X is the number of random spatio-temporal
## fields to apply. Omega1 is for the probability of
## occurrence, and Omega2 is for the density given occurrence,
## similarly for Epsilon.

## 0 = off
## "AR1" = AR1 process
## >1 = number of elements in a factor-analysis covariance
## "IID" = random effect following an IID distribution

FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)


## Autoregressive structure -----

## Control autoregressive structure for parameters over time
## Changing the settings here creates different
## autoregressive models for the intercept (Beta) and
## spatio-temporal process (Epsilon).

## 0 = each year is a fixed effect
## 1 = random effect
## 2 = random walk
## 3 = fixed effect that is constant over time
## 4 = AR1 process

RhoConfig <- c("Beta1"    = c(0, 1, 2, 3, 4)[4],
               "Beta2"    = c(0, 1, 2, 3, 4)[4],
               "Epsilon1" = c(0, 1, 2, 3, 4)[5],
               "Epsilon2" = c(0, 1, 2, 3, 4)[5])

## Correlated overdispersion -----

## Control correlated overdispersion among categories
## for each level of v_i, where eta1 is for encounter
## probability, and eta2 is for positive catch rates
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

## 0 = off,
## "AR1" = AR1 process,
## >0 = number of elements in a factor-analysis covariance
OverdispersionConfig <- c("eta1" = 0,
                          "eta2" = 0)


## Observation model -----
# Control observation model structure. The first
# component sets the distribution of the positive
# distribution component. ?VAST::make_data()

ObsModel <- c("PosDist" = 2, # Delta-Gamma; Alternative "Poisson-link delta-model" using log-link for numbers-density and log-link for biomass per number
              "Link"    = 1)

# Make settings
settings = make_settings(n_x = 100,
                         Region = "northwest_atlantic",
                         strata.limits = "EPU",
                         purpose = "index2",
                         FieldConfig = FieldConfig,
                         RhoConfig = RhoConfig,
                         ObsModel = ObsModel,
                         bias.correct = FALSE,
                         Options = c('treat_nonencounter_as_zero' = TRUE) )

settings$epu_to_use <- c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")

# Creating model formula
formula_use = ~ Season + Year_Cov

# Implement corner constraint for linear effect but not spatially varying effect:
# * one level for each term is 2 (just spatially varying) spatially varying, zero-centered linear effect on 1st linear predictor for category c
# * all other levels for each term is 3 (spatialy varying plus linear effect), spatially varying linear effect on 1st linear predictor for category c
X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1),
                            2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )
X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1),
                            2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )

#####
## Model fit -- make sure to use new functions
#####

## Spring model -----
working_dir <- here::here("analysis/vast_index/calfin_seasonal")

if(!dir.exists(working_dir)) {
  dir.create(working_dir, recursive  = TRUE)
}

fit_orig = fit_model(settings = settings,
                     Lat_i = samp_dat$Lat,
                     Lon_i = samp_dat$Lon,
                     t_i = samp_dat$year_season,
                     b_i =  as_units(samp_dat$abundance, "count"),
                     a_i = as_units(samp_dat$areaswept_km2, "km^2"),
                     epu_to_use = settings$epu_to_use,
                     working_dir = working_dir,
                     X1config_cp = X1config_cp_use,
                     X2config_cp = X2config_cp_use,
                     covariate_data = cov_dat,
                     X1_formula = formula_use,
                     X2_formula = formula_use,
                     X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE),
                                        Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
                     run_model = FALSE,
                     PredTF_i = samp_dat$Dummy,
                     Use_REML = FALSE,
                     getsd = TRUE,
                     test_fit = FALSE,
                     Options = c('treat_nonencounter_as_zero' = TRUE),
                     # Method = Method,
                     optimize_args = list("lower" = -Inf,
                                          "upper" = Inf))

# Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
Map_adjust = fit_orig$tmb_list$Map

# Pool variances for each term to a single value
Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]),
                                          nlevels(cov_dat$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season)+1]),
                                          nlevels(cov_dat$Year_Cov))))
Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]),
                                          nlevels(cov_dat$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season)+1]),
                                          nlevels(cov_dat$Year_Cov))))

# Fit final model with new mapping
fit = fit_model(settings = settings,
                Lat_i = samp_dat$Lat,
                Lon_i = samp_dat$Lon,
                t_i = samp_dat$year_season,
                b_i =  as_units(samp_dat$abundance, "count"),
                a_i = as_units(samp_dat$areaswept_km2, "km^2"),
                epu_to_use = settings$epu_to_use,
                working_dir = working_dir,
                X1config_cp = X1config_cp_use,
                X2config_cp = X2config_cp_use,
                covariate_data = cov_dat,
                X1_formula = formula_use,
                X2_formula = formula_use,
                X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE),
                                   Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
                newtonsteps = 1,
                PredTF_i = samp_dat$Dummy,
                Map = Map_adjust,
                Use_REML = FALSE,
                getsd = TRUE,
                test_fit = FALSE,
                Options = c('treat_nonencounter_as_zero' = TRUE),
                # Method = Method,
                optimize_args = list("lower" = -Inf,
                                     "upper" = Inf))

saveRDS(fit, file = "analysis/vast_index/calfin_seasonal/fit.rds")
#
plot( fit,
      projargs='+proj=natearth +lon_0=-68 +units=km',
      country = "united states of america",
      year_labels = yearseason_labels )

names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]

