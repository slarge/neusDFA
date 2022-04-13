# VAST attempt 1 univariate model as a script
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")

library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)

## For some reason I need to make sure Rtools has path properly set, else TMB won't compile
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

# Load Data -----

# spp_list <- c("calfin", "cham", "ctyp", "tlong")
season_list <- c("spring", "fall")

data("ecomon_epu")


## Number of stations per year
total_stations <- ecomon_epu %>%
  group_by(year) %>%
  summarize(total = n_distinct(id))

## Proportion of stations with positive tows per year. Cutoff used to "keep"
spps <- ecomon_epu %>%
  mutate(present = ifelse(abundance > 0,
                          1, 0)) %>%
  filter(year >= 1994) %>%
  left_join(total_stations) %>%
  group_by(year, season, spp) %>%
  summarize(positive_stations = sum(present, na.rm = TRUE)/total,
            keep = ifelse(positive_stations > .10,
                          1, 0)) %>%
  distinct()

## List of spp where the total number of years by season is greater than 5
## and each spp group has more than cutoff of positive tows
spp_list <- spps %>%
  filter(season %in% season_list) %>%
  group_by(spp, season) %>%
  summarize(count = sum(keep)) %>%
  filter(all(count > 5)) %>%
  select(spp) %>%
  distinct() %>%
  pull(spp)

#
# ggplot(ecomon_epu %>% filter(spp %in% spp_list,
#                              season %in% c("spring", "fall"),
#                              year > 1994,
#                              abundance > 0), aes(x = lon, y = lat, color = spp)) +
#   geom_point() +
#   facet_grid(season ~ year)


i <- "para"
j <- "spring"

for(i in spp_list){
  for(j in season_list){


    zoop_dat <- ecomon_epu %>%
      dplyr::filter(grepl(paste0("^",i), spp),
                    EPU %in% c("GB", "GOM", "MAB"),
                    season == j,
                    as.numeric(year) >= 1994) %>%
      dplyr::mutate(areaswept_km2 = 1,
                    catch_ab = abundance) %>%
      droplevels()


    # ggplot(zoop_dat, aes(x = lon, y = lat)) +
    #   geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
    #   geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = EPU, size = abundance), alpha = 0.5) +
    #   facet_wrap(~year) +
    #   NULL

    covariate_data <- with(zoop_dat, data.frame(Year = NA,
                                                Lat = lat,
                                                Lon = lon,
                                                bottom_depth = bottom_depth/100))

    # Vast Models -----

    ## Model settings in FishStatsUtils::make_settings()
    ## The following settings define whether to include spatial and
    ## spatio-temporal variation, whether its autocorrelated, and
    ## whether there is overdispersion. find more at ?VAST::make_data()

    # Model structure: presence/absence delta model (binomial and poisson)


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

    FieldConfig <- c("Omega1"   = "IID",
                     "Epsilon1" = "IID",
                     "Omega2"   = "IID",
                     "Epsilon2" = "IID")

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

    RhoConfig <- c("Beta1"    = c(0, 1, 2, 3, 4)[3],
                   "Beta2"    = c(0, 1, 2, 3, 4)[3],
                   "Epsilon1" = c(0, 1, 2, 3, 4)[3],
                   "Epsilon2" = c(0, 1, 2, 3, 4)[3])


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


    ObsModel <- c("PosDist" = 7, # Zero-inflated Poisson (1st linear predictor for logit-linked zero-inflation; 2nd linear predictor for log-linked conditional mean of Poisson)
                  "Link"    = 1) ## Might want this to be 0


    X1_formula <- ~ splines::bs(bottom_depth, knots = 3, intercept = FALSE)


    ## Make settings -----

    settings <- make_settings( n_x = 100,
                               Region = "northwest_atlantic",
                               strata.limits = "EPU",
                               # strata.limits = list('All_areas' = 1:1e5),# full area
                               purpose = "index2",
                               bias.correct = FALSE,
                               use_anisotropy = TRUE,
                               fine_scale = TRUE,
                               FieldConfig = FieldConfig,
                               RhoConfig = RhoConfig,
                               OverdispersionConfig = OverdispersionConfig
    )

    settings$epu_to_use <- c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")


    # Run model -----

    ## Spring model -----
    working_dir <- here::here(sprintf("analysis/vast_index/%s_%s/", i, j))

    if(!dir.exists(working_dir)) {
      dir.create(working_dir, recursive  = TRUE)
    }

    #run
    fit <- fit_model(settings = settings,
                     epu_to_use = settings$epu_to_use,
                     Lat_i = zoop_dat$lat,
                     Lon_i = zoop_dat$lon,
                     t_i = zoop_dat$year,
                     c_i = rep(0, nrow(zoop_dat)),
                     b_i= as_units(zoop_dat$abundance, "count"),
                     a_i = as_units(zoop_dat$areaswept_km2, "km^2"),
                     working_dir = working_dir,
                     covariate_data = covariate_data,
                     X1_formula = X1_formula,
                     X2_formula = ~ 0,
                     anisotropy = FALSE, # corresponds to ln_H_input params
                     test_fit = FALSE,
                     newtonsteps = 1,
                     getsd = TRUE,
                     run_model = TRUE,
                     Use_REML = TRUE,
                     knot_method = "grid",
                     getsd = TRUE,
                     # Method = Method,
                     optimize_args = list("lower" = -Inf,
                                          "upper" = Inf),
                     bias.correct.control = list("nsplit" = 100))



    # Plot results
    plot( fit,
          working_dir = paste0(working_dir, "/"))

  }
}

# yeah, two options: (1) rerun with newtonsteps=0 and getsd=FALSE and inspect the "parameter_estimates.txt" to learn more about what variance parameter is hitting its lower bound; (2) just ignore the "statistical pecadillo" of having a parameter at a bound, and try fit_model( ..., test_fit=FALSE) and see if the Hessian is sufficiently close to positive definite to proceed
# basically, VAST runs check_fit after running nlminb but before calculating SEs, and is conservative (risk-averse) in flagging potential statistical issues, including in this case that a parameter is approaching a bound.
# I could also walk through how to add if-then statements to automate model-changes in these cases ... its been on a back-burner for a while to expand TMBhelper::fit_tmb to automatically reduce models when parameters hit bounds


