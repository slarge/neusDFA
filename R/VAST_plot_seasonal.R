library(VAST)
library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)
library(ggrepel)

source("R/helper_functions.R")

fit_list <- list.files(path = here::here("analysis/vast_seasonal_index/"), pattern = "fit.rds", recursive = TRUE, full.names = TRUE)
fit_names <- gsub(pattern = ".*vast_seasonal_index/(.*)_seasonal/fit.rds$",  "\\1", fit_list)

# file.rename(list.files(path = here::here("analysis/vast_seasonal_index/"),
#                        pattern = ".*Kmeans", recursive = TRUE, full.names = TRUE),
#             stringr::str_replace(list.files(path = here::here("analysis/vast_seasonal_index/"),
#                                             pattern = ".*Kmeans", recursive = TRUE, full.names = TRUE),
#                                  pattern = "seasonalKmeans", "seasonal/Kmeans"))

spp_list <- list(`Calanus finmarchicus` = "calfin",
                 `Chaetognatha` = "chaeto",
                 `Centropages hamatus` = "cham",
                 `Clausocalanus arcuicornis` = "clauso",
                 `Centropages typicus` = "ctyp",
                `Euphausiacea` = "euph",
                "euph1",
                `Gastropoda` = "gas",
                `Hyperiidea` = "hyper",
                "larvaceans",
              `Metridia lucens` = "mlucens",
              `Oithona spp.` = "oithspp",
              `Paracalanus parvus` = "para",
              `Pseudocalanus spp.` = "pseudo",
              `Temora longicornis` = "tlong")





plot_seasonal <- function(spp = c( "calfin", "chaeto", "cham", "clauso", "ctyp", "euph",
                                   "gas", "hyper", "larvaceans", "mlucens", "oithspp", "para",
                                   "pseudo", "tlong")) {
  #

  # fit <- readRDS(fit_list[2])

  fit <- readRDS(grep(pattern = spp, fit_list, value = TRUE))

  d_gt <- data.frame(log(fit$Report$D_gct[,1,]),
                     fit$spatial_list$latlon_g) %>%
    dplyr::rename_with(~ yearseason_levels, .cols = !c("Lat", "Lon")) %>%
    tidyr::pivot_longer(cols =  !c("Lat", "Lon"),
                        names_to = "year_season",
                        values_to = "D") %>%
    mutate(D = as.numeric(D),
           season = factor(sub(".*?_", "", year_season), levels = c("winter", "spring", "summer", "fall")),
           year = as.numeric(sub("?_.*", "", year_season))) %>%
    sf::st_as_sf(coords = c("Lon","Lat")) %>%
    sf::st_set_crs(crs_epu)

  # names(grep(spp, spp_list, value = TRUE))
  spp_label = names(grep(spp, spp_list, value = TRUE))
  ## Density plot (by year)
  d_plot <- ggplot() +
    geom_sf(data = d_gt %>% filter(year %in% c(1995, 2000, 2005, 2010, 2015)), aes(color = D), alpha = 0.8) +
    geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
    geom_sf(data = ne_states, color = "grey40", size = 0.05) +
    coord_sf(crs = crs_epu, xlim = xlims, ylim = ylims) +
    facet_grid(year ~ season) +
    theme_bw() +
    scale_color_viridis_c() +
    labs(title = substitute(italic(x), list(x = spp_label)),
         subtitle = "Estimated density",
         x = "longitude",
         y = "latitude",
         color = expression(ln(number%.%km^{-2})))

  # d_plot

  ggsave(filename = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-density_plot.png", spp, spp)), d_plot)
  saveRDS(object = d_gt, file = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-d_gt.rds", spp, spp)))
  ## Spatio-temporal variation


  ## Index with Standard Error
  index_list <- list.files(path = here::here("analysis/vast_seasonal_index/"), pattern = "index_ctl_sd_array.rds", recursive = TRUE, full.names = TRUE)

  index_ctl_sd_array <- readRDS(grep(pattern = spp, index_list, value = TRUE))

  ## c, t, l, n
  epu <-  c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")
  dimnames(index_ctl_sd_array) <- list("c" = "c", "year_season" = yearseason_levels, "epu" = epu, "n" = paste0("X", 1:100))

  index_ctl_sd <- as.data.frame.table(index_ctl_sd_array) %>%
    mutate(season = factor(sub(".*?_", "", year_season), levels = c("winter", "spring", "summer", "fall")),
           year = as.numeric(sub("?_.*", "", year_season))) %>%
    group_by(year, epu, n) %>%
    summarize(annual_sum = sum(Freq, na.rm = TRUE)) %>% ## Annual value for each sample
    ungroup() %>%
    group_by(year, epu) %>%
    summarize(#annual_mean = mean(annual_sum, na.rm = TRUE), #
      sd = sd(annual_sum, na.rm = TRUE)/sqrt(100))#,
  # upper = annual_mean + sd,
  # lower = annual_mean - sd)


  index_ctl_array <- fit$Report$Index_ctl
  dimnames(index_ctl_array) <- list("c" = "c", "year_season" = yearseason_levels, "epu" = epu)

  index_ctl_annual <- as.data.frame.table(index_ctl_array) %>%
    mutate(#season = factor(sub(".*?_", "", year_season), levels = c("winter", "spring", "summer", "fall")),
      year = as.numeric(sub("?_.*", "", year_season)),
      Freq = as.numeric(Freq)) %>%
    group_by(year, epu) %>%
    summarize(est = mean(Freq, na.rm = TRUE)) %>%
    left_join(index_ctl_sd) %>%
    mutate(lower = est - sd,
           upper = est + sd)


  levels(index_ctl_annual$epu) <-list(`Georges Bank` = "Georges_Bank",
                                      `Mid-Atlantic Bight` = "Mid_Atlantic_Bight",
                                      `Gulf of Maine` = "Gulf_of_Maine" )


  p_index_ctl_annual <- ggplot(index_ctl_annual, aes(x = year, y = log(est), ymin = log(lower), ymax = log(upper), group = epu)) +
    geom_ribbon(aes(fill = epu), alpha = 0.3) +
    geom_line(aes(color = epu)) +
    geom_point(aes(color = epu)) +
    theme_bw() +
    labs(title = substitute(italic(x), list(x = spp_label)),
         subtitle = "Estimated annual density",
         x = "Year",
         y = expression(ln(number%.%km^{-2}))) +
    # facet_wrap(~epu) +
    NULL

  ggsave(filename = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-index_annual.png", spp, spp)), p_index_ctl_annual)
  saveRDS(object = index_ctl_annual, file = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-index_ctl_annual.rds", spp, spp)))

  ## Seasonal
  index_ctl_season <- as.data.frame.table(index_ctl_array) %>%
    mutate(season = factor(sub(".*?_", "", year_season), levels = c("winter", "spring", "summer", "fall")),
           year = as.numeric(sub("?_.*", "", year_season)),
           Freq = as.numeric(Freq)) %>%
    group_by(year, epu) %>%
    summarize(est = mean(Freq, na.rm = TRUE)) %>%
    left_join(index_ctl_sd) %>%
    mutate(lower = est - sd,
           upper = est + sd)


  TmbData <- fit$tmb_list$Obj$env$data
  sdreport <- fit$parameter_estimates$SD

  SD <- TMB::summary.sdreport(sdreport)
  SD_estimate <- TMB:::as.list.sdreport(sdreport, what = "Estimate", report = TRUE)
  SD_stderr <- TMB:::as.list.sdreport(sdreport, what = "Std. Error", report = TRUE)

  names(TmbData)[grepl('_gc|_gct', x=names(fit$Report))]
  epu <-  c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")

  index_ctl <- array(NA, dim = c(unlist(TmbData[c("n_t", "n_l")]), 2), dimnames = list("year_season" = yearseason_levels,
                                                                                       "epu" = c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight"),
                                                                                       c("Estimate", "Std. Error")))

  index_ctl[] = SD[which(rownames(SD) == "ln_Index_ctl"),
                   c("Estimate", "Std. Error")]

  ln_index_ctl <- data.frame(yearseason_levels, index_ctl, row.names = NULL) %>%
    rename_all(~ stringr::str_replace_all(.,"\\.\\.","")) %>%
    tidyr::pivot_longer(cols = -yearseason_levels, names_sep = "\\.", names_to = c("EPU", "variable")) %>%
    tidyr::pivot_wider(names_from = variable, values_from = value) %>%
    mutate(EPU = factor(EPU, levels = c( "Gulf_of_Maine", "Georges_Bank", "Mid_Atlantic_Bight")),
           upper = Estimate + StdError,
           lower = Estimate - StdError,
           season = factor(sub(".*?_", "", yearseason_levels), levels = c("winter", "spring", "summer", "fall")),
           year = as.numeric(sub("?_.*", "", yearseason_levels)),
           year_season = as.numeric(case_when(grepl("_spring", yearseason_levels) ~ year + .2,
                                              grepl("_summer", yearseason_levels) ~ year + .5,
                                              grepl("_fall", yearseason_levels) ~ year + .8,
                                              TRUE ~ year)))

  levels(ln_index_ctl$EPU) <- list(`Gulf of Maine` = "Gulf_of_Maine",
                                   `Georges Bank` = "Georges_Bank",
                                   `Mid-Atlantic Bight` = "Mid_Atlantic_Bight")

  p_index_ctl_seasonal <- ggplot(data = ln_index_ctl, aes(group = EPU)) +
    geom_ribbon(aes(x = season, ymin = lower, ymax = upper, fill = EPU), alpha = 0.4) +
    geom_line(aes(x = season, y = Estimate, color = EPU)) +
    geom_point(aes(x = season, y = Estimate, color = EPU)) +
    labs(title = substitute(italic(x), list(x = spp_label)),
         subtitle = "Estimated seasonal density",
         x = "Year",
         y = expression(ln(numbers%.%km^{-2}))) +
    theme_bw() +
    facet_wrap(~as.factor(year))


  p_index_ctl_seasonal_yr <- ggplot(data = ln_index_ctl) +
    geom_ribbon(aes(x = year_season, ymin = lower, ymax = upper, fill = EPU), alpha = 0.4, show.legend = FALSE) +
    geom_line(aes(x = year_season, y = Estimate, color = EPU), show.legend = FALSE) +
    labs(title = substitute(italic(x), list(x = spp_label)),
         subtitle = "Estimated seasonal density",
         x = "Year",
         y = expression(ln(abundance%.%km^{-2}))) +
    theme_bw() +
    facet_wrap(~EPU, nrow = 3)


  ggsave(filename = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-index_seasonal.png", spp, spp)), p_index_ctl_seasonal)
  ggsave(filename = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-index_seasonal_yr.png", spp, spp)), p_index_ctl_seasonal_yr)

  saveRDS(object = ln_index_ctl, file = here::here(sprintf("analysis/vast_seasonal_index/%s_seasonal/%s-ln_index_ctl.rds", spp, spp)))


}


purrr::map(fit_names[11:14], plot_seasonal)



# # Similar process, but for the observations
# yearseason_i = apply(zoop_dat[,c("year","season")], MARGIN = 1, FUN = paste, collapse = "_")
# yearseason_i = factor(yearseason_i, levels = yearseason_levels)
#
# # Add the year_season factor column to our sampling_data data set
# zoop_dat$year_season = yearseason_i
# zoop_dat$season = factor(zoop_dat$season, levels = season_set)

# mdl <- make_map_info(Region = fit$settings$Region,
#                      spatial_list = fit$spatial_list,
#                      Extrapolation_List = fit$extrapolation_list)
#
#
# ak_map <- subset(map_data("world"), region=='USA' & subregion %in% c("Massachusettes", "Maine", "New York"))
#
# ## Have to duplicate it for each year so can facet below
# ak_map <- cbind(ak_map[rep(1:nrow(ak_map), times=nyrs),],
#                 Year=rep(yearseason_levels, each=nrow(ak_map)))
#
#
# gmap <- ggplot(ak_map, aes(x = long, y = lat, group = group)) +
#   geom_polygon(fill="black", colour = "white") +
#   scale_color_viridis_c(option = "magma") +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         panel.spacing.x=unit(0, "lines"),
#         panel.spacing.y=unit(0, "lines"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank() ) +
#   coord_cartesian(xlim=mdl$Xlim, ylim=mdl$Ylim)

## Below shows to you get the model estimate of density, D_gct,
## for each grid (g), category (c; not used here single
## univariate); and year (t); and link it spatially to a lat/lon
## extrapolation point.  You can do this for any _gct or _gc
## variable in the Report.
# names(fit$Report)[grepl('_gc|_gct', x=names(fit$Report))]
#
#
# D_gt <- fit$Report$D_gct[,1,] # drop the category
# dimnames(D_gt) <- list(cell=1:nrow(D_gt), year_season = yearseason_levels)
# ## tidy way of doing this, reshape2::melt() does
# ## it cleanly but is deprecated
# D_gt <- D_gt %>% as.data.frame() %>%
#   tibble::rownames_to_column(var = "cell") %>%
#   tidyr::pivot_longer(-cell, names_to = "year_season", values_to='D') %>%
#   mutate(D = as.numeric(D),
#          year_season = as.factor(year_season),
#          season = sub(".*?_", "", year_season),
#          year = as.numeric(sub("?_.*", "", year_season)),
#          index = as.numeric(year_season)-1)
#
# D <- merge(D_gt, mdl$PlotDF, by.x='cell', by.y='x2i')
#
# test_dat <-  D %>% filter(year == 2017) %>%
#   mutate(log_d = log(D)) %>%
#   select(Lat, Lon, log_d, season, year) %>%
#   distinct()
# head(test_dat)
# # g <-
# # gmap +
#
#   ggplot() +
#   geom_point(data = test_dat, aes(Lon, Lat, color = log_d, group=NULL),
#              size=.3, stroke=0,shape=16) +
#              ## These settings are necessary to avoid
#              ## overlplotting which is a problem here. May need
#              ## to be tweaked further.
#   facet_wrap(~season)
#
#



# results = plot_results(settings = settings,
#                        fit = fit,
#                        plot_set = 8,
#                        # plot_set = c(3, 11, 13, 14),
#                        working_dir = plot_dir,
#                        # category_names = "piscivores",
#                        year_labels = fit$year_labels,
#                        check_residuals = TRUE)


# ne_strata <- sf::read_sf("analysis/data/shapefiles/BTS_Strata.shp") %>%
#   sf::st_transform(crs = crs)

#
# td <- te %>%
#   sfc_as_cols(names = c("lon", "lat")) %>%
#   sf::st_drop_geometry()
#
#
# ggplot() +
#   geom_sf(data = te %>% filter(year %in% 2010:2015), aes(color = condition)) +
#   geom_sf(data = ne_countries, color = "grey40", size = 0.05) +
#   geom_sf(data = ne_states, color = "grey40", size = 0.05) +
#   coord_sf(crs = crs, xlim = c(-70.5, -65.5), ylim = c(39.5, 43.5)) +
#   scale_color_viridis_c() +
#   facet_wrap(~year)


spp_df <- stack(spp_list) %>%
  select(spp = values,
         label2 = ind) %>%
  mutate(label2 = sprintf('italic("  %s")', label2))

index_ctl_annual_list <- data.frame(path = list.files(path = here::here("analysis/vast_seasonal_index/"),
                                                      pattern = ".*index_ctl_annual.rds$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(spp = gsub(pattern = ".*vast_seasonal_index/(.*)_seasonal(.*)",  "\\1", path),
         dat = purrr::map(path, ~readRDS(file = .x))) %>%
  tidyr::unnest(cols = c(dat)) %>%
  left_join(spp_df) %>%
  mutate(label = if_else(year == max(year),
                  as.character(spp),
                  NA_character_))


p1 <- index_ctl_annual_list %>%
  filter(epu %in% c("Georges Bank", "Gulf of Maine")) %>%
  ggplot(aes(x = year, y = est, ymin = lower, ymax = upper, fill = spp, group = spp)) +
  geom_ribbon(alpha = 0.1) +
  geom_line(aes(color = spp)) +
  geom_point(aes(color = spp)) +
  labs(x = "Year",
       y = expression(Density~ln(abundance~km^{-2}))) +
  theme_minimal() +
  scale_x_continuous(limits = c(1994, 2024))+
  theme(legend.position = "none") +
  facet_wrap(~epu, ncol = 1)

p1 +
  geom_text_repel(aes(label = gsub("^.*$", " ", label)), # This will force the correct position of the link's right end.
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = "grey",
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 2,
                  # nudge_y = 10,
                  force = 2,
                  hjust = 0,
                  direction = "y",
                  na.rm = TRUE,
                  xlim = c(2017, Inf),
                  ylim = c(0, Inf),
                  parse = TRUE,
  ) +
  geom_text_repel(data = . %>% filter(!is.na(label)),
                  aes(label = paste0("  ", label), color = spp),
                  segment.alpha = 0, ## This will 'hide' the link
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.color = "grey",
                  box.padding = 0.1,
                  point.padding = 0.6,
                  nudge_x = 2,
                  # nudge_y = 1000,
                  force = 2,
                  hjust = 0,
                  direction = "y",
                  na.rm = TRUE,
                  xlim = c(2017, Inf),
                  ylim = c(0, Inf),
                  parse = TRUE)




ggplot(data = ln_index_ctl_annual) +
  geom_ribbon(aes(x = year, ymin = lower, ymax = upper, fill = EPU), alpha = 0.1) +
  geom_line(aes(x = year, y = est, color = EPU)) +
  labs(x = "Year",
       y = expression(Density~ln(abundance~km^{-2}))) +
  theme_bw() +
  facet_wrap(~EPU)

ggplot(data = ln_index_ctl) +
  geom_ribbon(aes(x = year_season, ymin = lower, ymax = upper, fill = EPU), alpha = 0.1) +
  geom_line(aes(x = year_season, y = Estimate, color = EPU)) +
  labs(x = "Year",
       y = expression(Density~ln(abundance~km^{-2}))) +
  theme_bw()



# Load Data -----
# data("ecomon_epu")
# set.seed(1234)
sp <- "calfin"

zoop_dat <- ecomon_epu %>%
  dplyr::filter(grepl(sprintf("^%s", sp) , spp),
                EPU %in% c("GB", "GOM", "MAB"),
                as.numeric(year) >= 1994) %>%
  dplyr::mutate(areaswept_km2 = 1,
                catch_ab = abundance) %>%
  # group_by(year, season) %>%
  # slice_sample(prop = .5) %>%
  droplevels() %>%
  data.frame()

# ggplot(zoop_dat, aes(x = lon, y = lat)) +
#   geom_point(data = zoop_dat %>% filter(abundance == 0), color = "black", fill = "black", shape = 21) +
#   geom_point(data = zoop_dat %>% filter(abundance > 0), aes(color = season, size = abundance), alpha = 0.5) +
#   facet_wrap(~year) +
#   labs(title = i) +
#   NULL


# Load data and quick exploration of structure
# Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
year_set = sort(unique(zoop_dat[,'year']))
season_set = c("winter", "spring", "summer", "fall")

# Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
yearseason_grid = expand.grid("season" = season_set, "year" = year_set)
yearseason_levels = apply(yearseason_grid[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_labels = round(yearseason_grid[,'year'] + (as.numeric(factor(yearseason_grid[,'season'], levels = season_set))-1)/length(season_set), digits=1)


plot( fit,
      projargs='+proj=natearth +lon_0=-68 +units=km',
      country = "united states of america",
      year_labels = yearseason_labels )



####
