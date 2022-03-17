library(mgcv)
library(gratia)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)

data("ecomon_long")

sfc_as_cols <- function(x, names = c("x","y")) {
  stopifnot(inherits(x,"sf") && inherits(sf::st_geometry(x),"sfc_POINT"))
  ret <- sf::st_coordinates(x)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

calfin_dat <- ecomon_long %>%
  rename(count = abundance,
         taxa = spp) %>%
  mutate(jday = as.numeric(format(date, "%j")),
         year = as.numeric(format(date, "%Y")),
         occurence = ifelse(count == 0,
                            0,
                            1)) %>%
  filter(grepl("^calfin", taxa),
         year %in% c(1977:2017)) %>%
  na.omit(occurence)

# calfin_dat$year <- as.factor(calfin_dat$year)
# calfin_dat$occurence <- as.factor(calfin_dat$occurence)
## Ecomon delta model
## Pyr(lat, lon) = yr + s(lat, lon) + jday + jday^2

m1 <- gam(occurence ~ as.factor(year) + jday + I(jday^2) + s(lat, lon, bs = "ds"),
          data = calfin_dat,
          method = "REML",
          family = binomial(link = "logit"))

crs_epu <- 4269 # NAD83 https://epsg.org/crs_4269/NAD83.html
# crs_epu <- 9311 # NAD27 https://epsg.org/crs_9311/NAD27-US-National-Atlas-Equal-Area.html

epu <- ecodata::epu_sf %>%
  st_transform(crs = st_crs(crs_epu))
# st_crs(epu) <-
st_crs(epu) <- crs_epu


## Load EPU shapefiles as a spatialpolygonsdataframe
xmin = -76
xmax = -66
ymin = 36
ymax = 45
xlims <- c(xmin, xmax)
ylims <- c(ymin, ymax)


#Get base map and set CRS
coast <- rnaturalearth::ne_countries(scale = 10,
                      continent = "North America",
                      returnclass = "sf") %>%
  sf::st_transform(crs = st_crs(crs_epu))

ggplot()+
  # geom_sf(data = coast) +
  geom_sf(data = gb)

ggplot(gb, aes(x = long, y = lat)) +
  geom_path() +
  geom_point(data = fdepth, aes(x = os_x, y = os_y, colour = depth)) +
  coord_fixed() + ylab("Northing") + xlab("Easting") +
  scale_color_viridis()
ggplot()+
  geom_sf(data = coast) +
  geom_sf(data = epu) +
  coord_sf(xlim = xlims, ylim = ylims)

# appraise(m1)
# draw(m1)
# summary(m1)
#
# pred_dat <-expand.grid(
#   year = c(1977:2017),
#   jday = seq(from = quantile(as.numeric(calfin_dat$jday), probs = 0.2),
#              to = quantile(as.numeric(calfin_dat$jday), probs = 0.8), by = 30),
#   lat = seq(from = quantile(calfin_dat$lat, probs = c(0.2)),
#                    to = quantile(calfin_dat$lat, probs = c(0.8)), length = 20),
#   lon = seq(from = quantile(calfin_dat$lon, probs = c(0.2)),
#             to = quantile(calfin_dat$lon, probs = c(0.8)), length = 20))

point.sf <- epu %>%
  filter(EPU != "SS") %>%
  st_union() %>%
  st_as_sf() %>%
  st_transform(9311) %>%
  st_sample(size = 100, type = "regular") %>%
  st_transform(crs_epu)

# ggplot()+
#   geom_sf(data = coast) +
#   geom_sf(data = epu) +
#   geom_sf(data = point.sf) +
#   coord_sf(xlim = xlims, ylim = ylims)

pred_dat <- point.sf %>%
  st_cast("MULTIPOINT") %>%
  st_cast("POINT") %>%
  st_coordinates() %>%
  data.frame() %>%
  dplyr::select(lat = Y, lon = X) %>%
  crossing(year = as.factor(unique(calfin_dat$year)),
           jday = 120) %>%
  data.frame()

str(pred_dat)
str(calfin_dat)

cal_pred <- predict.gam(m1, newdata = pred_dat, type = "response", se.fit = TRUE)

tt <- cbind(pred_dat,
            fit = cal_pred$fit)

ggplot(tt, aes(x = lat, y = lon, color = fit)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~year)


calfin_count <- calfin_dat %>%
  filter(occurence != 0)

## Cyr(lat, lon) = yr + s(lat, lon) + jday
c1 <- gam(count ~ as.factor(year) + s(lat, lon, bs = "ds") + jday + I(jday^2),
          data = calfin_count,
          method = "REML",
          family = poisson(link = "log"))

# appraise(c1)
# draw(c1)
# summary(c1)

ggplot(cal_abundance, aes(x = lat, y = lon, color = count)) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~year)


cal_pred_count <- predict.gam(c1, newdata = pred_dat, type = "response", se.fit = TRUE)

cal_abundance <- cbind(pred_dat,
                       occurence = cal_pred$fit,
                       count = cal_pred_count$fit) %>%
  mutate(abundance = occurence * count,
         year = as.numeric(as.character(year))) %>%
  group_by(year) %>%
  summarize(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  right_join(data.frame(year = 1977:2017))


ggplot(cal_abundance, aes(x = year, y = mean_abundance)) +
  geom_point() +
  geom_line()





