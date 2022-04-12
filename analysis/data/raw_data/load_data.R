library(dplyr)
library(tidyr)
library(marmap)

## Right now the most up to date survdat is on ECSA. Soon this will migrate to ecodata
load(url("https://github.com/NOAA-EDAB/ECSA/blob/master/data/Survdat.RData?raw=true"))
usethis::use_data(survdat, overwrite = TRUE)

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Download and reformat EcoMon data ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# download.file(url = "https://www.nodc.noaa.gov/archive/arc0143/0187513/1.1/data/0-data/EcoMon_Plankton_Data_v3_5.csv",
#               destfile = "analysis/data/raw_data/EcoMon_Plankton_data_v3_5.csv")

ecomon_long <- readr::read_csv("analysis/data/raw_data/EcoMon_Plankton_Data_v3_5.csv") %>%
  mutate(volume = as.double(volume_100m3),
         date = as.Date(date, format="%d-%b-%y")) %>%
  select(-ends_with("_10m2"),
         -ends_with("_10m2x"),
         -starts_with("volume_"),
         cavoli_100m3 = cavoli_100m3x) %>%
  pivot_longer(cols = c(ends_with("_100m3"), "clauso"),
               names_to = "spp",
               values_to = "abundance")

usethis::use_data(ecomon_long, overwrite = TRUE)


# These oblique tows include a measure of total volume swept, and we divide the total number of
# zoop by volume swept and then multiply by the seafloor depth at the beginning of the tow to
# obtain vertically integrated numbers-density.
# Using vertically integrated numbers-density as response variable then allows us to predict
# vertically integrated densities across a standard survey area, where the sum across this survey
# area represents a prediction of vertical and spatially integrated abundance in numbers.


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
## Prepare EcoMon data for VAST ##
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

# download.file(url = "https://www.nodc.noaa.gov/archive/arc0143/0187513/1.1/data/0-data/EcoMon_Plankton_Data_v3_5.csv",
#               destfile = "analysis/data/raw_data/EcoMon_Plankton_data_v3_5.csv")

data("ecomon_long")
#

# prepare bottom depth
xlims <- c(-77, -65)
ylims <- c(35, 45)
res <- 1
bath_filename <- sprintf("marmap_coord_%s;%s;%s;%s_res_%s.csv",
                         xlims[1], ylims[1], xlims[2], ylims[2], res)

if(!bath_filename %in% list.files("analysis/data/raw_data")){
  nesbath <- marmap::getNOAA.bathy(lon1 = xlims[1], lon2 = xlims[2],
                                   lat1 = ylims[1], lat2 = ylims[2],
                                   resolution = res,
                                   keep = TRUE) %>%
    marmap::as.raster()

  file.copy(bath_filename, "analysis/data/raw_data")
  file.remove(bath_filename)
} else {
  nesbath <- marmap::read.bathy(sprintf("analysis/data/raw_data/%s", bath_filename), header = T) %>%
    marmap::as.raster()
}

ecomon_format <- ecomon_long %>%
  mutate(marmap_depth = raster::extract(nesbath, y = cbind(.$lon, .$lat)) * -1,
         depth = ifelse(depth == 9999,
                        marmap_depth,
                        depth),
         spp = gsub("_100m3", "", spp),
         date = as.Date(date, format = "%d-%b-%y"),
         id = as.factor(paste0(cruise_name, "_", station)),
         vessel = as.factor(substr(cruise_name, start = 1, stop = 2)),
         areaswept_km2 = 1,
         day = as.numeric(strftime(date, format = "%j")),
         year = as.numeric(strftime(date, format = "%Y")),
         season = case_when(day %in% 1:90 ~ "winter",
                            day %in% 91:181 ~ "spring",
                            day %in% 182:273 ~ "summer",
                            day %in% 274:365 ~ "fall",
                            TRUE ~ NA_character_)) %>%
  dplyr::filter(grepl("6B3", zoo_gear),
                !grepl("6B5", ich_gear))

crs_epu <- 4269 # NAD83 https://epsg.org/crs_4269/NAD83.html
# crs_epu <- 9311 # NAD27 https://epsg.org/crs_9311/NAD27-US-National-Atlas-Equal-Area.html

epu <- ecodata::epu_sf %>%
  sf::st_transform(crs = sf::st_crs(crs_epu))

sf::sf_use_s2(FALSE)

## Post stratify data according to EPUs
ecomon_epu <- ecomon_format %>%
  sf::st_as_sf(coords = c("lon","lat"), crs = crs_epu) %>%
  sf::st_join(ecodata::epu_sf) %>%
  sfc_as_cols(names = c("lon", "lat")) %>%
  sf::st_drop_geometry() %>%
  select(id,
         spp,
         EPU,
         day,
         season,
         bottom_depth = depth,
         vessel,
         abundance,
         areaswept_km2,
         lat,
         lon,
         year,
         zoo_gear,
         ich_gear) %>%
  data.frame()

usethis::use_data(ecomon_epu, overwrite = TRUE)

