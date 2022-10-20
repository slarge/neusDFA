# remotes command to get TMB from GitHub
# remotes::install_github("kaskr/adcomp/TMB", force = TRUE)
## https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Steps-to-install-TMB
# install.packages("INLA", repos = c(getOption("repos"),
#                                    INLA = "https://inla.r-inla-download.org/R/stable"),
#                  dep = TRUE)
# remotes::install_github("james-thorson/VAST")

library(TMB)
library(VAST)
library(dplyr)
# library(sf)
# library(INLA)
library(ggplot2)

## For some reason I need to make sure Rtools has path properly set, else TMB won't compile
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

##############
## Get data ##
##############
# Date = Sys.Date()
Date = "2019-02-22"
## Change to your target directory
# RootDir = "~/"
# RootDir = "C:/Users/scott.large/Documents/projects/neusDFA/analysis/guts/"
RootDir = paste0(getwd(), "/projects/neusDFA/analysis/guts/")
RunDir = paste0(RootDir,Date, "/")
dir.create(RunDir)

dat_url <- "https://github.com/NOAA-EDAB/ECSA/blob/master/data/allfhsg.RData?raw=true"
# Windows
# load(url(survdat_url))
# other OS's
download.file(dat_url, paste0(RunDir, "allfhsg.Rdata"), mode = "wb")
load(paste0(RunDir, "allfhsg.Rdata"))

sp_list <- c(73)
# year_range <- 1998:2017

## General mapping parameters
xmin = -77
xmax = -65
ymin = 35
ymax = 45

xlims <- c(xmin, xmax)
ylims <- c(ymin, ymax)
crs <-  "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Download data layers
## 1) Strata
epu = c("Georges_Bank", "Gulf_of_Maine")

## 2) North America layer
ne_countries <- rnaturalearth::ne_countries(scale = 10,
                                            continent = "North America",
                                            returnclass = "sf") %>%
  sf::st_transform(crs = crs)

## 3) State layer
ne_states <- rnaturalearth::ne_states(country = "united states of america",
                                      returnclass = "sf") %>%
  sf::st_transform(crs = crs)

## 4) select epu to work with
strata_grid <- sf::st_read("analysis/data/raw_data/EPU_NOESTUARIES.shp",
                           quiet = TRUE) %>%
  mutate(EPU = dplyr::case_when(EPU == "GB" ~ "Georges_Bank",
                                EPU == "MAB" ~ "Mid_Atlantic_Bight",
                                EPU == "GOM" ~ "Gulf_of_Maine",
                                EPU == "SS" ~ "Scotian_Shelf",
                                TRUE ~ "Other")) %>%
  filter(EPU %in% epu)


## 5) Select and rename survey data
dat_df <- allfhsg %>%
  dplyr::filter(svspp %in% sp_list,
         # year %in% year_range,
         season == "SPRING",
         pynam != 'BLOWN',
         pynam != 'PRESERVED',
         pynam != ' ',
         pynam != 'EMPTY',
         purcode == 10) %>%
  dplyr::select(pynam,
         Year = year,
         Catch_KG = pyamtv,
         Lat = declat,
         Lon = declon) %>%
  dplyr::distinct(.keep_all = TRUE)


## Select the most frequently observed prey
# Link and Garrison (2002a)
# Fogarty et al. 1991
other_hakes <- c("UROPHYCIS TENUIS", "UROPHYCIS CHUSS", "UROPHYCIS SP", "UROPHYCIS CHESTERI", "UROPHYCIS REGIA")
silver_hake <- c("MERLUCCIUS BILINEARIS")
redfish <- "SEBASTES FASCIATUS"
cancer_spp <- c("CANCER BOREALIS", "CANCRIDAE", "CANCER IRRORATUS", "CANCER SP")
atl_herring <- c("CLUPEA HARENGUS", "CLUPEIDAE")
atl_mackerel <- c("SCOMBER SCOMBRUS", "SCOMBRIDAE")
crangon <- c("CRANGON SEPTEMSPINOSA", "CRANGONIDAE")
euphausids <- "EUPHAUSIIDAE"
sandlance <- "AMMODYTES SP"
amphipod <- c("AMPHIPODA", "CAPRELLIDAE")
gadids <- c("GADIDAE", "GADIFORMES", "GADUS MORHUA")
cephalopods <- c("CEPHALOPODA", "OCTOPODA", "ILLEX SP", "LOLIGO SP", "SEPIOLIDAE", "LOLIGO PEALEII", "BATHYPOLYPUS ARCTICUS")
pout <- "MACROZOARCES AMERICANUS"
shrimp <- c("CRUSTACEA SHRIMP", "CRANGON SEPTEMSPINOSA", "DECAPODA SHRIMP", "GAMMARIDEA")

diets <- dat_df %>%
  select(pynam) %>%
  distinct(.keep_all = TRUE) %>%
  mutate(spp = case_when(pynam %in% other_hakes~ "other_hakes",
                         pynam %in% silver_hake ~ "silver_hake",
                         pynam %in% cancer_spp ~ "cancer_spp",
                         pynam %in% atl_herring ~ "atl_herring",
                         pynam %in% atl_mackerel ~ "atl_mackerel",
                         pynam %in% crangon ~ "crangon",
                         pynam %in% redfish ~ "redfish",
                         pynam %in% euphausids ~ "euphausids",
                         pynam %in% sandlance ~ "sandlance",
                         pynam %in% amphipod ~ "amphipod",
                         pynam %in% gadids ~ "gadids",
                         pynam %in% cephalopods ~ "cephalopods",
                         pynam %in% pout ~ "ocean_pout",
                         pynam %in% shrimp ~ "shrimp",
                         TRUE ~ NA_character_)) %>%
  filter(!is.na(spp))

#
# td <- diets %>%
#   group_by(pynam, spp, analcom3, pycomnam2) %>%
#   summarize(count = n()) %>%
#   arrange(-count)
#
#

gut_df <- diets %>%
  left_join(dat_df)  %>%
  # dplyr::filter(pynam %in% unique(diets$spp),
  #               !is.na(spp)) %>%
  dplyr::mutate(AreaSwept_km2 = 0.01,
                Lon = -Lon,
                Vessel = 0,
                spp = as.factor(gsub("\\s", "_", spp)),
                Catch_KG = log(Catch_KG + 1)) %>%
  select(-pynam)

## 6) Select survey data within EPU polygon
gut_sf <- sf::st_as_sf(gut_df, coords = c("Lon", "Lat"), crs = crs)
strata_grid <- sf::st_as_sf(strata_grid, crs = crs)

gut_grid <- sf::st_join(gut_sf, strata_grid, st_within) %>%
  dplyr::filter(!is.na(EPU))

# add coordinates to dataframe
gut_grid$Lon <- sf::st_coordinates(gut_grid)[,1] # get coordinates
gut_grid$Lat <- sf::st_coordinates(gut_grid)[,2] # get coordinates

# test plot to make sure everything looks correct
epu_sp_plot <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = ne_countries, color = "grey60", size = 0.25) +
  ggplot2::geom_sf(data = ne_states, color = "grey60", size = 0.05) +
  ggplot2::geom_point(data = gut_grid, ggplot2::aes(fill = spp, size = Catch_KG, x = Lon, y = Lat), shape = 21, alpha = 0.5) +
  ggplot2::geom_sf(data = strata_grid, size = 0.05, alpha = .1, color = "grey40") +
  viridis::scale_fill_viridis(discrete = TRUE) +
  ggplot2::coord_sf(crs = crs, xlim = xlims, ylim = ylims) +
  ggthemes::theme_map()
epu_sp_plot

## 7) Create Data_Geostat df
Data_Geostat <- gut_grid %>%
  select(spp,
         Year,
         Catch_KG,
         Lat,
         Lon,
         AreaSwept_km2,
         Vessel) %>%
  st_set_geometry((NULL))

## 8) Start VAST prep
Region = "northwest_atlantic"
strata.limits = "EPU"
Version = get_latest_version( package="VAST")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Build objects related to spatial information ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#"Grid" is a 2D AR1 process, and "Mesh" is the SPDE (stochastic partial differential equation) method with geometric anisotropy
Method = c("Grid", "Mesh", "Spherical_mesh")[2]

# determines spatial resolution when Method="Grid"
grid_size_km = c(10, 50, 100)[1]
# Specify number of stations (a.k.a. "knots")
n_x = c(250, 500, 1000, 1100)[3]
# Parameters for k-means
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )

# Create data frame of spatial domain
## SL: I customized the function from FishStatUtils to select EPUs
NWA_Extrapolation_Data_Fn <- function (strata.limits = NULL, zone = NA, EPU = NA, ...) {
  if (is.null(strata.limits)) {
    strata.limits = list(All_areas = 1:1e+05)
  }
  message("Using strata ", strata.limits)
  utils::data(northwest_atlantic_grid, package = "FishStatsUtils")
  Data_Extrap <- northwest_atlantic_grid
  if(!is.na(EPU) && EPU %in% unique(Data_Extrap$EPU)){
    Data_Extrap <- Data_Extrap[Data_Extrap$EPU %in% EPU,]
    Data_Extrap <- droplevels(Data_Extrap)
  }
  Area_km2_x = Data_Extrap[, "Area_in_survey_km2"]
  Tmp = cbind(BEST_DEPTH_M = 0, BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
  if (length(strata.limits) == 1 && strata.limits[1] == "EPU") {
    a_el = matrix(NA, nrow = nrow(Data_Extrap),
                  ncol = length(unique(Data_Extrap[, "EPU"])),
                  dimnames = list(NULL, unique(Data_Extrap[, "EPU"])))
    for (l in 1:ncol(a_el)) {
      a_el[, l] = ifelse(Data_Extrap[, "EPU"] == unique(Data_Extrap[, "EPU"])[l], Area_km2_x, 0)
    }
  }
  else {
    a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap),
                                ncol = length(strata.limits), dimnames = list(NULL,
                                                                              names(strata.limits))))
    for (l in 1:ncol(a_el)) {
      a_el[, l] = ifelse(Data_Extrap[, "stratum_number"] %in%
                           strata.limits[[l]], Area_km2_x, 0)
    }
  }
  tmpUTM = Convert_LL_to_UTM_Fn(Lon = Data_Extrap[, "Lon"],
                                Lat = Data_Extrap[, "Lat"], zone = zone)
  Data_Extrap = cbind(Data_Extrap, Include = 1)
  Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
  Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM,
                                                                    "zone"), flip_around_dateline = FALSE, Area_km2_x = Area_km2_x)
  return(Return)
}

Extrapolation_List <- NWA_Extrapolation_Data_Fn(strata.limits = strata.limits, EPU = epu)

# spatio-temporal parameter estimation
Spatial_List = make_spatial_info( grid_size_km=grid_size_km,
                                  n_x=n_x,
                                  Method=Method,
                                  Lon=Data_Geostat[,'Lon'],
                                  Lat=Data_Geostat[,'Lat'],
                                  Extrapolation_List=Extrapolation_List,
                                  randomseed=Kmeans_Config[["randomseed"]],
                                  nstart=Kmeans_Config[["nstart"]],
                                  iter.max=Kmeans_Config[["iter.max"]],
                                  Save_Results=TRUE,
                                  DirPath=RunDir)

# Add knots to Data_Geostat
# Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )
Data_Geostat <- readRDS(paste0(RootDir, "/", Date, "/vast_cod_guts.rds"))
Spatial_List <- readRDS(paste0(RootDir, "/", Date, "/Spatial_List.rds"))

## Spatial (omega) and spatio-temporal (epsilon) factors used for each component
## Omega1 == encounter probability and omega2 == positive catch rates
## 0 = off, "AR1" = AR1 process, and >0 is the number of elements in a factor-analysis covariance
FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=0, "Epsilon2"=0)
## Intercepts (Beta1 and Beta2) or spatio-temporal variation (Epsilon1 and Epsilon2)
## is structured among time intervals (0: each year as fixed effect; 1: each year as random following IID distribution;
## 2: each year as random following a random walk; 3: constant among years as fixed effect; 4: each year as random following AR1 process)
RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=0)

OverdispersionConfig = c("Eta1"=0, "Eta2"=0)

## I'm not 100% certain what I should put here...
ObsModel = c(2,1)
##

Options =  c("Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE)

PredTF_i = rep(0, nrow(Data_Geostat))
Vars2Correct = c()


# save.image(file = "neus_vast.rda")

# Assemble data
TmbData <- Data_Fn("Version" = Version,
                   "FieldConfig" = FieldConfig,
                   "OverdispersionConfig" = OverdispersionConfig,
                   "RhoConfig" = RhoConfig,
                   "ObsModel" = ObsModel,
                   "c_i" = as.numeric(Data_Geostat[,'spp'])-1,
                   "b_i" = Data_Geostat[,'Catch_KG'],
                   "a_i" = Data_Geostat[,'AreaSwept_km2'],
                   "v_i" = as.numeric(Data_Geostat[,'Vessel'])-1,
                   "s_i" = Data_Geostat[,'knot_i']-1,
                   "t_i" = Data_Geostat[,'Year'],
                   "a_xl" = Spatial_List$a_xl,
                   "MeshList" = Spatial_List$MeshList,
                   "GridList" = Spatial_List$GridList,
                   "Method" = Spatial_List$Method,
                   "Options" = Options,
                   "PredTF_i" = PredTF_i)


# Build model
# dyn.unload( paste0(RunDir,"/",TMB::dynlib(Version)) )
TmbList = Build_TMB_Fn("TmbData"=TmbData,
                       "Version"= Version,
                       "RhoConfig"=RhoConfig,
                       "loc_x"=Spatial_List$loc_x,
                       "Method"=Method,
                       "RunDir"=RunDir)
Obj = TmbList[["Obj"]]

# Estimate parameters
Opt = TMBhelper::Optimize( obj=Obj,
                           startpar=Obj$par+0.01*runif(length(Obj$par)),
                           lower=TmbList[["Lower"]],
                           upper=TmbList[["Upper"]],
                           getsd=TRUE,
                           newtonsteps=1,
                           savedir=RunDir,
                           bias.correct=ifelse(length(Vars2Correct)>0,TRUE,FALSE),
                           bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1,
                                                     vars_to_correct=Vars2Correct) )


# Save results
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(RunDir,"Save.RData"))

## Plotting functions from Lewis
# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))


Index = plot_biomass_index( cex.main = 0.75, oma = c(1,1.5,0,0.5), DirName=RunDir, TmbData=TmbData, Sdreport=Opt[["SD"]],
                            Year_Set=Year_Set, Years2Include=Years2Include, strata_names=strata.limits,
                            use_biascorr=TRUE, category_names=levels(Data_Geostat[,'spp']) )

plot_maps(plot_set=c(3),
          MappingDetails=MapDetails_List[["MappingDetails"]],
          Report=Report,
          Sdreport=Opt$SD,
          PlotDF=MapDetails_List[["PlotDF"]],
          MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
          Xlim=MapDetails_List[["Xlim"]],
          Ylim=MapDetails_List[["Ylim"]],
          FileName=RunDir,
          Year_Set=Year_Set,
          Years2Include=Years2Include,
          Rotate=MapDetails_List[["Rotate"]],
          Cex=MapDetails_List[["Cex"]],
          Legend=MapDetails_List[["Legend"]],
          zone=MapDetails_List[["Zone"]],
          mar=c(0,0,2,0), oma=c(3.5,3.5,0,0),
          cex=1.8, category_names=levels(Data_Geostat[,'spp']))


Factor_list = Plot_factors( Report=Report, ParHat=Obj$env$parList(),
                            Data=TmbData, RotationMethod="Varimax",
                            SD=Opt$SD, mapdetails_list=MapDetails_List,
                            Year_Set=Year_Set, category_names=levels(Data_Geostat[,'spp']), plotdir=RunDir )
plot_anisotropy( FileName=paste0(RunDir,"Aniso.png"), Report=Report, TmbData=TmbData )

Range = max(abs( c(as.vector(Factor_list$Rotated_loadings$Omega1), as.vector(Factor_list$Rotated_loadings$Epsilon1)) ))
# Spatial
plot( 1, type="n", xlim=c(-1,1)*Range, ylim=c(-1,1)*Range, xaxt="n", main="Spatial", xlab="", ylab="")
text( x=Factor_list$Rotated_loadings$Omega1[,1], y=Factor_list$Rotated_loadings$Omega1[,2], labels=1:nrow(Factor_list$Rotated_loadings$Omega1) )
abline( h=0, lty="dotted" )
abline( v=0, lty="dotted" )
# Spatio-temporal
plot( 1, type="n", xlim=c(-1,1)*Range, ylim=c(-1,1)*Range, xaxt="n", main="Spatio-temporal", xlab="", ylab="")
text( x=Factor_list$Rotated_loadings$Epsilon1[,1], y=Factor_list$Rotated_loadings$Epsilon1[,2], labels=1:nrow(Factor_list$Rotated_loadings$Epsilon1) )
abline( h=0, lty="dotted" )
abline( v=0, lty="dotted" )
axis(1)
mtext( side=1:2, text=c("Factor 1", "Factor 2"), outer=TRUE, line=c(1,0) )


Psi_rot = Factor_list$Rotated_factors[["Omega1"]]
Psi_rot = ifelse( abs(Psi_rot)>4, sign(Psi_rot)*4, Psi_rot )
zlim = c(-1,1) * max(abs(Psi_rot[1:TmbData$n_x,,]))
Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"),alpha=0.2)
ThorsonUtilities::save_fig( paste0(RunDir,"Ordination_Omega"), width=MapDetails_List$MapSizeRatio['Width(in)']*2+0.5, height=MapDetails_List$MapSizeRatio['Height(in)']*1+1 )
par( mfcol=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(2.5,4.5,2,0), mgp=c(1.75,0.25,0), tck=-0.02 )
for( colI in 1:2 ){
  PlotMap_Fn( MappingDetails=list("worldHires",NULL), zlim=zlim, Mat=Psi_rot[,,1][,colI,drop=FALSE], PlotDF=MapDetails_List$PlotDF, MapSizeRatio=MapDetails_List$MapSizeRatio, Xlim=MapDetails_List$Xlim, Ylim=MapDetails_List$Ylim, FileName="", Year_Set="", Rescale=FALSE, Rotate=MapDetails_List$Rotate, Format="", Res=MapDetails_List$Res, zone=Extrapolation_List$zone, Cex=0.15, textmargin="", add=TRUE, pch=15, Legend=list("use"=FALSE), plot_legend_fig=FALSE )
  mtext( side=3, text=paste0("Factor ",colI), line=0.5, font=2 )
  if( colI==1 ) axis(2)
  axis(1)
  if( colI==2 ){
    FishStatsUtils:::smallPlot( FishStatsUtils:::Heatmap_Legend(colvec=Col(50), heatrange=zlim, dopar=FALSE), x=MapDetails_List$Legend$x, y=MapDetails_List$Legend$y, mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.2, font=2 )  #
  }
}
mtext( side=1:2, text=c("Longitude","Latitude"), outer=TRUE, line=c(1,3) )
dev.off()


Psi_rot = Factor_list$Rotated_factors[["Epsilon1"]]
Psi_rot = ifelse( abs(Psi_rot)>4, sign(Psi_rot)*4, Psi_rot )
zlim = c(-1,1) * max(abs(Psi_rot[1:TmbData$n_x,,]))
Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"),alpha=0.2)
ThorsonUtilities::save_fig( paste0(RunDir,"Ordination_Epsilon"), width=MapDetails_List$MapSizeRatio['Width(in)']*4+0.5, height=MapDetails_List$MapSizeRatio['Height(in)']*2+1 )
par( mfrow=c(2,4), mar=c(0.5,0.5,0.5,0.5), oma=c(2.5,3,2,1.5), mgp=c(1.75,0.25,0), tck=-0.02 )
for( rowI in 1:2 ){
  for( colI in 1:4 ){
    tI = c(1,6,11,20)[colI]
    PlotMap_Fn( MappingDetails=list("worldHires",NULL), zlim=zlim, Mat=Psi_rot[,,tI][,rowI,drop=FALSE],
                PlotDF=MapDetails_List$PlotDF,
                MapSizeRatio=MapDetails_List$MapSizeRatio,
                Xlim=MapDetails_List$Xlim,
                Ylim=MapDetails_List$Ylim,
                FileName="",
                Year_Set="",
                Rescale=FALSE,
                Rotate=MapDetails_List$Rotate,
                Format="",
                Res=MapDetails_List$Res,
                zone=Extrapolation_List$zone, Cex=0.2,
                textmargin="", add=TRUE, pch=15, Legend=list("use"=FALSE), plot_legend_fig=FALSE )
    if( colI==4 ) mtext( side=4, text=paste0("Factor ",rowI), line=0.5 )
    if( colI==1 ) axis(2)
    if( rowI==2 ) axis(1)
    if( rowI==1 ) mtext( side=3, text=paste0("Year ",Year_Set[tI]), line=0.5, font=2 )
    if( colI==4 & rowI==2 ){
      FishStatsUtils:::smallPlot( FishStatsUtils:::Heatmap_Legend(colvec=Col(50), heatrange=zlim, dopar=FALSE), x=MapDetails_List$Legend$x, y=MapDetails_List$Legend$y, mar=c(0,0,0,0), mgp=c(2,0.5,0), tck=-0.2, font=2 )  #
    }
  }}
mtext( side=1:2, text=c("Longitude","Latitude"), outer=TRUE, line=c(1,1) )
dev.off()





