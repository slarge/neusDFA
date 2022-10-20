library(TMB)
library(VAST)
library(dplyr)

## Might need to make sure Rtools has path properly set
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

Version = get_latest_version( package="VAST")

Method = c("Grid", "Mesh", "Spherical_mesh")[2]
grid_size_km = 50
n_x = c(50, 100, 250, 500, 1000, 2000)[1] # Number of stations
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )


FieldConfig = c(Omega1 = 3, Epsilon1 = 3,
                Omega2 = 3, Epsilon2 = 3)
RhoConfig = c(Beta1 = 0, Beta2 = 0,
              Epsilon1 = 0, Epsilon2 = 0)
OverdispersionConfig = c(Vessel = 0, VesselYear = 0)
ObsModel = c(2, 0)

Options = c(SD_site_density = 0, SD_site_logdensity = 0,
            Calculate_Range = 1, Calculate_evenness = 0,
            Calculate_effective_area = 1, Calculate_Cov_SE = 0,
            Calculate_Synchrony = 0, Calculate_Coherence = 0)


strata.limits <- data.frame(STRATA = "All_areas")

Region = "Eastern_Bering_Sea"
Species_set = c("Atheresthes stomias","Gadus chalcogrammus","Hippoglossoides elassodon")

Date = Sys.Date()

RootDir = "C:/Users/scott.large/Documents/projects/neusDFA/analysis/"
DateFile = paste0(RootDir,Date, "/")
dir.create(DateFile)

Record = list(Version = Version,
              Method = Method,
              grid_size_km = grid_size_km,
              n_x = n_x,
              FieldConfig = FieldConfig,
              RhoConfig = RhoConfig,
              OverdispersionConfig = OverdispersionConfig,
              ObsModel = ObsModel,
              Kmeans_Config = Kmeans_Config,
              Region = Region,
              Species_set = Species_set,
              strata.limits = strata.limits)
save(Record, file = file.path(DateFile, "Record.RData"))
capture.output(Record, file = paste0(DateFile, "Record.txt"))

DF = FishData::download_catch_rates(survey = Region,
                                    species_set = Species_set)

Data_Geostat = data.frame(spp = DF[, "Sci"],
                          Year = DF[, "Year"],
                          Catch_KG = DF[, "Wt"],
                          AreaSwept_km2 = 0.01,
                          Vessel = 0,
                          Lat = DF[, "Lat"],
                          Lon = DF[, "Long"])

Data_Geostat = cbind(Data_Geostat, knot_i = Spatial_List$knot_i)

Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )

# Determine location of knots
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
                                  DirPath=DateFile)


TmbData = Data_Fn(Version = Version,
                  FieldConfig = FieldConfig,
                  OverdispersionConfig = OverdispersionConfig,
                  RhoConfig = RhoConfig,
                  ObsModel = ObsModel,
                  c_i = as.numeric(Data_Geostat[,"spp"]) - 1,
                  b_i = Data_Geostat[, "Catch_KG"],
                  a_i = Data_Geostat[, "AreaSwept_km2"],
                  v_i = as.numeric(Data_Geostat[,"Vessel"]) - 1,
                  s_i = Data_Geostat[, "knot_i"] - 1,
                  t_i = Data_Geostat[, "Year"],
                  a_xl = Spatial_List$a_xl,
                  MeshList = Spatial_List$MeshList,
                  GridList = Spatial_List$GridList,
                  Method = Spatial_List$Method,
                  Options = Options)


TmbList = Build_TMB_Fn(TmbData = TmbData,
                       RunDir = DateFile,
                       Version = Version,
                       RhoConfig = RhoConfig,
                       loc_x = Spatial_List$loc_x,
                       Method = Method)
Obj = TmbList[["Obj"]]

# Opt = TMBhelper::Optimize(obj = Obj,
#                           lower = TmbList[["Lower"]],
#                           upper = TmbList[["Upper"]],
#                           getsd = TRUE,
#                           newtonsteps = 1,
#                           savedir = DateFile,
#                           bias.correct = FALSE)

Opt = TMBhelper::Optimize( obj=Obj,
                           lower=TmbList[["Lower"]],
                           upper=TmbList[["Upper"]],
                           getsd=TRUE,
                           savedir=DateFile,
                           bias.correct=TRUE,
                           bias.correct.control=list(sd=FALSE,
                                                     split=NULL,
                                                     nsplit=1,
                                                     vars_to_correct="Index_cyl"),
                           newtonsteps=1 )



Report = Obj$report()
Save = list("Opt"=Opt,
            "Report"=Report,
            "ParHat"=Obj$env$parList(Opt$par),
            "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))

plot_data(Extrapolation_List=Extrapolation_List, Spatial_List=Spatial_List, Data_Geostat=Data_Geostat, PlotDir=DateFile )

Enc_prob = plot_encounter_diagnostic( Report=Report, Data_Geostat=Data_Geostat, DirName=DateFile)


Q = plot_quantile_diagnostic( TmbData=TmbData, Report=Report, FileName_PP="Posterior_Predictive",
                              FileName_Phist="Posterior_Predictive-Histogram",
                              FileName_QQ="Q-Q_plot", FileName_Qhist="Q-Q_hist", DateFile=DateFile)

# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))


plot_residuals(Lat_i=Data_Geostat[,'Lat'],
               Lon_i=Data_Geostat[,'Lon'],
               TmbData=TmbData, Report=Report,
               Q=Q,
               savedir=DateFile,
               MappingDetails=MapDetails_List[["MappingDetails"]],
               PlotDF=MapDetails_List[["PlotDF"]],
               MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
               Xlim=MapDetails_List[["Xlim"]],
               Ylim=MapDetails_List[["Ylim"]],
               FileName=DateFile,
               Year_Set=Year_Set,
               Years2Include=Years2Include,
               Rotate=MapDetails_List[["Rotate"]],
               Cex=MapDetails_List[["Cex"]],
               Legend=MapDetails_List[["Legend"]],
               zone=MapDetails_List[["Zone"]],
               mar=c(0,0,2,0),
               oma=c(3.5,3.5,0,0),
               cex=1.8)


Opt$AIC

plot_anisotropy( FileName=paste0(DateFile,"Aniso.png"), Report=Report, TmbData=TmbData )


Cov_List = Summarize_Covariance( Report=Report,
                                 ParHat=Obj$env$parList(),
                                 Data=TmbData, SD=Opt$SD,
                                 plot_cor=FALSE,
                                 category_names=levels(Data_Geostat[,'spp']),
                                 plotdir=DateFile,
                                 plotTF=FieldConfig,
                                 mgp=c(2,0.5,0),
                                 tck=-0.02,
                                 oma=c(0,5,2,2) )


plot_maps(plot_set=c(3),
          MappingDetails=MapDetails_List[["MappingDetails"]],
          Report=Report,
          Sdreport=Opt$SD,
          PlotDF=MapDetails_List[["PlotDF"]],
          MapSizeRatio=MapDetails_List[["MapSizeRatio"]],
          Xlim=MapDetails_List[["Xlim"]],
          Ylim=MapDetails_List[["Ylim"]],
          FileName=DateFile,
          Year_Set=Year_Set,
          Years2Include=Years2Include,
          Rotate=MapDetails_List[["Rotate"]],
          Cex=MapDetails_List[["Cex"]],
          Legend=MapDetails_List[["Legend"]],
          zone=MapDetails_List[["Zone"]],
          mar=c(0,0,2,0),
          oma=c(3.5,3.5,0,0),
          cex=1.8,
          category_names=levels(Data_Geostat[,'spp']))


Index = plot_biomass_index( DirName=DateFile,
                            TmbData=TmbData,
                            Sdreport=Opt[["SD"]],
                            Year_Set=Year_Set,
                            Years2Include=Years2Include,
                            strata_names=strata.limits[,1],
                            use_biascorr=TRUE,
                            category_names=levels(Data_Geostat[,'spp']) )

pander::pandoc.table( Index$Table[,c("Category","Year","Estimate_metric_tons","SD_mt")] )



plot_range_index(Report=Report, TmbData=TmbData, Sdreport=Opt[["SD"]], Znames=colnames(TmbData$Z_xm), PlotDir=DateFile, category_names=levels(Data_Geostat[,'spp']), Year_Set=Year_Set)


Plot_Overdispersion( filename1=paste0(DateDir,"Overdispersion"), filename2=paste0(DateDir,"Overdispersion--panel"), Data=TmbData, ParHat=ParHat, Report=Report, ControlList1=list("Width"=5, "Height"=10, "Res"=200, "Units"='in'), ControlList2=list("Width"=TmbData$n_c, "Height"=TmbData$n_c, "Res"=200, "Units"='in') )


Plot_factors( Report=Report, ParHat=Obj$env$parList(), Data=TmbData, SD=Opt$SD, mapdetails_list=MapDetails_List, Year_Set=Year_Set, category_names=levels(DF[,'Sci']), plotdir=DateFile )
