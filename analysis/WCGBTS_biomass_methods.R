
# remotes::install_github("james-thorson/VAST", force = TRUE)
# remotes::install_github("lewisab/FishData", force = TRUE)
# remotes::install_github("james-thorson/FishData", force = TRUE)

# RootDir = "C:/Users/lewis.barnett/Documents/Postdoc/Eric postdoc/Distribution Shift PFMC/"
RootDir = "C:/Users/scott.large/Documents/projects/neusDFA/analysis/"

library(TMB)
library(VAST)
library(dplyr)

## Might need to make sure Rtools has path properly set
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")


# Global memory
Version = get_latest_version(package="VAST")  # Checks for latest version
Date = Sys.Date()

# Choose example
Example = "Ordination"

# Settings for each example
if( Example == "Ordination" ){
  Method = "Mesh"
  n_x = 500   # Specify number of stations (a.k.a. "knots")
  Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )
  ## Spatial (omega) and spatio-temporal (epsilon) factors used for each component
  ## Omega1 == zero inflation model and omega2 == the count-data model
  FieldConfig = c("Omega1"=2, "Epsilon1"=2, "Omega2"=0, "Epsilon2"=0)
  RhoConfig = c("Beta1"=3, "Beta2"=3, "Epsilon1"=4, "Epsilon2"=0)
  OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
  ObsModel = c(2,1)
  Options =  c("Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE)
  strata.limits <- data.frame('STRATA'="All_areas")
  Survey = "WCGBTS"
  Region = "California_Current"
  Species_set = c("Anoplopoma fimbria",
                  "Hydrolagus colliei", "Microstomus pacificus",
                  "Eopsetta jordani", "Sebastolobus alascanus",
                  "Parophrys vetulus", "Sebastolobus altivelis",
                  "Raja rhina")
  Vars2Correct = c()
}

# Create directory for results
RunDir = paste0(RootDir,Date, "/")
dir.create(RunDir)

# Save a record of settings
Record = list("Version"=Version,"Method"=Method,"grid_size_km"=50,"n_x"=n_x,"FieldConfig"=FieldConfig,"RhoConfig"=RhoConfig,"OverdispersionConfig"=OverdispersionConfig,"ObsModel"=ObsModel,"Kmeans_Config"=Kmeans_Config,"Region"=Region,"Survey"=Survey,"Species_set"=Species_set,"strata.limits"=strata.limits)
save( Record, file=file.path(RunDir,"Record.RData"))
capture.output( Record, file=paste0(RunDir,"Record.txt"))

# Download public data
DF = FishData::download_catch_rates(survey=Survey, species_set=Species_set, localdir=RunDir)

load(paste0(RunDir, "WCGBTS_download.RData"))
DF = Downloaded_data %>%
  rename(Sci = scientific_name,
         Year = `date_dim$year`,
         Wt = cpue_kg_per_ha_der,
         Lat = latitude_dd,
         Long = longitude_dd) %>%
  filter(!is.na(Wt))

#DF = FishData::download_catch_rates(survey=Survey, species_set=Species_set, localdir=RunDir)
# DF = droplevels(filter(DF, Sci != "Merluccius_productus", Sci != "Luidia_foliolata")) # remove hake from species list
Data_Geostat = data.frame( "spp"=DF[,"Sci"], "Year"=DF[,"Year"], "Catch_KG"=DF[,"Wt"], "AreaSwept_km2"=0.01, "Vessel"=0, "Lat"=DF[,"Lat"], "Lon"=DF[,"Long"] )
PredTF_i = rep(0, nrow(Data_Geostat))

# Create data frame of spatial domain
Extrapolation_List = make_extrapolation_info( Region=Region, strata.limits=strata.limits )

# Determine location of knots
Spatial_List = make_spatial_info( grid_size_km=50, n_x=n_x, Method=Method, Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], Extrapolation_List=Extrapolation_List, randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=RunDir, Save_Results=TRUE )

# Add knots to Data_Geostat
Data_Geostat = cbind( Data_Geostat, "knot_i"=Spatial_List$knot_i )

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
TmbList = Build_TMB_Fn("TmbData"=TmbData, "RunDir"=RunDir, "Version"= Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method )
Obj = TmbList[["Obj"]]

# Estimate parameters
Opt = TMBhelper::Optimize( obj=Obj, startpar=Obj$par+0.01*runif(length(Obj$par)), lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], getsd=TRUE, newtonsteps=1, savedir=RunDir, bias.correct=ifelse(length(Vars2Correct)>0,TRUE,FALSE), bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct=Vars2Correct) )

# Save results
Report = Obj$report()
Save = list("Opt"=Opt, "Report"=Report, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(RunDir,"Save.RData"))

# Get region-specific settings for plots
MapDetails_List = make_map_info( "Region"=Region, "NN_Extrap"=Spatial_List$PolygonList$NN_Extrap, "Extrapolation_List"=Extrapolation_List )
# Decide which years to plot
Year_Set = seq(min(Data_Geostat[,'Year']),max(Data_Geostat[,'Year']))
Years2Include = which( Year_Set %in% sort(unique(Data_Geostat[,'Year'])))

if( Example == "Ordination" ){
  Index = plot_biomass_index( cex.main = 0.75, oma = c(1,1.5,0,0.5), DirName=RunDir, TmbData=TmbData, Sdreport=Opt[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, strata_names=strata.limits[,1], use_biascorr=TRUE, category_names=levels(Data_Geostat[,'spp']) )
  plot_maps(plot_set=c(3), MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, Sdreport=Opt$SD, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=RunDir, Year_Set=Year_Set, Years2Include=Years2Include, Rotate=MapDetails_List[["Rotate"]], Cex=MapDetails_List[["Cex"]], Legend=MapDetails_List[["Legend"]], zone=MapDetails_List[["Zone"]], mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), cex=1.8, category_names=levels(Data_Geostat[,'spp']))
  Factor_list = Plot_factors( Report=Report, ParHat=Obj$env$parList(), Data=TmbData, RotationMethod="Varimax", SD=Opt$SD, mapdetails_list=MapDetails_List, Year_Set=Year_Set, category_names=levels(DF[,'Sci']), plotdir=RunDir )
  plot_anisotropy( FileName=paste0(RunDir,"Aniso.png"), Report=Report, TmbData=TmbData )

  # Names
  ThorsonUtilities::save_fig( paste0(RunDir,"Ordination_summary"), width=3, height=6 )
    par( mfrow=c(2,1), mar=c(0.5,1.5,2,0.5), oma=c(2,1,0,0), mgp=c(1.75,0.25,0), tck=-0.02 )
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
  dev.off()
  Codes = cbind( 1:nrow(Factor_list$Rotated_loadings$Epsilon1), rownames(Factor_list$Rotated_loadings$Epsilon1) )
  write.csv( Codes, file=paste0(RunDir,"Ordination_summary_codes.csv"), row.names=FALSE )

  # Omega
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

  # Epsilon
  Psi_rot = Factor_list$Rotated_factors[["Epsilon1"]]
  Psi_rot = ifelse( abs(Psi_rot)>4, sign(Psi_rot)*4, Psi_rot )
  zlim = c(-1,1) * max(abs(Psi_rot[1:TmbData$n_x,,]))
  Col = colorRampPalette(colors=c("darkblue","blue","lightblue","lightgreen","yellow","orange","red"),alpha=0.2)
  ThorsonUtilities::save_fig( paste0(RunDir,"Ordination_Epsilon"), width=MapDetails_List$MapSizeRatio['Width(in)']*4+0.5, height=MapDetails_List$MapSizeRatio['Height(in)']*2+1 )
    par( mfrow=c(2,4), mar=c(0.5,0.5,0.5,0.5), oma=c(2.5,3,2,1.5), mgp=c(1.75,0.25,0), tck=-0.02 )
    for( rowI in 1:2 ){
    for( colI in 1:4 ){
      tI = c(1,6,11,15)[colI]
      PlotMap_Fn( MappingDetails=list("worldHires",NULL), zlim=zlim, Mat=Psi_rot[,,tI][,rowI,drop=FALSE], PlotDF=MapDetails_List$PlotDF, MapSizeRatio=MapDetails_List$MapSizeRatio, Xlim=MapDetails_List$Xlim, Ylim=MapDetails_List$Ylim, FileName="", Year_Set="", Rescale=FALSE, Rotate=MapDetails_List$Rotate, Format="", Res=MapDetails_List$Res, zone=Extrapolation_List$zone, Cex=0.2, textmargin="", add=TRUE, pch=15, Legend=list("use"=FALSE), plot_legend_fig=FALSE )
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
}


