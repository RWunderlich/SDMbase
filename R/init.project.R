# SDMbase, init.project.R
# Requiring external packages
require(raster)
require(sp)
require(SDMTools)
require(rgdal)


### Basic module to set-up a project, download, import, align and resample layers including those required for a projection
# At the moment only pre-prepared asc files are read in, next version will be able to download the whole dataset from the worldclim website.


# Setting directories, file types and extent
writeLines("Please enter path of maxent.jar in the format /home/USER/Desktop/Maxent/ press-ENTER")
pathmaxent <- scan(what = "character", nmax = 1)

writeLines("Please enter full path for the current project in the format /home/USER/Desktop/example/ press-ENTER")
pathproject <- scan(what = "character", nmax = 1)
dir.create(path = pathproject)
setwd(dir = pathproject)
dir.create(path = "./layers/")
dir.create(path = "./locations/")
writeLines("Does the project include a projection in time or space? 2=Space, 1=Time, 0=No press-ENTER")
projectprojections <- scan(what = "integer", nmax = 1)
if (projectprojections == 1)
  {
  dir.create(path = "./projection/")
  writeLines("Please enter path of the projection layers in the format /home/USER/Desktop/Projectionlayers/ press-ENTER")
  pathprojections <- scan(what = "character", nmax = 1)
  }

if (projectprojections == 2)
  {
  dir.create(path = "./projection/")
  }

writeLines("Please enter path of the global bioclim[1:19], elevation and other.asc files in the format /home/USER/Desktop/Layers/ press-ENTER")
pathlayers <- scan(what = "character", nmax = 1)
if (projectprojections == 2) {pathprojections <- pathlayers}

writeLines("Please enter xmin, xmax, ymin and ymax of the project extent in the format 1.2345,123.456,-12.345,12.345 press-ENTER")
extentproject <- extent(as.numeric(scan(what = "double", nmax = 4, dec = ".", sep = ",")))
if (projectprojections == 2)
  {
  writeLines("Please enter projection extent in the format 1.2345,123.456,-12.345,12.345 press-ENTER")
  extentprojectprojection <- extent(as.numeric(scan(what = "double", nmax = 4, dec = ".", sep = ",")))
  }

# Setting default projection for existing layers to longlat WGS84 and setting equal area projection. The same original projection is assumed (longlat) for both layers and projection layers.
bioclimcrs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
projectcrs <- CRS("+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

writeLines("Please enter complete file path of the locations.csv file in the format /home/USER/Desktop/Locations/ press-ENTER")
pathlocations <- scan(what = "character", nmax = 1)
setwd(pathlocations)
projectlocations <- read.csv(file = "locations.csv")
templocations <- matrix(ncol = 2, nrow = 357)
templocations[,1] <- projectlocations$decimallongitude
templocations[,2] <- projectlocations$decimallatitude
templocations <- SpatialPoints(coords = templocations, proj4string = bioclimcrs)
templocationsrp <- spTransform(x = templocations, CRSobj = projectcrs)
projectlocations[,3] <- templocationsrp@coords[,1]
projectlocations[,2] <- templocationsrp@coords[,2]

# Setting layer names and cropping layers' extent to project extent and reprojecting to eaqual area global grid (epsg3975)
setwd(dir = pathlayers)
biofnames <- list.files(path = ".", pattern = "^bio.*\\.asc") # before running bio1 to bio9 should be renamed to bio01 to bio09 to avoid wrong order
elevationfnames <- "topo_Elevation.asc"
bio01 <- projectRaster(crop(x = raster(biofnames[1], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio02 <- projectRaster(crop(x = raster(biofnames[2], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio03 <- projectRaster(crop(x = raster(biofnames[3], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio04 <- projectRaster(crop(x = raster(biofnames[4], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio05 <- projectRaster(crop(x = raster(biofnames[5], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio06 <- projectRaster(crop(x = raster(biofnames[6], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio07 <- projectRaster(crop(x = raster(biofnames[7], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio08 <- projectRaster(crop(x = raster(biofnames[8], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio09 <- projectRaster(crop(x = raster(biofnames[9], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio10 <- projectRaster(crop(x = raster(biofnames[10], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio11 <- projectRaster(crop(x = raster(biofnames[11], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio12 <- projectRaster(crop(x = raster(biofnames[12], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio13 <- projectRaster(crop(x = raster(biofnames[13], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio14 <- projectRaster(crop(x = raster(biofnames[14], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio15 <- projectRaster(crop(x = raster(biofnames[15], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio16 <- projectRaster(crop(x = raster(biofnames[16], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio17 <- projectRaster(crop(x = raster(biofnames[17], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio18 <- projectRaster(crop(x = raster(biofnames[18], crs = bioclimcrs), y = extentproject), crs = projectcrs)
bio19 <- projectRaster(crop(x = raster(biofnames[19], crs = bioclimcrs), y = extentproject), crs = projectcrs)
topo_elevation <- projectRaster(crop(x = raster(elevationfnames, crs = bioclimcrs), y = extentproject), crs = projectcrs)

# As above for projections
if (projectprojections == 1)
  {
  setwd(dir = pathprojections)
  ptbiofnames <- list.files(path = ".", pattern = "^bio.*\\.asc")
  ptelevationfnames <- "topo_Elevation.asc"
  pbio01 <- projectRaster(crop(x = raster(ptbiofnames[1], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio02 <- projectRaster(crop(x = raster(pbiofnames[2], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio03 <- projectRaster(crop(x = raster(pbiofnames[3], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio04 <- projectRaster(crop(x = raster(pbiofnames[4], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio05 <- projectRaster(crop(x = raster(pbiofnames[5], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio06 <- projectRaster(crop(x = raster(pbiofnames[6], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio07 <- projectRaster(crop(x = raster(pbiofnames[7], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio08 <- projectRaster(crop(x = raster(pbiofnames[8], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio09 <- projectRaster(crop(x = raster(pbiofnames[9], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio10 <- projectRaster(crop(x = raster(pbiofnames[10], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio11 <- projectRaster(crop(x = raster(pbiofnames[11], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio12 <- projectRaster(crop(x = raster(pbiofnames[12], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio13 <- projectRaster(crop(x = raster(pbiofnames[13], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio14 <- projectRaster(crop(x = raster(pbiofnames[14], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio15 <- projectRaster(crop(x = raster(pbiofnames[15], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio16 <- projectRaster(crop(x = raster(pbiofnames[16], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio17 <- projectRaster(crop(x = raster(pbiofnames[17], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio18 <- projectRaster(crop(x = raster(pbiofnames[18], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  pbio19 <- projectRaster(crop(x = raster(pbiofnames[19], crs = bioclimcrs), y = extentproject), crs = projectcrs)
  ptopo_elevation <- projectRaster(crop(x = raster(pelevationfnames, crs = bioclimcrs), y = extentproject), crs = projectcrs)
  }
if (projectprojections == 2)
  {
  setwd(dir = pathprojections)
  pbiofnames <- list.files(path = ".", pattern = "^bio.*\\.asc")
  pelevationfnames <- "topo_Elevation.asc"
  pbio01 <- projectRaster(crop(x = raster(pbiofnames[1], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio02 <- projectRaster(crop(x = raster(pbiofnames[2], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio03 <- projectRaster(crop(x = raster(pbiofnames[3], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio04 <- projectRaster(crop(x = raster(pbiofnames[4], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio05 <- projectRaster(crop(x = raster(pbiofnames[5], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio06 <- projectRaster(crop(x = raster(pbiofnames[6], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio07 <- projectRaster(crop(x = raster(pbiofnames[7], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio08 <- projectRaster(crop(x = raster(pbiofnames[8], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio09 <- projectRaster(crop(x = raster(pbiofnames[9], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio10 <- projectRaster(crop(x = raster(pbiofnames[10], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio11 <- projectRaster(crop(x = raster(pbiofnames[11], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio12 <- projectRaster(crop(x = raster(pbiofnames[12], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio13 <- projectRaster(crop(x = raster(pbiofnames[13], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio14 <- projectRaster(crop(x = raster(pbiofnames[14], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio15 <- projectRaster(crop(x = raster(pbiofnames[15], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio16 <- projectRaster(crop(x = raster(pbiofnames[16], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio17 <- projectRaster(crop(x = raster(pbiofnames[17], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio18 <- projectRaster(crop(x = raster(pbiofnames[18], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  pbio19 <- projectRaster(crop(x = raster(pbiofnames[19], crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  ptopo_elevation <- projectRaster(crop(x = raster(pelevationfnames, crs = bioclimcrs), y = extentprojectprojection), crs = projectcrs)
  }


# Setting NA values, resample elevation and compute its slope and aspect. Also compute H (Shannon Entropy), E (Pielou Evenness), Q (Q Statistic) and NM (Fager's Number of Moves)
NAvalue(bio01) <- -9999
NAvalue(bio02) <- -9999
NAvalue(bio03) <- -9999
NAvalue(bio04) <- -9999
NAvalue(bio05) <- -9999
NAvalue(bio06) <- -9999
NAvalue(bio07) <- -9999
NAvalue(bio08) <- -9999
NAvalue(bio09) <- -9999
NAvalue(bio10) <- -9999
NAvalue(bio11) <- -9999
NAvalue(bio12) <- -9999
NAvalue(bio13) <- -9999
NAvalue(bio14) <- -9999
NAvalue(bio15) <- -9999
NAvalue(bio16) <- -9999
NAvalue(bio17) <- -9999
NAvalue(bio18) <- -9999
NAvalue(bio19) <- -9999

# For projections
NAvalue(pbio01) <- -9999
NAvalue(pbio02) <- -9999
NAvalue(pbio03) <- -9999
NAvalue(pbio04) <- -9999
NAvalue(pbio05) <- -9999
NAvalue(pbio06) <- -9999
NAvalue(pbio07) <- -9999
NAvalue(pbio08) <- -9999
NAvalue(pbio09) <- -9999
NAvalue(pbio10) <- -9999
NAvalue(pbio11) <- -9999
NAvalue(pbio12) <- -9999
NAvalue(pbio13) <- -9999
NAvalue(pbio14) <- -9999
NAvalue(pbio15) <- -9999
NAvalue(pbio16) <- -9999
NAvalue(pbio17) <- -9999
NAvalue(pbio18) <- -9999
NAvalue(pbio19) <- -9999
NAvalue(ptopo_elevation) <- -9999

# Terrain analysis
topo_elevation <- resample(x = topo_elevation, y = bio01, method = "bilinear")
topo_slope <- terrain(x = topo_elevation, opt = "slope", unit = "degrees", neighbors=8)
topo_aspect <- terrain(x = topo_elevation, opt = "aspect", unit = "degrees", neighbors=8)
vaspectrcl <- c(67,112,1,112,157,2,157,202,3,202,247,4,247,292,5,292,337,6,337,360,7,0,23,7,0,0,NA,23,67,8) # after reclassification aspect is categorical: 1 N, 2 NW, 3 W...
aspectrcl <- matrix(data = vaspectrcl, ncol = 3, byrow = TRUE)
topo_aspect <- reclassify(topo_aspect,rcl = aspectrcl)
NAvalue(topo_slope) <- -9999
NAvalue(topo_aspect) <- -9999

# For projections
if ((projectprojections == 1) | (projectprojections == 2))
  {
  ptopo_elevation <- resample(x = ptopo_elevation, y = pbio01, method = "bilinear")
  ptopo_slope <- terrain(x = ptopo_elevation, opt = "slope", unit = "degrees", neighbors=8)
  ptopo_aspect <- terrain(x = ptopo_elevation, opt = "aspect", unit = "degrees", neighbors=8)
  ptopo_aspect <- reclassify(ptopo_aspect,rcl = aspectrcl)
  NAvalue(ptopo_slope) <- -9999
  NAvalue(ptopo_aspect) <- -9999
  }

# Save all reprojected rasters as .asc files and the reprojected locations as locations.csv in the project layers and locations folder, respectively
setwd(dir = pathproject)
write.asc(x = asc.from.raster(bio01), "./layers/bio01.asc")
write.asc(x = asc.from.raster(bio02), "./layers/bio02.asc")
write.asc(x = asc.from.raster(bio03), "./layers/bio03.asc")
write.asc(x = asc.from.raster(bio04), "./layers/bio04.asc")
write.asc(x = asc.from.raster(bio05), "./layers/bio05.asc")
write.asc(x = asc.from.raster(bio06), "./layers/bio06.asc")
write.asc(x = asc.from.raster(bio07), "./layers/bio07.asc")
write.asc(x = asc.from.raster(bio08), "./layers/bio08.asc")
write.asc(x = asc.from.raster(bio09), "./layers/bio09.asc")
write.asc(x = asc.from.raster(bio10), "./layers/bio10.asc")
write.asc(x = asc.from.raster(bio11), "./layers/bio11.asc")
write.asc(x = asc.from.raster(bio12), "./layers/bio12.asc")
write.asc(x = asc.from.raster(bio13), "./layers/bio13.asc")
write.asc(x = asc.from.raster(bio14), "./layers/bio14.asc")
write.asc(x = asc.from.raster(bio15), "./layers/bio15.asc")
write.asc(x = asc.from.raster(bio16), "./layers/bio16.asc")
write.asc(x = asc.from.raster(bio17), "./layers/bio17.asc")
write.asc(x = asc.from.raster(bio18), "./layers/bio18.asc")
write.asc(x = asc.from.raster(bio19), "./layers/bio19.asc")
write.asc(x = asc.from.raster(topo_elevation), "./layers/topo_elevation.asc")
write.asc(x = asc.from.raster(topo_slope), "./layers/topo_slope.asc")
write.asc(x = asc.from.raster(topo_aspect), "./layers/topo_aspect.asc")
write.csv(x = projectlocations, file = "./locations/locations.csv")

# For projections
setwd(dir = pathproject)
write.asc(x = asc.from.raster(pbio01), "./projection/bio01.asc")
write.asc(x = asc.from.raster(pbio02), "./projection/bio02.asc")
write.asc(x = asc.from.raster(pbio03), "./projection/bio03.asc")
write.asc(x = asc.from.raster(pbio04), "./projection/bio04.asc")
write.asc(x = asc.from.raster(pbio05), "./projection/bio05.asc")
write.asc(x = asc.from.raster(pbio06), "./projection/bio06.asc")
write.asc(x = asc.from.raster(pbio07), "./projection/bio07.asc")
write.asc(x = asc.from.raster(pbio08), "./projection/bio08.asc")
write.asc(x = asc.from.raster(pbio09), "./projection/bio09.asc")
write.asc(x = asc.from.raster(pbio10), "./projection/bio10.asc")
write.asc(x = asc.from.raster(pbio11), "./projection/bio11.asc")
write.asc(x = asc.from.raster(pbio12), "./projection/bio12.asc")
write.asc(x = asc.from.raster(pbio13), "./projection/bio13.asc")
write.asc(x = asc.from.raster(pbio14), "./projection/bio14.asc")
write.asc(x = asc.from.raster(pbio15), "./projection/bio15.asc")
write.asc(x = asc.from.raster(pbio16), "./projection/bio16.asc")
write.asc(x = asc.from.raster(pbio17), "./projection/bio17.asc")
write.asc(x = asc.from.raster(pbio18), "./projection/bio18.asc")
write.asc(x = asc.from.raster(pbio19), "./projection/bio19.asc")
write.asc(x = asc.from.raster(ptopo_elevation), "./projection/topo_elevation.asc")
write.asc(x = asc.from.raster(ptopo_slope), "./projection/topo_slope.asc")
write.asc(x = asc.from.raster(ptopo_aspect), "./projection/topo_aspect.asc")

# Determine the size of the rasters to adjust the downstream processes
projectdim <- max(c(ncell(bio01), ncell(pbio1)))
