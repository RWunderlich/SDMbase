require(sp)
require(raster)
run <- function(pathmaxent="/home/affu/Desktop/Maxent/", pathproject="/home/affu/Desktop/example/", doproject=0,
pathlayers="/home/affu/Desktop/Layers/", pathplayers="/home/affu/Desktop/Layers/Players", inputfiles="asc",
dobgmanip=FALSE, updownbgmanip=0, levelbgmanip=300, inputCRS=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"),
outputCRS=CRS("+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"),
projectextent=c(60,110,8.5,46), projectpextent=c(), pathlocations="/home/affu/Desktop/Locations/")
  {
  # Creating project folder structure
  setwd(dir = pathproject)
  dir.create(path = "./layers/")
  dir.create(path = "./locations/")
  dir.create(path = "./results/")
  if ((doproject == 1) | (doproject == 2))
    {
      dir.create(path = "./projection/")
    }
  # Reading and reprojecting location data
  setwd(dir = pathlocations)
  projectlocations <- read.csv(file = "locations.csv")
  templocations <- matrix(ncol = 2, nrow = 357)
  templocations[,1] <- projectlocations$decimallongitude
  templocations[,2] <- projectlocations$decimallatitude
  templocations <- SpatialPoints(coords = templocations, proj4string = inputCRS)
  templocationsrp <- spTransform(x = templocations, CRSobj = outputCRS)
  projectlocations[,3] <- templocationsrp@coords[,1]
  projectlocations[,2] <- templocationsrp@coords[,2]
  # Reading in, cropping to extent, reprojecting to equal area
  setwd(dir = pathlayers)
  layerfnames <- c(list.files(path = ".", pattern = "^bio.*\\.asc"), "topo_Elevation.asc")
  layerstack <- stack()
  NAvalue(layerstack) <- -9999
  for (i in 1:(length(layerfnames)-1)
    {
    layerstack[i] <- projectRaster(crop(x = raster(layerfnames[i], crs = inputCRS), y = projectextent), crs = outputCRS)
    }
    layerstack[20] <- resample(x = raster(layerfnames[20]), y = layerstack, method = "bilinear")
  # Insert code to read in projection layers if.. for i in...
  if (doproject == 1)
    {
    setwd(dir = pathplayers)
    playerfnames <- c(list.files(path = ".", pattern = "^bio.*\\.asc"), "topo_Elevation.asc")
    playerstack <- stack()
    NAvalue(playerstack) <- -9999
    for (i in 1:(length(playerfnames)-1)
      {
      playerstack[i] <- projectRaster(crop(x = raster(playerfnames[i], crs = inputCRS), y = projectextent), crs = outputCRS)
      }
    playerstack[20] <- resample(x = raster(playerfnames[20]), y = playerstack, method = "bilinear")    
    }

  if (doproject == 2)
    {
    setwd(dir = pathlayers)
    playerfnames <- c(list.files(path = ".", pattern = "^bio.*\\.asc"), "topo_Elevation.asc")
    playerstack <- stack()
    NAvalue(playerstack) <- -9999
    for (i in 1:(length(playerfnames)-1)
      {
      playerstack[i] <- projectRaster(crop(x = raster(playerfnames[i], crs = inputCRS), y = projectpextent), crs = outputCRS)
      }
    playerstack[20] <- resample(x = raster(playerfnames[20]), y = playerstack, method = "bilinear")
    }

  # Terrain analysis
  layerstack[21] <- terrain(x = layerstack[20], opt = "slope", unit = "degrees", neighbors=8)
  layerstack[22] <- terrain(x = layerstack[20], opt = "aspect", unit = "degrees", neighbors=8)

  if ((doproject == 1) | (doproject == 2))
    {
    playerstack[21] <- terrain(x = playerstack[20], opt = "slope", unit = "degrees", neighbors=8)
    playerstack[22] <- terrain(x = playerstack[20], opt = "aspect", unit = "degrees", neighbors=8)
    }

  # Background manipulation if required
  if dobgmanip==TRUE
    {
    for (i in 1:(length(layerstack))
      {
      if (updownbgmanip == 0)
        {
        layerstack[i] <- layerstack[i] * (layerstack[20] < levelbgmanip)
        }
      if (updownbgmanip == 1)
        {
        layerstack[i] <- layerstack[i] * (layerstack[20] > levelbgmanip)
        }
      }
    if ((doproject == 1) | (doproject == 2))
      {
      for (i in 1:(length(playerstack))
        {
        if (updownbgmanip == 0)
          {
          playerstack[i] <- playerstack[i] * (playerstack[20] < levelbgmanip)
          }
        if (updownbgmanip == 1)
          {
          playerstack[i] <- playerstack[i] * (playerstack[20] > levelbgmanip)
          }
        }
      }
    }

  # Writing the .csv and .asc output files
  setwd(dir = pathproject)
  write.csv(x = projectlocations, file = "./locations/locations.csv", sep = ",")
  for (i in (1:9))
    {
    write.asciigrid(x = layerstack[i], file = "./layers/bio0,i,\.asc", na.value = -9999)
    }
  for (i in (10:19))
    {
    write.asciigrid(x = layerstack[i], file = "./layers/bio,i,\.asc", na.value = -9999)
    }
  write.asciigrid(x = layerstack[20], file = "./layers/topo_Elevation.asc", na.value = -9999)
  write.asciigrid(x = layerstack[21], file = "./layers/topo_Slope.asc", na.value = -9999)
  write.asciigrid(x = layerstack[22], file = "./layers/topo_Aspect.asc", na.value = -9999)

  # Writing projection layers if they exist
  if ((doproject == 1)|(doproject == 2))
    {
    for (i in (1:9))
      {
      write.asciigrid(x = playerstack[i], file = "./projection/bio0,i,\.asc", na.value = -9999)
      }
    for (i in (10:19))
      {
      write.asciigrid(x = layerstack[i], file = "./projection/bio,i,\.asc", na.value = -9999)
      }
    write.asciigrid(x = layerstack[20], file = "./projection/topo_Elevation.asc", na.value = -9999)
    write.asciigrid(x = layerstack[21], file = "./projection/topo_Slope.asc", na.value = -9999)
    write.asciigrid(x = layerstack[22], file = "./projection/topo_Aspect.asc", na.value = -9999)
    }
  }
