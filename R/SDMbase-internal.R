#' @keywords internal 

# =======================================================================================
# Requirements
# =======================================================================================
require(sp)
require(raster)
require(rgdal)
require(rgbif)
require(KernSmooth)
require(doParallel) #use multi-core capabilities to speed up big computations

# =======================================================================================
# Functions
# =======================================================================================

# read_asc
.read_asc <- function(file, header = 6) {
  
  # read the header first
  metadata <- scan(file, what = "character", sep = "\n", 
                   nlines = header, quiet = TRUE)
  
  # extract info from metadata
  ncol <- as.integer(rev(unlist(strsplit(
    metadata[grep("ncols", metadata,ignore.case=TRUE)], " ")))[1])
  nrow <- as.integer(rev(unlist(strsplit(
    metadata[grep("nrows", metadata,ignore.case=TRUE)], " ")))[1])
  nodata <- as.integer(rev(unlist(strsplit(
    metadata[grep("NODATA_value", metadata,ignore.case=TRUE)], " ")))[1])
  
  # then read the layer
  layer <- scan(file, what = double(), skip = header, quiet = TRUE)
  
  # conform all nodata value to -9999
  if (nodata != -9999) {
    layer[layer == nodata] <- -9999
    metadata[grep("NODATA_value", metadata)] <- "NODATA_value -9999"
  }
  
  # change the data into matrix
  layer <- matrix(data = layer, nrow = nrow, ncol = ncol, byrow = TRUE)
  
  return(list(metadata = metadata, layer = layer))
}

# =======================================================================================

# write_asc
.write_asc <- function(file, metadata, layer, prefix) {
  file.name <- paste0(dirname(file), paste0("/", prefix, "_"), 
                      basename(file))
  cat(metadata, paste0(" ", apply(layer, 1, paste, collapse = " ")), 
      sep = "\n", file = file.name)
  cat("\nCompleted.", prefix, "result is saved at:\n", file.name, "\n")
}

# =======================================================================================

# normalize
.normalize <- function(layer, invert) {
  nonnaindex <- which(layer != -9999)
  layer[nonnaindex] <- (layer[nonnaindex] - min(layer[nonnaindex])) / 
    diff(range(layer[nonnaindex]))
  # invert the values?
  if (invert == TRUE) {
    layer[nonnaindex] <- abs(layer[nonnaindex] - 1)
  }
  return(layer)
}

# =======================================================================================

# alternative normalisation using existing packages x is raster, y is desired levels ie. the multiplicator after 0..1
.norm.raster <- function(x, y) round(((x - cellStats(x, min)) / diff(cellStats(x, range))) * y)

# moving window algorithm
.moving_window <- function(layer, size, FUN_list, keep) {
  # expand layer edges by window size - 1 / 2... 
  # ...so moving windows can move along edges
  expand <- (size - 1) / 2
  layer_expand <- matrix(-9999, dim(layer)[1] + expand * 2, 
                         dim(layer)[2] + expand * 2)
  x_inner_idx <- (expand + 1):(dim(layer)[1] + (expand))
  y_inner_idx <- (expand + 1):(dim(layer)[2] + (expand))
  layer_expand[x_inner_idx, y_inner_idx] <- layer
  
  # initialize
  result <- matrix(-9999, dim(layer)[1], dim(layer)[2])
  result <- rep(list(result), length(FUN_list))
  names(result) <- names(FUN_list)
  
  # loop thru
  for (i in x_inner_idx) {
    # progress output in console
    cat("\rProcessing raster ln:", i - expand, "/", nrow(layer))
    flush.console()      
    
    for(j in y_inner_idx) {
      focal <- layer_expand[i, j]
      # proceed only if the focal cell is not NA
      if (focal != -9999) {
        window <- layer_expand[(i - expand):(i + expand), 
                               (j - expand):(j + expand)]
        # proceed if only windows contain more values than specified threshold
        if (length(window != -9999) > keep) {
          i_idx <- i - expand
          j_idx <- j - expand
          no_na_win <- window[window != -9999]
          # ---------------calculation-----------------
          for (k in 1: length(FUN_list)) {
            result[[k]][i_idx, j_idx] <- FUN_list[[k]](no_na_win)
          }               
          # ---------------calculation-----------------
        }
      }
    }
  }
  return(result)
}

# =======================================================================================

# turn window into presence matrix
.presence_matrix <- function(window, nodata) {
  # convert into presence matrix
  x <- summary(as.factor(window))
  # remove no data cell
  nodata.idx <- which(names(x) == as.character(nodata))
  if (length(nodata.idx) > 0)
    x <- x[-nodata.idx]
  return(x)
}

# =======================================================================================

# calculate Shanon diversity
.ShanonH <- function(window, nodata = -9999) {
  x <- .presence_matrix(window, nodata = nodata)
  # total individuals
  N <- sum(x)
  # calculate p (proportion of data point belongs to i-th type over all data points)
  p <- x / N
  H <- sum(-p * log(p))
  return(H)
}

# =======================================================================================

# calculate Pielou's Evenness
#.PielouJ <- function(window, nodata = -9999) {
#  H <- .ShanonH(window = window, nodata = nodata)
#  x <- .presence_matrix(window, nodata = nodata)
#  S <- length(x) # total type number
#  J <- H / log(S)
#  return(J)
#}

# =======================================================================================

# calculate Fager's number of moves
# based on Fager(1972) p. 300
.NMS <- function(window, nodata = -9999) {
  x <- .presence_matrix(window, nodata = nodata)
  N <- sum(x)
  S <- length(x)
  R <- rank(-x, na.last = NA, ties.method = "average")
  NMS <- ((N * (S + 1)) / 2) - sum(R * x) # corrected formula, now values are zero or positive
  return(NMS)
}

# =======================================================================================

# calculate Q statistics
.Qstat <- function(window, nodata = -9999) {
  x <- .presence_matrix(window, nodata = nodata)
  S <- length(x) # total type number
  Q1 <- ceiling(S / 4)
  Q3 <- ceiling(3 * S / 4)
  R1 <- sort(x)[Q1]
  R2 <- sort(x)[Q3]
  sum_n_r <- sum(which(R1 < sort(x) & sort(x) < R2))
  nR1 <- sort(x) == sort(x)[Q1]
  nR2 <- sort(x) == sort(x)[Q3]
  # follow equation from Magurran (2004) p. 103
  if (log(R2/R1) > 0) # report Q only if divisor is not 0, divide 0 doesnt make sense
    Q <- sum(nR1/2, sum_n_r, nR2/2) / (log(R2/R1))
  else
    Q <- NA
  return(Q)
}

# =======================================================================================
# The remainder of this file contains functions used to prepare the data for annd
# to run maxent
# =======================================================================================

# retrieve zipped bioclim data and unzip it
.retrieveBio <- function(rBurl1 = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio1-9_30s_bil.zip", rBurl2 = "http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio10-19_30s_bil.zip", rBmethod = "internal", bioclimraw = "./SDMbaserawraster/") {
  dir.create(path = rBdest)
  download.file(url = rBurl1,   destfile = paste0(bioclimraw,"bio1-9_30s_bil.zip"), method =  rBmethod)
  download.file(url = rBurl2,   destfile = paste0(bioclimraw,"bio10-19_30s_bil.zip"), method =  rBmethod)
  unzip(zipfile = paste0(bioclimraw,"bio1-9_30s_bil.zip"), exdir = bioclimraw)
  unzip(zipfile = paste0(bioclimraw,"bio10-19_30s_bil.zip"), exdir = bioclimraw)
  file.remove(paste0(bioclimraw,"bio1-9_30s_bil.zip"), paste0(bioclimraw,"bio10-19_30s_bil.zip"))
}

# =======================================================================================
# Read in, crop, reproject, mask and write the raster
.prepareRaster <- function(bioclimraw = "./SDMbaserawraster/", bioclimtype = "bil", bioclimCRSdef = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0", outputCRSdef = "+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0",projectextent = projectextent,maskfolder = "./", maskfile = "fname",domask = FALSE, maskCRS = NULL, resolution = 1000, layerfolder = "./example/Layers/") {

  inputCRS <- CRS(bioclimCRSdef)
  outputCRS <- CRS(outputCRSdef)
  if (bioclimtype == "bil") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.bil")
  }
    if (bioclimtype == "grd") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.grd")
  }
    if (bioclimtype == "asc") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.asc")
  }

 # Correct the order of files
  layerfnames <- layerfnames[order(as.numeric(sub("([0-9]*).*", "\\1", layerfnames)))]
  
 # Calculate!
 layerstack <- stack() 
 for (i in 1:length(layerfnames)) {
    layerstack <- addLayer(layerstack, projectRaster(crop(x = raster(layerfnames[i], crs = inputCRS), y = projectextent), crs = outputCRS, res = resolution))
  }
 # Mask
  if (domask == TRUE) {
    coastlinein <- readOGR(dsn = maskfolder, layer = maskfile, p4s = maskCRS)
    coastlineout <- spTransform(x = coastlinein, CRSobj = outputCRS)
    layerstack <- mask(layerstack, coastlineout, updatevalue = NA)
  }
  NAvalue(layerstack) <- -9999

 # Write the asc files
  for (i in (1:length(layerfnames))) {
    writeRaster(x = paste0("layerstack$bio",i), file = paste0(layerfolder, "bio", i), format = "ascii", NAflag = -9999)
  }
}

# =======================================================================================
# Read in, crop, reproject, mask and write the raster with 2 or more cores 
.prepareRasterMC <- function(pRsource = "./SDMbaserawraster/", pRtype = "bil", pRinputCRSdef = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0", pRoutputCRSdef = "+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0", pRmasksource = "./", pRmaskfname = "fname", pRmask = FALSE, pRresolution = 1000, pRdest = "./example/Layers/", ncores = 2) {
    
  inputCRS <- CRS(bioclimCRSdef)
  outputCRS <- CRS(outputCRSdef)
  if (bioclimtype == "bil") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.bil")
  }
  if (bioclimtype == "grd") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.grd")
  }
  if (bioclimtype == "asc") {
    layerfnames <- list.files(path = bioclimraw, pattern = "^bio.*\\.asc")
  }

 # Correct the order of files
  layerfnames <- layerfnames[order(as.numeric(sub("([0-9]*).*", "\\1", layerfnames)))]

 # Setup multicore cluster 
  tempclust <- makeCluster(ncores)
  registerDoParallel(tempclust)
  layerstack <- stack()
 # Running multicore code
  foreach (i = 1:19, .packages = 'raster') %dopar% {
    paste0("bio",i) <- projectRaster(crop(x = raster(layerfnames[i], crs = inputCRS), y = projectextent), crs = outputCRS, res = resolution)
  }
  layerstack <- stack(bio1:bio19)  
    if (domask == TRUE) {
        coastlinein <- readOGR(dsn = maskfolder, layer = maskfile, p4s = inputCRS)
        coastlineout <- spTransform(x = coastlinein, CRSobj = outputCRS)
        layerstack <- mask(layerstack, coastlineout,updatevalue = NA)
    }
    NAvalue(layerstack) <- -9999
    # Write the asc files
    foreach (i = 1:length(layerfnames)) %dopar% {
        writeRaster(x = paste0("layerstack$bio",i), file = paste0(pRdest, "bio", i), format = "ascii", NAflag = -9999)
    }
    stopCluster(tempclust)
}

# =======================================================================================
# Reproject locations and write locations file if required. file is required to be without any header and in the format: species, declon, declat

.rLocations <- function (locationfolder = "./", locationfile = "locations.csv", inputCRSdef = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0", outputCRSdef = "+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0", dowritefile = TRUE) {
  inputCRS <- CRS(inputCRSdef)
  outputCRS <- CRS(outputCRSdef)
  rLlocations <- base::unique(read.csv(file = paste0(locationfolder, locationfile), header = FALSE))
  templocations <- matrix(ncol = 2, nrow = length(rLlocations$V3)) #make sure this is the correct length maybe better data@ncols?
  templocations[,1] <- rLlocations$V2 # make sure it is $V2 not [,2]
  templocations[,2] <- rLlocations$V3
  templocations <- SpatialPoints(coords = templocations, proj4string = inputCRS)
  templocationsrp <- spTransform(x = templocations, CRSobj = outputCRS)
  rLlocations[,2] <- templocationsrp@coords[,1]
  rLlocations[,3] <- templocationsrp@coords[,2]
  if (dowritefile == TRUE) {
    write.csv(x = rLlocations, file = paste0(locationfolder, "projlocations.csv"), sep = ",", quote = FALSE)
  }
}

# =======================================================================================
# Create a bias file for Maxent based on presence, absence and presence and absence of similar genera to factor out sampling bias

.createBiasfile <- function(locationfolder = "./Locations/", biaslocations = "samplinglocations.csv", cBinputCRSdef = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0", cBoutputCRSdef = "+init=epsg:3975 +proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0", cBbwidth = 100000, cBrange = NULL, layerfolder = "./", layerfile = "template.asc", cBsimilar = NULL, cBgeometry = NULL) {
  inputCRS <- CRS(cBinputCRSdef)
  outputCRS <- CRS(cBoutputCRSdef)
  cBrecords <- base::unique(read.csv(paste0(locationfolder, biaslocations))) # file containing all presences, absences and more data of similar species. genus, declon, declat
  if (cBsimilar != NULL) {
  cBkeys <- sapply(cBsimilar, function(x) name_suggest(x, rank = "genus")$key[1], USE.NAMES=FALSE)
  cBaddrecords <- occ_search(taxonKey = cBkeys, hasCoordinate = TRUE, hasGeospatialIssue = FALSE, geometry = cBgeometry, year = paste0("1959,", format(Sys.Date(), "%Y")), fields = "minimal", return = "data")
  }
  # how to combine two data frames? rbind!
  cBrecords <- rbind(matrix(cBrecords), matrix(cBaddrecords))
  templocations <- matrix(ncol = 2, nrow = length(cBrecords$V3)) # again is this length correct?
  templocations[,1] <- cBrecords[,2]
  templocations[,2] <- cBrecords[,3]
  templocations <- SpatialPoints(coords = templocations, proj4string = inputCRS)
  templocations <- spTransform(x = templocations, CRSobj = outputCRS)

  template <- raster(x = paste0(layerfolder, layerfile), crs = outputCRS)
  cBrange <- list( c(xyFromCell(template, cell = c(ncell(template) - (ncol(template)-1)))[1], xyFromCell(template, cell = c(ncell(template) - (ncol(template)-1)))[2]), c(xyFromCell(template, cell = ncol(template))[1], xyFromCell(template, cell = ncol(template))[2]))
  est <- bkde2D(templocations@coords, bandwidth=c(cBwidth,cBwidth), gridsize = c(paste0(template@ncols,"L"), paste0(template@nrows,"L")), range.x = cBrange, truncate = FALSE)
  est.raster <- raster(list(x=est$x1,y=est$x2,z=est$fhat))
  est.raster@crs <- outputCRS
  est.raster@extent <- template@extent
  est.raster <- (est.raster - cellStats(est.raster, min))/(cellStats(est.raster,max)-cellStats(est.raster,min)) # normalization or feature scaling LATER make sure to reuse Jin's code above!
  est.raster <- (est.raster * 4) + 1 # scale with 4 and add 1 to fullfill maxent demands of non zero non negative
  est.raster <- mask(x = est.raster, mask = template) #template is any raster ready for maxent use
  writeRaster(x = est.raster, filename = paste0(locationfolder,"biasraster.asc"), format = "ascii", NAflag = -9999)
}

# =======================================================================================
# Run maxent within R
.Maxent <- function() {
maxentcommand <- paste("java",,,,,,,,,,,)
system(command)
}
# =======================================================================================
# End of file
# =======================================================================================
