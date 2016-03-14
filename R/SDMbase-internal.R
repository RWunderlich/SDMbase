#' @keywords internal 

# read_asc
.read_asc <- function(file, header = 6) {
  
  # read the header first
  metadata <- scan(file, what = "character", sep = "\n", 
                   nlines = header, quiet = TRUE)
  
  # extract info from metadata
  ncol <- as.integer(rev(unlist(strsplit(
    metadata[grep("ncols", metadata)], " ")))[1])
  nrow <- as.integer(rev(unlist(strsplit(
    metadata[grep("nrows", metadata)], " ")))[1])
  nodata <- as.integer(rev(unlist(strsplit(
    metadata[grep("NODATA_value", metadata)], " ")))[1])
  
  # then read the layer
  layer <- scan(file, what = integer(), skip = header, quiet = TRUE)
  
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
.PielouJ <- function(window, nodata = -9999) {
  H <- .ShanonH(window = window, nodata = nodata)
  x <- .presence_matrix(window, nodata = nodata)
  S <- length(x) # total type number
  J <- H / log(S)
  return(J)
}

# =======================================================================================

# calculate Fager's number of moves
# based on Fager(1972) p. 300
.NMS <- function(window, nodata = -9999) {
  x <- .presence_matrix(window, nodata = nodata)
  N <- sum(x)
  S <- length(x)
  R <- order(x, decreasing = TRUE) 
  NMS <- N * (S + 1) / 2 - sum(R * x)  
}

# =======================================================================================

## calculate Q statistics
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
