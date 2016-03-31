layer_calc <- function(file, header = 6, size = 3, keep,
                       calc = c("cat", "cont"),
                       normalize = TRUE, invert = FALSE) {   
  # argument check
  calc <- match.arg(calc)
  if (missing(keep)) # default keep to size if not specified
    keep <- size
  
  # read the file
  dat <- .read_asc(file, header = header)
  layer <- dat$layer
  
  # ================ cont start ===============================
  if (calc == "cont") {

    # idm is calculated using glcm package directly... 
    # ... using its built-in window algorithm 
    require(glcm)
    idm_result <- glcm(layer, n_grey = length(levels(as.factor(c(layer)))) - 1, 
                       window = c(size, size), statistics = "homogeneity", 
                       na_val = -9999, na_opt = "center")
    
    # moving window calc
    FUN_list <- list(sd = sd, 
                     vmr = function(x) {var(x) / prod(x) ^ (1/length(x))},
                     mad = mad)
    result <- .moving_window(layer = layer, size = size, 
                             FUN_list = FUN_list, keep = keep)
    
    # normalize the result
    if (normalize == TRUE) 
      lapply(result, .normalize, invert = invert)
    
    # export
    for (i in 1:length(FUN_list)) {
      .write_asc(file = file, metadata = dat$metadata, layer = result[[i]], 
                 prefix = names(FUN_list)[i])
    }
    
    # export idm
    .write_asc(file = file, metadata = dat$metadata, layer = idm_result, 
               prefix = "idm_")
  }
  # ================ cont ends ===============================

  # ================ cat starts ==============================
  if (calc == "cat") {  
    
    FUN_list <- list(H = .ShanonH, J = .PielouJ, 
                     HxJ = function(x) .ShanonH(x) * .PielouJ(x), 
                     NMS = .NMS, Q = .Qstat)
    result <- .moving_window(layer = layer, size = size, 
                             FUN_list = FUN_list, keep = keep)

    if (normalize == TRUE) 
      result <- lapply(result, .normalize, invert = invert)
    
    # export
    for (i in 1:length(FUN_list)) {
      .write_asc(file = file, metadata = dat$metadata, layer = result[[i]], 
                 prefix = names(FUN_list)[i])
    }
  } 
  # ================ cat ends ================================
}
