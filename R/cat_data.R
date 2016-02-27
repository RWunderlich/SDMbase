#' Calculate diversity indices for layer with categorical data
#' 
#' @description This is an internal function used by \code{\link{layer_calc}},
#'   intended to be used in a moving window.
#' @param window matrix of any size containing characters or integers 
#'   representing different classes in the layer.
#' @param nodata character or integer giving value used for no data.
#' @return 
#' A list of the following indices are returned:
#'  \item{H}{Shanon's diversity}
#'  \item{J}{Pielou's evenness}
#'  \item{Q}{Q statistics}
#'  \item{NMS}{Fager's number of moves}
#' @seealso
#'  \code{\link{layer_calc}}
#' @references 
#' \url{http://www.wcsmalaysia.org/analysis/diversityIndexMenagerie.htm#Qstatistic}
#' 
#' Magurran, A. E. (2004). \emph{Measuring biological diversity}. UK: Blackwell
#' Publishing.
#' 
#' Fager, E. W. (1972). Diversity: a sampling study. \emph{American Naturalist},
#' 293-310.
#' @export

cat_data <- function(window, nodata = -9999) {
  
  # convert into presence matrix
  x <- summary(as.factor(window))
  
  # remove no data cell
  nodata.idx <- which(names(x) == as.character(nodata))
  if (length(nodata.idx) > 0)
    x <- x[-nodata.idx]
  
  # total individuals
  N <- sum(x)
  
  # calculate p (proportion of data point belongs to i-th type over all data points)
  p <- x / N
  
  ## calculate Shanon diversity
  H <- sum(-p * log(p))
  
  ## calculate Pielou's evenness
  S <- length(x) # total type number
  J <- H / log(S)
  
  ## calculate Fager's number of moves
  # based on Fager(1972) p. 300
  R <- order(x, decreasing = TRUE) 
  NMS <- N * (S + 1) / 2 - sum(R * x)  
  
  ## calculate Q statistics
  Q1 <- ceiling(S/4)
  Q3 <- ceiling(3*S/4)
  R1 <- sort(x)[Q1]
  R2 <- sort(x)[Q3]
  sum_n_r <- sum(which(R1 < sort(x) & sort(x) < R3))
  nR1 <- sort(x) == sort(x)[Q1]
  nR2 <- sort(x) == sort(x)[Q3]
  # follow equation from Magurran (2004) p. 103
  if (log(R2/R1) > 0) # report Q only if divisor is not 0, divide 0 doesnt make sense
    Q <- sum(nR1/2, sum_n_r, nR2/2) / (log(R2/R1))
  else
    Q <- NA
  return(list(H = H, J = J, Q = Q, NMS = NMS))
}
