#' INTERNAL: data management not for user
#'
#' @description
#' This function gets a matrix need to draw 2SSP for interval data.
#'
#' @param skm passed from TwoSSPicens
#' @return skm: a matrix with event times, estimated survival probs and indicator for ties.
#'
#'@keywords internal
#'@noRd

getIntSKM <- function(skm){
  skm <- skm[order(skm[,1], skm[,3]),]

  ties_check <- unique(table(skm[,1]))

  if (length(ties_check) > 1) {
    ties_times = skm[duplicated(skm[,1]),1]
    ties_ind <- rep(0, nrow(skm))
    ties_ind[which(skm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(skm))}

  skm = cbind(skm, ties_ind)

  mskm <- min(skm[,2])

  skm <- as.matrix(distinct(data.frame(skm)))

  return(skm)
}
