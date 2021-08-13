#' INTERNAL: data management not for user
#'
#' @description
#' This function gets AKM survival estimators for each group for right censored data.
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param weights Vector with the weights for each subject (i).
#'
#' @return A list including:
#' \itemize{
#'  \item skm: a matrix with event times, estimated survival probs and indicator for ties.
#'  \item mskm: the lowest value for estimated survival.
#' }
#'
#'@import survival
#'@importFrom dplyr distinct
#'@importFrom RISCA ipw.survival
#'@keywords internal
#'@noRd
#'
#'

getAKMtab <- function(time, event, group, weights) {

  akm <- ipw.survival(times = time, failures = event,
                                  variable = group, weights = weights)

  skm <- akm[[1]][,-c(2:3)]
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

  return(list(skm = skm, mskm = mskm,
              adjkm = akm[[1]]))

}
