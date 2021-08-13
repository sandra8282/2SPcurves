#' INTERNAL: data management not for user
#'
#' @description
#' This function gets NPMLE survival estimators for each group for doubly censored data.
#'
#' @param dat1 Numeric or character vector of subject's unique identifier (i).
#' @param dat0 Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param dtype Vector with the lengths of time spent in event of Type I for individual i in episode j.
#'
#' @return A list including:
#' \itemize{
#'  \item skm: a matrix with event times, estimated survival probs and indicator for ties.
#'  \item mskm: the lowest value for estimated survival.
#' }
#'
#'@import survival
#'@importFrom dplyr distinct
#'@keywords internal
#'@noRd
#'
#'

getDKMtab <- function(dat0, dat1, maxtime, ctime, dtype) {

  if (dtype == "DIcens"){
  #get control results
  res0 <- DGLalg(o_l = dat0$o_l, o_r = dat0$o_r,
                 e_l = dat0$e_l, e_r = dat0$e_r, maxtime=maxtime, ctime=ctime)
  #get control results
  res1 <- DGLalg(o_l = dat1$o_l, o_r = dat1$o_r,
                 e_l = dat1$e_l, e_r = dat1$e_r, maxtime=maxtime, ctime=ctime)

  #set up skm
  dskm0 <- cbind(res0$gap[ , c(1, 4)], type = rep(0, nrow(res0$gap)))
  dskm1 <- cbind(res1$gap[ , c(1, 4)], type = rep(1, nrow(res1$gap)))

  }

  dskm <- rbind(dskm0, dskm1)
  dskm <- dskm[order(dskm[,1], dskm[,3]),]
  ties_check <- unique(table(dskm[,1]))

  if (length(ties_check) > 1) {
    ties_times = dskm[duplicated(dskm[,1]),1]
    ties_ind <- rep(0, nrow(dskm))
    ties_ind[which(dskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(dskm))}

  dskm = cbind(dskm, ties_ind)

  dmskm <- min(dskm[,2])

  dskm <- as.matrix(distinct(data.frame(dskm)))

  return(list(dskm = dskm, dmskm = dmskm,
              npmle0 = res0, npmle1 = res1))



}
