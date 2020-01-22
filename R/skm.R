#' INTERNAL: data management not for user
#'
#' @description
#' This function gets KM survival estimators for each group.
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#'
#' @return A list including:
#' \itemize{
#'  \item skm: a matrix with event times, estimated survival probs and indicator for ties.
#'  \item mskm: the lowest value for estimated survival.
#' }
#'
#'@import survival
#'@keywords internal
#'
getKMtab <- function(time, event, group) {

  #Data mgmt
  km_placebo <- survfit(Surv(time, event) ~ 1,
                        subset=(group==0), type='kaplan-meier')
  km_drug <- survfit(Surv(time, event) ~ 1,
                     subset=(group==1), type='kaplan-meier')

  skm_p <- cbind(time = summary(km_placebo)$time,
                 surv = summary(km_placebo)$surv)
  skm_p <- cbind(skm_p, type = rep(0, nrow(skm_p)))
  skm_d <- cbind(time = summary(km_drug)$time,
                 surv = summary(km_drug)$surv)
  skm_d <- cbind(skm_d, type = rep(1, nrow(skm_d)))

  skm <- rbind(skm_p, skm_d)
  skm <- skm[order(skm[,1]),]

  ties_check <- unique(table(skm[,1]))

  if (length(ties_check) > 1) {
    ties_times = as.integer(names(which(table(skm[,1])>1)))
    ties_ind <- rep(NA, nrow(skm))
    ties_ind[which(skm[,1] %in% ties_times)]=1
    ties_ind[-which(skm[,1] %in% ties_times)]=0
  } else {ties_ind <- rep(0, nrow(skm))}

  skm = cbind(skm, ties_ind)

  mskm <- min(skm[,2])

  return(list(skm = skm, mskm = mskm,
              km_placebo = km_placebo, km_drug = km_drug))

}
