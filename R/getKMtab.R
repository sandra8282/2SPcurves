#' INTERNAL: data management not for user
#'
#' @description
#' This function gets KM survival estimators for each group for right censored data.
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
#'@importFrom dplyr distinct
#'@keywords internal
#'@noRd
#'
#'

getKMtab <- function(time, event, group) {

  #Data mgmt
  check0 = check1 = 1
  if (length(unique(event[group==0]))==1) {
    if (unique(event[group==0])==0) {check0=0}}

  if (length(unique(event[group==1]))==1) {
    if (unique(event[group==1])==0) {check1=0}}

  if (check0==1) {
      km_placebo <- survfit(Surv(time, event) ~ 1,
                        subset=(group==0), type='kaplan-meier')
      skm_p <- cbind(time = summary(km_placebo)$time, surv = summary(km_placebo)$surv)
      skm_p <- cbind(skm_p, type = rep(0, nrow(skm_p)))

  } else {
     km_placebo <- NA
     skm_p <- cbind(time = max(time[group==0]), surv = 1, type=0)
  }

  if (check1==1) {
    km_drug <- survfit(Surv(time, event) ~ 1,
                       subset=(group==1), type='kaplan-meier')
    skm_d <- cbind(time = summary(km_drug)$time, surv = summary(km_drug)$surv)
    skm_d <- cbind(skm_d, type = rep(1, nrow(skm_d)))

  } else {
    km_drug <- NA
    skm_d <- cbind(time = max(time[group==1]), surv = 1, type=1)
  }

  skm <- rbind(skm_p, skm_d)
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
              km_placebo = km_placebo, km_drug = km_drug))

}
