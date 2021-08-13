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

getMKMtab <- function(entry, exit, event, group, c, a) {

  #Data mgmt
  dat <- data.frame(entry, exit, group, event)
  dat0 <- dat[dat$group==0,];  dat1 <- dat[dat$group==1,];
  n0 = nrow(dat0); n1=nrow(dat1);
  qcomp0 <- c*(n0^a); qcomp1 <- c*(n1^a);

  #fit and calculate modified prod est for not exposed
  km.fit0<-survfit(Surv(entry, exit, event, type='counting') ~ 1, data=dat0, type='kaplan-meier')
  km0 <- summary(km.fit0)
  mkm0 <- data.frame(time=km0$time, nevent = km0$n.event, nrisk = km0$n.risk)
  Ind <- ifelse(mkm0$nrisk >= qcomp0, 1, 0); probs <- 1-(mkm0$nevent/mkm0$nrisk)*Ind;
  mkm0$surv <- cumprod(probs)

  #fit and calculate modified prod est for exposed
  km.fit1<-survfit(Surv(entry, exit, event, type='counting') ~ 1, data=dat1, type='kaplan-meier')
  km1 <- summary(km.fit1);
  mkm1 <- data.frame(time=km1$time, nevent = km1$n.event, nrisk = km1$n.risk)
  Ind <- ifelse(mkm1$nrisk >= qcomp1, 1, 0); probs <- 1-(mkm1$nevent/mkm1$nrisk)*Ind;
  mkm1$surv <- cumprod(probs)

  mskm0 <- cbind(time = mkm0$time, surv = mkm0$surv)
  mskm0 <- cbind(mskm0, type = rep(0, nrow(mskm0)))
  mskm1 <- cbind(time = mkm1$time, surv = mkm1$surv)
  mskm1 <- cbind(mskm1, type = rep(1, nrow(mskm1)))

  mskm <- rbind(mskm0, mskm1)
  mskm <- mskm[order(mskm[,1], mskm[,3]),]

  ties_check <- unique(table(mskm[,1]))

  if (length(ties_check) > 1) {
    ties_times = mskm[duplicated(mskm[,1]),1]
    ties_ind <- rep(0, nrow(mskm))
    ties_ind[which(mskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(mskm))}

  mskm = cbind(mskm, ties_ind)
  mmskm <- min(mskm[,2])

  mskm <- as.matrix(distinct(data.frame(mskm)))

  return(list(mskm = mskm, mmskm = mmskm,
              km0 = km0, km1 = km1, mkm0 = mkm0, mkm1 = mkm1))

}
