getKM <- function(x,d){
  tempsumm <- summary(survfit(Surv(x, d) ~ 1, type='kaplan-meier'))
  tempsumm <- cbind(time = tempsumm$time, surv = tempsumm$surv)
  return(tempsumm)
}

ljoinf <- function(alist, adf, B){
  for (b in 1:B){
    ind <- which(adf$time %in% alist[[b]][,1])
    adf[ind,b+1] <- alist[[b]][,2]
  }
  return(adf)
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)

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
getKMtabRecurrent <- function(dat) {

  B <- 1000
  #Take out censored observations for subjects with more than 1 gap time
  outind <- which(dat$delta==1 & dat$mi == dat$episode_j)
  dat_fs <- dat[-outind, ]

  #set-up matrix for samples
  unique_i <- dat_fs[dat_fs$episode_j==1, c(1,2,4,6)]
  samples <- cbind(unique_i,
                   matrix(rep(rep(NA, nrow(unique_i)), B), nrow=nrow(unique_i)))

  #
  for (i in (samples$id)){
    temp <- subset(dat_fs, dat_fs$id == i)
    ind <- which(samples$id==i)
    if (nrow(temp)==1){
      samples[ind, 5:ncol(samples)] <- temp$time
    } else{samples[ind, 5:ncol(samples)] <- sample(temp$time, B, replace = TRUE)}
  }

  sample0 <- samples[samples$group==0,]; sample1 <- samples[samples$group==1,];

  #control
  Sb0 <- apply(sample0[,5:ncol(sample0)], 2, function(x) getKM(x, sample0$delta))
  Sb0df <- data.frame(time = seq(0,
                                 max(colMax(sample0[,5:ncol(sample0)])),
                                 min(sapply(Sb0, function(x) min(diff(x[,1]))))))
  Sb0df <- cbind(Sb0df, matrix(rep(rep(NA, nrow(Sb0df)), B), ncol=B))
  Sb0df <- ljoinf(Sb0, Sb0df, B); Sb0df[1, 2:ncol(Sb0df)] <- 1;
  Sb0_t <- data.frame(time = Sb0df$time,
                      surv = rowMeans(Sb0df[,-1], na.rm = TRUE),
                      group = 0)

  #treatment
  Sb1 <- apply(sample1[,5:ncol(sample1)], 2, function(x) getKM(x, sample1$delta))
  Sb1df <- data.frame(time = seq(0,
                                 max(colMax(sample1[,5:ncol(sample1)])),
                                 min(sapply(Sb1, function(x) min(diff(x[,1]))))))
  Sb1df <- cbind(Sb1df, matrix(rep(rep(NA, nrow(Sb1df)), B), ncol=B))
  Sb1df <- ljoinf(Sb1, Sb1df); Sb1df[1, 2:ncol(Sb1df)] <- 1;
  Sb1_t <- data.frame(time = Sb1df$time,
                      surv = rowMeans(Sb1df[,-1], na.rm = TRUE),
                      group = 1)

  skm <- na.omit(rbind(Sb0_t, Sb1_t))
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
              MWCRsurv_placebo = na.omit(Sb0_t),
              MWCRsurv_treatmen = na.omit(Sb1_t)))



}
