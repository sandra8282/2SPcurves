base_surv_weibull <- function(lambda, alpha, times){
  s0 <- exp(-lambda*(times)^alpha)
  return(s0)
}

getSKM4fit <- function(time, fitsurv, group) {
  fitskm <- cbind(time, fitsurv, group)
  fitskm <- fitskm[!duplicated(fitskm),]
  ties_check <- unique(table(fitskm[,1]))
  if (length(ties_check) > 1) {
    ties_times = fitskm[duplicated(fitskm[,1]),1]
    ties_ind <- rep(0, nrow(fitskm))
    ties_ind[which(fitskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(fitskm))}
  fitskm = cbind(fitskm, ties_ind)
  fitskm = fitskm[order(fitskm[,1],fitskm[,3]),]
}

#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param dist String indicating the distribution.
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object for the treatment group.
#'  \item A survfit object for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @importFrom stats na.omit
#' @import survival
#' @importFrom stats cor
#'
#' @export

ROCparametric <- function(time, event, group, method="weibull") {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  if (!is.character(method)) {stop("Argument dist must be a string.")}

  if (method=="cox") {
    coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")
    baseline <- c(forplot[,1], seq(forplot[nrow(forplot),1], 0, length.out = round(nrow(forplot)/4, 0)))
    paramsurv <- cbind(baseline, baseline^exp(coxfit$coefficients))
  }

  if (method == "loglogistic") {
    fit <- survreg(Surv(time, event) ~ group, dist="loglogistic")
    lambda1 = exp(-unname(fit$coefficients[1])/fit$scale)
    beta1 = -unname(fit$coefficients[2])/fit$scale
    alpha1 = 1/fit$scale
    fitsurv1 = (1 + lambda1*exp(beta1*group)*(time^alpha1))^(-1)
    fitskm1 <- getSKM4fit(time, fitsurv1, group)
    forplotfit1 = get4plot(fitskm1)
    newsurv <- seq( min(fitskm1[,1]), 0.00000000001, length.out = round(nrow(forplotfit1)/10, 0) )
    time0 <- ((1 + lambda1*exp(0))*newsurv)^(-alpha1)
    time1 <- ((1 + lambda1*exp(beta1))*newsurv)^(-alpha1)
    newskm1 <- getSKM4fit(time = c(time0, time1), fitsurv = c(newsurv, newsurv),
                          group = c(rep(0, length(newsurv)), rep(1, length(newsurv)))
                          )
    new4plot <- get4plot(newskm1)





  }

  fit2 <- survreg(Surv(time, event) ~ group, dist="weibull")
  lambda2 = exp(-unname(fit2$coefficients[1])/fit2$scale)
  beta2 = -unname(fit2$coefficients[2])/fit2$scale
  alpha2 = 1/fit2$scale
  theta = -unname(fit2$coefficients[2])
  fitsurv2 = base_surv_weibull(lambda2, alpha2, times = time*exp(theta*group))
  fitskm2 <- getSKM4fit(time, fitsurv2, group)
  forplotfit2 = get4plot(fitskm2)





  return(parametricfit = fit)

}
