base_surv_weibull <- function(lambda, alpha, times){
  s0 <- exp(-lambda*(times)^alpha)
  return(s0)
}

getSKM4fit <- function(time, fitsurv, group) {
  fitskm <- cbind(time, fitsurv, group)
  fitskm <- fitskm[!duplicated(fitskm),]
  ties_check <- unique(table(fitskm[,1]))
  if (length(ties_check) > 1) {
    ties_times = as.integer(names(which(table(fitskm[,1])>1)))
    ties_ind <- rep(NA, nrow(fitskm))
    ties_ind[which(fitskm[,1] %in% ties_times)]=1
    ties_ind[-which(fitskm[,1] %in% ties_times)]=0
  } else {ties_ind <- rep(0, nrow(fitskm))}
  fitskm = cbind(fitskm, ties_ind)
  fitskm = fitskm[order(fitskm[,1],fitskm[,3]),]
}
#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
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
#' @keywords internal
#' @noRd

ROCcompare <- function(time, event, group) {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")

  fit <- survreg(Surv(time, event) ~ group, dist="loglogistic")
  lambda1 = exp(-unname(fit$coefficients[1])/fit$scale)
  beta1 = -unname(fit$coefficients[2])/fit$scale
  alpha1 = 1/fit$scale
  fitsurv1 = (1 + lambda*exp(beta*group)*(time^alpha))^(-1)
  fitskm1 <- getSKM4fit(time, fitsurv1, group)
  forplotfit1 = get4plot(fitskm1)

  fit2 <- survreg(Surv(time, event) ~ group, dist="weibull")
  lambda2 = exp(-unname(fit2$coefficients[1])/fit2$scale)
  beta2 = -unname(fit2$coefficients[2])/fit2$scale
  alpha2 = 1/fit2$scale
  theta = -unname(fit$coefficients[2])
  fitsurv2 = base_surv_weibull(lambda, alpha, times = time*exp(theta*group))
  fitskm2 <- getSKM4fit(time, fitsurv2, group)
  forplotfit2 = get4plot(fitskm2)

  sum_sqrres = 0

  plot(NULL, type="n", xlab="", ylab="", las=1,
       xlim=c(0,1), ylim = c(0, 1)) #to make tight axis: xaxs="i", yaxs="i"
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
        cex.main = 1)

  points(forplot[,1], forplot[,2])
  lines(forplot[,1], forplot[,1]^exp(coxfit$coefficients), col="blue")
  lines(forplotfit1[,1], forplotfit1[,2], col="red")
  lines(forplotfit2[,1], forplotfit2[,2], col="green")
  abline(c(0,1), col = "red", lty=2)

  legend("bottomright", legend=c("Cox PH", "Log-logistic", "Weibull"), lty = 1,
         col=c("blue", "red", "green"),
         cex=0.9, bty = "n", xjust = 1, yjust = 0, y.intersp = 0.9)

  #correlations and SSR

  text(x=0.99, y=0.1,
       labels = paste("rho = ", round(HRcheck, 4), sep=""),
       pos=2)

  return(list(KMres = KMres, SSR = SSR, rho = HRcheck))

}
