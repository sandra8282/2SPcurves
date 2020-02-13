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

ROCparametric <- function(time, event, group, dist="Weibull") {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  dat <- cbind(time, event, group)
  if (!is.character(dist)) {stop("Argument dist must be a string.")}

  fit <- survreg(Surv(time, event) ~ group, dist=dist)
  lambda = exp(-unname(fit$coefficients[1])/fit$scale)
  beta = -unname(fit$coefficients[2])/fit$scale
  alpha = 1/fit$scale

  if (dist=="loglogistic") {
    fitsurv = (1 + lambda*exp(beta*group)*(time^alpha))^(-1)
  } else {
    theta = -unname(fit$coefficients[2])
    fitsurv = base_surv_weibull(lambda, alpha, times = time*exp(theta*group))
  }

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
  fitskm = fitskm[order(fitskm[,1], fitskm[,3]),]

  forplotfit = get4plot(fitskm)

  sum_sqrres = 0

  plot(NULL, type="n", xlab="", ylab="", las=1,
       xlim=c(0,1), ylim = c(0, 1)) #to make tight axis: xaxs="i", yaxs="i"
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
        cex.main = 1)

  forplotfit <- forplotfit[-which(forplotfit[,1]<min(forplot[,1])),]
  points(forplot[,1], forplot[,2])
  lines(forplotfit[,1], forplotfit[,2], col="blue")
  abline(c(0,1), col = "red", lty=2)

  #correlations and SSR
  text(x=0.99, y=0.4,
       labels = dist,
       pos=2)
  text(x=0.99, y=0.3,
       labels = paste("lambda = ", round(lambda, 2), sep=""),
       pos=2)
  text(x=0.99, y=0.2,
       labels = paste("beta = ", round(beta, 2), sep=""),
       pos=2)
  text(x=0.99, y=0.1,
       labels = paste("alpha = ", round(alpha, 2), sep=""),
       pos=2)

  return(parametricfit = fit)

}
