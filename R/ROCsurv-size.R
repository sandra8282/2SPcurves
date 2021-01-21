fnpBiomark <- function(n, Za, Zb, theta1){
  lt <- theta1*(1-theta1)
  bt <- 2*(theta1^2)/(1+theta1) + theta1/(2 - theta1) - 2*theta1^2
  trm1 <- (Za)*sqrt((2*n + 1)/(12*n^2))
  trm2 <- (Zb)*sqrt((lt + (n-1)*bt)/(n^2))
  (trm1 + trm2)^2 - (theta1-0.5)^2
}

#' Nonparametric sample size or power based on the Area Under the ROC Curve.
#'
#' @description
#' This function finds the needed sample size to .
#'
#'AUC, sig.level, power, variance, alternative = "two.sided"
#' @param AUC numeric value for the size of the AUC = P(T1 >  T0) that is expected. Must be between 0.55 and 1.0.
#' @param sig.level Numeric value for the type 1 Error Rate. Default is 0.05.
#' @param power Numeric value for the target power. Must be between 0.2 and 1. Default is 0.8.
#' @param alternative Character string. Options: "two.sided" or "one.sided".
#' @param method Character string. Options: "biomarker" or "survival".
#' @param control.prop The proportion of participants in the control group.
#' @param event.rate Anticipated event rate.
#' @param cens.rate Anticipated censoring rate.
#'
#' @return The sample size for each of the two trial arms (n) or power.
#'
#' @details
#' Methods.
#'
#' @importFrom rootSolve uniroot.all
#'
#' @export
#'

ROCsurvSize <- function(AUC, control.prop, event.rate, cens.rate=0, sig.level=0.05, power=0.8,
                          alternative = "two.sided", method="survival"){

  Za <- ifelse(alternative == "two.sided", qnorm(1-sig.level/2), qnorm(1-sig.level))
  Zb <- qnorm(power)

  if (method=="survival"){
    #Survival Specific
    if (missing(AUC)) {stop("Must provide an AUC for the alternate hypothesis.")
      } else {HR = AUC/(1-AUC)}
    if (missing(control.prop)) {stop("Must provide control.prop")}
    if (missing(event.rate)) {stop("Must provide event.rate")}

    P0 <- control.prop; P1 <- 1-control.prop; d <- event.rate;

    numerator <- Za/4 + Zb*sqrt((exp(2*log(HR)))/((1+HR)^4))
    denom <- 2*(1-cens.rate)*((AUC-0.5)^2)*d*P0*P1
    res <- ((numerator)^2) / denom
  } else {
    #Biomarker/General based on Hanley Variance
    res <- uniroot.all(fnpBiomark, c(0, 1.79e+300), maxiter = 10000, tol = 0.00001,
                       Za = Za, Zb = Zb, theta1=AUC)}
  return(res)
}


