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
#' @param AUC Numeric value for the size of the AUC = P(T1 >  T0) that is expected. Must be between 0.55 and 1.0.
#' @param HR Only needed is AUC cannot be provided. Numeric value for the size of the hazard ratio.
#' @param sig.level Numeric value for the type 1 Error Rate. Default is 0.05.
#' @param power Numeric value for the target power. Must be between 0.2 and 1. Default is 0.8.
#' @param alternative Character string. Options: "two.sided" or "one.sided".
#' @param method Character string for method to be used to obtain sample size or power. Options: "parametric binormal", "nonparametric"
#' @param type Character string. Options: "biomarker" or "survival".
#' @param variance Required only when method = "parametric binormal". Numeric value for the variance to use when calculating sample size and power.
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

ROCSSnonparam <- function(AUC, HR, ctrl.event.rate, trt.event.rate, cens.rate=0, sig.level=0.05, power=0.8,
                          alternative = "two.sided", method="survival"){

  Za <- ifelse(alternative == "two.sided", qnorm(1-sig.level/2), qnorm(1-sig.level))
  Zb <- qnorm(power)

  if (method=="survival"){
    #Survival Specific
    if (missing(AUC)) {
      if (missing(HR)) {stop("Must provide either AUC or HR.")} else {AUC=HR/(HR+1)}
      } else {HR = AUC/(1-AUC)}
    if (missing(ctrl.event.rate)) {stop("Must provide ctrl.event.rate")}
    p0tinv <- 1/ ctrl.event.rate
    if (missing(trt.event.rate)) {p1tinv <- 1/(ctrl.event.rate/HR)
    } else{ p1tinv <- 1/trt.event.rate}
    sqrt1 <- sqrt(p0tinv/8)
    sqrt2 <- sqrt((exp(2*log(HR))*(p1tinv + p0tinv))/((1+HR)^4))
    denom <- (1-cens.rate)*((AUC-0.5)^2)
    res <- (Za*sqrt1 + Zb*sqrt2)^2 / denom
  } else {
    #Biomarker/General as in Obuchowski paper
    res <- uniroot.all(fnpBiomark, c(0, 1.79e+300), maxiter = 10000, tol = 0.00001,
                       Za = Za, Zb = Zb, theta1=AUC)}
  return(res)
}

