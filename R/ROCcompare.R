getSKM4fit <- function(time, fitsurv, group) {
  fitskm <- cbind(time, fitsurv, group)
  fitskm <- fitskm[!duplicated(fitskm),]
  ties_check <- unique(table(fitskm[,1])) #how many times are repeated
  if (length(ties_check) > 1) {
    ties_times = fitskm[duplicated(fitskm[,1]),1]
    ties_ind <- rep(0, nrow(fitskm))
    ties_ind[which(fitskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(fitskm))}
  fitskm = cbind(fitskm, ties_ind)
  fitskm <- fitskm[order(fitskm[,1], fitskm[,3]),]
  return(fitskm)
}

getmat4cor <- function(forplot, forplotfit){
  all <- matrix(nrow = nrow(forplot), ncol=4)
  all[1,] <- cbind(forplot[1, 1:2], forplotfit[1,])
  for (q in 2:nrow(forplot)) {
    absdiff = (abs(forplotfit[,1]-forplot[q,1]))
    index <- which(absdiff==min(absdiff))
    all[q,] <- cbind(forplot[q,1:2], forplotfit[index[1],])
  }
  return(na.omit(all))
}
#' ROC when survival goes to 0 for either group
#'
#' @param time passed from ROCsurv.
#' @param event passed from ROCsurv.
#' @param group passed from ROCsurv.
#' @param silent passed from ROCsurv.
#' @param abtwc passed from ROCsurv.
#' @param xlab passed from ROCsurv.
#' @param ylab passed from ROCsurv.
#' @param main passed from ROCsurv.
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
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @import survival
#' @importFrom stats cor
#' @importFrom pathmapping CreateMap
#'
#' @keywords internal
#' @noRd

ROCcompare <- function(time, event, group, silent, abtwc, xlab, ylab, main, cex.axis,
                       cex.lab, lty, label.inset, label.cex, lwd) {

  d <- c("lognormal", "loglogistic")
  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  skm <- skm[,-3]
  forplot = get4plot(skm)

  #obtain coefficients and parameters to calculate survival fcts
  lambda = beta = alpha = mu = gamma = sigma = c()
  for (i in 1:2) {
    fit = survreg(Surv(time, event) ~ group, dist=d[i])
    tmu = fit$coefficients[1];  tgamma = fit$coefficients[2]; tsigma = fit$scale;
    mu = c(mu, tmu); gamma = c(gamma, tgamma); sigma = c(sigma, tsigma);
    l = exp(-unname(tmu/tsigma)); b = (-unname(tgamma)/tsigma); a = 1/tsigma;
    lambda = c(lambda, l); beta = c(beta, b); alpha = c(alpha, a);
  }

  #calculate survival fcts
  surv_lognorm = 1 - pnorm((log(time)- mu[1] - gamma[1]*group)/sigma[1])
  surv_loglogis = (1 + lambda[2]*exp(beta[2]*group)*(time^alpha[2]))^(-1)
  fitsurv <- list(surv_lognorm, surv_loglogis)

  forplotfit = list()
  for (i in 1:2) {
    fitskm <- getSKM4fit(time, fitsurv[[i]], group)
    forplotfit[[i]] = get4plot(fitskm)
  }

  forplotfit[[1]][,2] = 1 - pnorm(qnorm(1 - forplotfit[[1]][,1]) - gamma[1]/sigma[1] )
  forplotfit[[2]][,2] = 1 / (1 + (exp(- gamma[2]/sigma[2]) * (1/forplotfit[[2]][,1] - 1)))
  names(forplotfit) <- d

  #Get cox model
  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")

  #correlations, SSR, ABTC
  rho = SSR = areaBTWcurves = rep(0,3)
  cox_surv1 <- forplot[,1]^exp(1/coxfit$coefficients)
  rho[1] <- cor(forplot[,2], cox_surv1)
  resid <- forplot[,2] - cox_surv1
  SSR[1] = sum(resid^2)
  forplot <- cbind(forplot, cox = cox_surv1)
  if (abtwc==TRUE) {
      invisible(capture.output(out <- CreateMap(forplot[,c(1,2)], forplot[,c(1,3)],
                      plotgrid=F, verbose=F, insertopposites=F)))
      areaBTWcurves[1] <- out$deviation
  }


  for (i in 2:3) {
    comparetofit <- getmat4cor(forplot, forplotfit[[i-1]])
    rho[i] <- cor(comparetofit[,2], comparetofit[,4])
    temp <- comparetofit[,2] - comparetofit[,4]
    SSR[i] <- sum(temp^2)
    if (abtwc==TRUE){
      invisible(capture.output(out <- CreateMap(forplot[,c(1,2)], forplotfit[[i-1]],
                                                plotgrid=F, verbose=F, insertopposites=F)))
      areaBTWcurves[i] <- out$deviation

    }
  }

  areaBTWcurves = as.numeric(areaBTWcurves)
  names(rho) = names(SSR) = names(areaBTWcurves) = c("Cox", d)

  bestindex = which(rho == max(rho))
  best = names(rho)[bestindex]
  bestindex2 = which(SSR == min(SSR))
  best2 = names(SSR)[bestindex2]
  bestindex3 = which(areaBTWcurves == min(areaBTWcurves))
  best3 = names(areaBTWcurves)[bestindex3]

  coefficients = rbind(mu, gamma, sigma)
  rownames(coefficients) <- c("(Intercept)", "group", "scale")
  colnames(coefficients) <- d

  if (silent==FALSE){
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    points(forplot[,1], forplot[,2], col = "grey", cex = 1)
    lines(forplot[,1], forplot[,1]^exp(coxfit$coefficients), lty=lty[1], lwd = lwd)
    lines(forplotfit[[1]][,1], forplotfit[[1]][,2], lty=lty[2], lwd = lwd)
    lines(forplotfit[[2]][,1], forplotfit[[2]][,2], lty=lty[3], lwd = lwd)
    abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)
    legend("bottomright", legend= c("Cox PHM", "Lognormal", "Loglogistic"),
           lty = lty, cex=label.cex, xjust = 1, yjust = 0, bg = "white",
           bty='n', seg.len = 0.8, x.intersp=0.9, y.intersp = 0.85)
  }

  if (abtwc==TRUE) {
    return(list(HR = exp(coxfit$coefficients),
                BestRHO = best, BestSSR = best2, BestAreaBTWcurves = best3,
                RHO = rho, SSR = SSR, areaBTWcurves = areaBTWcurves,
                coefficients = coefficients))
  } else {return(list(KMres = KMres,  HR = exp(coxfit$coefficients),
                      BestRHO = best, BestSSR = best2,
                      RHO = rho, SSR = SSR, coefficients = coefficients))}


}
