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

  d <- c("lognormal", "loglogistic", "spline")
  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  #obtain coefficients and parameters to calculate AFT survival fcts
  lambda = beta = alpha = mu = gamma = sigma = c()
  for (i in 1:2) {
    fit = survreg(Surv(time, event) ~ group, dist=d[i])
    tmu = fit$coefficients[1];  tgamma = fit$coefficients[2]; tsigma = fit$scale;
    mu = c(mu, tmu); gamma = c(gamma, tgamma); sigma = c(sigma, tsigma);
    l = exp(-unname(tmu/tsigma)); b = (-unname(tgamma)/tsigma); a = 1/tsigma;
    lambda = c(lambda, l); beta = c(beta, b); alpha = c(alpha, a);
  }

  #Get cox model
  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")
  #calculate survival fcts

  u = forplot[,1]
  forplotfit = data.frame(u, u^exp(coxfit$coefficients),
                           1 - pnorm(qnorm(1 - u) - gamma[1]/sigma[1]),
                           1 / (1 + (exp(- gamma[2]/sigma[2]) * (1/u - 1))))
  names(forplotfit) =  names(forplotfit2) = c("u", "Cox", d)

  for (i in 1:3) {
    rho[i] <- cor(forplotfit[,i+1], forplot[,2])
    temp <- forplot[,2] - forplotfit[,i+1]
    SSR[i] <- sum(temp^2)
    if (abtwc==TRUE){
      invisible(capture.output(out <- CreateMap(forplotfit[, c(1,forplotfit[,i+1])],
                                                forplot,
                                                plotgrid=F, verbose=F, insertopposites=F)))
      areaBTWcurves[i] <- out$deviation

    }
  }

  areaBTWcurves = as.numeric(areaBTWcurves)
  names(rho) = names(SSR) = names(areaBTWcurves) = c("Cox", d)


  u = c(seq(0, 1, 0.0001), 1)
  forplotfit2 = data.frame(u, u^exp(coxfit$coefficients),
                          1 - pnorm(qnorm(1 - u) - gamma[1]/sigma[1]),
                          1 / (1 + (exp(- gamma[2]/sigma[2]) * (1/u - 1))))

  #correlations, SSR, ABTC
  rho = SSR = areaBTWcurves = rep(0,3)

  # a=data.table(forplotfit)
  # a[,merge:=u]
  # b=data.table(forplot)
  # b[,merge:=forplot[,1]]
  # setkeyv(a,c('merge'))
  # setkeyv(b,c('merge'))
  # MergedForplot=a[b,roll='nearest']
  # setorder(MergedForplot, -x)
  # MergedForplot = as.data.frame(MergedForplot)

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
    lines(forplotfit[,1], forplotfit[,2], lty=lty[1], lwd = lwd, col = "darkred")
    lines(forplotfit[,1], forplotfit[,3], lty=lty[2], lwd = lwd, col = "blue")
    lines(forplotfit[,1], forplotfit[,4], lty=lty[3], lwd = lwd, col = "darkgreen")
    abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)
    legend("bottomright", legend= c("Cox PHM", "Lognormal", "Loglogistic"),
           lty = lty, col = c("darkred", "blue", "darkgreen"), cex=label.cex,
           xjust = 1, yjust = 0, bg = "white", bty='n', seg.len = 0.8,
           x.intersp=0.9, y.intersp = 0.85)
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
