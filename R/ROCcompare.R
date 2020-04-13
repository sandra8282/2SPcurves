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
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @import survival
#' @importFrom stats cor
#'
#' @keywords internal
#' @noRd

ROCcompare <- function(time, event, group) {

  d <- c("lognormal", "loglogistic")
  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
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

  plot(NULL, type="n", las=1,
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab="Control Group Survival", ylab="Treatment Group Survival",
       cex.axis = 1.25, cex.lab = 1.25)
  points(forplot[,1], forplot[,2], col = "grey50", cex = 0.75)
  lines(forplot[,1], forplot[,1]^exp(coxfit$coefficients), col="darkred")
  lines(forplotfit[[1]][,1], forplotfit[[1]][,2], col="darkblue")
  lines(forplotfit[[2]][,1], forplotfit[[2]][,2], col="darkgreen")
  abline(c(0,1), col = "grey", lty=2)
  legend("bottomright", legend= c("Cox PHM", "Lognormal", "Loglogistic"),
         lty = 1, col=c("darkred", "darkblue", "darkgreen"),
         cex=1, bty = "n", xjust = 1, yjust = 0, y.intersp = 1)

  #correlations and SSR
  rho = SSR = areaBTWcurves = rep(0,3)
  cox_surv1 <- forplot[,1]^exp(coxfit$coefficients)
  rho[1] <- cor(forplot[,2], cox_surv1)
  resid <- forplot[,2] - cox_surv1
  SSR[1] = sum(resid^2)
  forplot <- cbind(forplot, cox = cox_surv1)
  invisible(capture.output(out <- pathmapping::CreateMap(forplot[,c(1,2)], forplot[,c(1,3)],
                                                         plotgrid=F, verbose=F, insertopposites=F)))
  areaBTWcurves[1] <- out$deviation

  for (i in 2:3) {
    comparetofit <- getmat4cor(forplot, forplotfit[[i-1]])
    rho[i] <- cor(comparetofit[,2], comparetofit[,4])
    temp <- comparetofit[,2] - comparetofit[,4]
    SSR[i] <- sum(temp^2)
    invisible(capture.output(
      out <- pathmapping::CreateMap(forplot[,c(1,2)], forplotfit[[i-1]],
                                    plotgrid=F, verbose=F, insertopposites=F)))
    areaBTWcurves[i] <- out$deviation
  }
  areaBTWcurves = as.numeric(areaBTWcurves)
  names(rho) = names(SSR) = names(areaBTWcurves) = c("Cox", d)
  # plot(resid[,1], col = "darkred", pch = 20)
  # points(resid[,2], col = "darkblue", pch = 20)
  # points(resid[,3], col = "darkgreen", pch = 20)
  # abline(h=0, col="black", lty = 2)
  # legend("topleft", legend= c("Cox PHM", "Lognormal", "Loglogistic"),
  #        col=c("darkred", "darkblue", "darkgreen"), pch = 20,
  #        cex=1, bty = "n", xjust = 1, yjust = 0, y.intersp = 1)

  #find best
  bestindex = which(rho == max(rho))
  best = names(rho)[bestindex]
  bestindex2 = which(SSR == min(SSR))
  best2 = names(SSR)[bestindex2]
  bestindex3 = which(areaBTWcurves == min(areaBTWcurves))
  best3 = names(areaBTWcurves)[bestindex3]

  coefficients = rbind(mu, gamma, sigma)
  rownames(coefficients) <- c("(Intercept)", "group", "scale")
  colnames(coefficients) <- d

  return(list(KMres = KMres,  HR = exp(coxfit$coefficients),
              BestRHO = best, BestSSR = best2, BestAreaBTWcurves = best3,
              RHO = rho, SSR = SSR, areaBTWcurves = areaBTWcurves,
              coefficients = coefficients))

}
