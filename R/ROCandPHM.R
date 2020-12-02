#' ROC when survival goes to 0 for either group
#'
#' @param time passed from ROCsurv.
#' @param event passed from ROCsurv.
#' @param group passed from ROCsurv.
#' @param silent passed from ROCsurv.
#' @param abtwc passed from ROCsurv.
#' @param xlab passed from ROCsurv
#' @param ylab passed from ROCsurv
#' @param main passed from ROCsurv
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
#' @importFrom pathmapping CreateMap
#' @importFrom utils capture.output
#'
#' @keywords internal
#' @noRd

ROCandPHM <- function(time, event, group, silent, abtwc, xlab, ylab, main, cex.axis,
                      cex.lab, lty, label.inset, label.cex, lwd) {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")

  #correlations, SSR and area between curves
  cox_surv1 <- forplot[,1]^exp(coxfit$coefficients) #surv1 = surv0^HR
  rho <- cor(forplot[,2], cox_surv1)
  resid <- forplot[,2] - cox_surv1
  SSR <- sum(resid^2)
  forplot <- cbind(forplot, cox = cox_surv1)
  if (abtwc == TRUE){
      invisible(capture.output(out <- CreateMap(forplot[,c(1,2)], forplot[,c(1,3)],
                                            plotgrid=F, verbose=F, insertopposites=F)))
      areaBTWcurves <- out$deviation
  }

  if(silent==FALSE){
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
  lines(forplot[,1], forplot[,2], col="black", lty=lty[1], lwd = lwd)
  lines(forplot[,1], forplot[,1]^exp(coxfit$coefficients), lty=lty[2], lwd = lwd)
  abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)
  legend("topleft", c("KMROC", "Cox ROC"), lty = lty,
         inset=label.inset, cex=label.cex, bg = "white", bty='n', seg.len = 0.8,
         x.intersp=0.9, y.intersp = 0.85, lwd = lwd)

  text(x=0.99, y=0.25,
       labels = bquote(hat(rho) == .(round(rho, 4))),
       pos=2, cex=label.cex)
  text(x=0.99, y=0.15,
       labels = paste("SSR = ", round(SSR, 4), sep=""),
       pos=2, cex=label.cex)
    if (abtwc == TRUE){
      text(x=0.99, y=0.05,
           labels = paste("ABTC = ", round(areaBTWcurves, 4), sep=""),
           pos=2, cex=label.cex)
    }


  }

  if (abtwc == TRUE){
    res <- list(KMres = KMres, SSR = SSR, rho = rho, areaBTWcurves = areaBTWcurves)
    } else {res <- list(KMres = KMres, SSR = SSR, rho = rho)}

  return(res)

}

