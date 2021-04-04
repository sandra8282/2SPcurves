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
#' @import survival
#' @importFrom stats cor
#' @importFrom pathmapping CreateMap
#' @importFrom utils capture.output
#' @importFrom data.table data.table
#' @importFrom data.table setkeyv
#' @importFrom data.table setorder
#' @importFrom data.table :=
#'
#' @keywords internal
#' @noRd

ROCandPHM <- function(time, event, group, silent, abtwc, xlab, ylab, main, cex.axis,
                      cex.lab, lty, label.inset, label.cex, lwd) {

  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")
  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  #correlations, SSR and area between curves
  cx = c(seq(0, 1, 0.0001), 1)
  cy <- cx^(exp(coxfit$coefficients)) #surv1 = surv0^HR

  a=data.table(cx, cy)
  a[,merge:=cx]
  b=data.table(forplot)
  b[,merge:=forplot[,1]]
  setkeyv(a,c('merge'))
  setkeyv(b,c('merge'))
  MergedForplot=a[b,roll='nearest']
  setorder(MergedForplot)

  rho <- cor(MergedForplot$cy, MergedForplot$y)
  resid <-  MergedForplot$cy - MergedForplot$y
  SSR <- sum(resid^2)
  if (abtwc == TRUE){
      invisible(capture.output(out <- CreateMap(MergedForplot[,c(1,2)],
                                                MergedForplot[,c(4,5)],
                                            plotgrid=F, verbose=F, insertopposites=F)))
      areaBTWcurves <- out$deviation
  }

  if(silent==FALSE){
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    lines(forplot[,1], forplot[,2], col="black", lty=lty[1], lwd = lwd)
    lines(cx, cy, lty=lty[2], lwd = lwd)
    abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)
    legend("topleft", c("KM-Based Curve", "Cox-Based Curve"), lty = lty,
           inset=label.inset, cex=label.cex, bg = "white", bty='n', seg.len = 0.8,
           x.intersp=0.9, y.intersp = 0.85, lwd = lwd)

  text(x=0.99, y=0.25,
       labels = bquote(hat(rho) == .(round(rho, 3))),
       pos=2, cex=label.cex)
  text(x=0.99, y=0.15,
       labels = paste("SSR = ", round(SSR, 3), sep=""),
       pos=2, cex=label.cex)
    if (abtwc == TRUE){
      text(x=0.99, y=0.05,
           labels = paste("ABTC = ", round(areaBTWcurves, 3), sep=""),
           pos=2, cex=label.cex)
    }
  }

  if (abtwc == TRUE){
    res <- list(SSR = SSR, rho = rho, areaBTWcurves = areaBTWcurves)
    } else {res <- list(SSR = SSR, rho = rho)}

  return(res)

}

