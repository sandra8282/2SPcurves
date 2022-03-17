#' Check Proportional Hazards for a given model
#'
#' @param time Time to event or censoring.
#' @param event An indicator vector with values of 1 for individuals who had the event occur or 0 if the participant was censored.
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param changepoints A vector with times where we will cut as part of a time-varying covariate.
#' @param silent Logical argument, FALSE indicates the user wants plots and TRUE indicates no plots only AUC calculations (default is FALSE).
#' @param xlab String argument for the horizontal axis label of the ROC curve.
#' @param ylab String argument for the vertical axis label of the ROC curve.
#' @param main String argument for the title of the ROC curve.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling the inset of the legend.
#' @param legend.cex Optional graphical parameter for magnification of the legend's text.
#' @param lty Optional graphical parameter to set the type of line to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.

#'
#' @return A plot of the curve and object containing:
#' \itemize{
#'  \item A cox model object based on the change points.
#'  \item The area under the curve for the given plot.
#' }
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @importFrom stats na.omit
#' @import survival
#' @importFrom stats cor
#' @importFrom stats approxfun
#' @importFrom stats integrate
#' @importFrom pathmapping CreateMap
#' @importFrom utils capture.output
#' @importFrom stats coef
#' @importFrom stats formula
#'
#' @export

CPCM2SSP <- function(time, event, group, changepoints, silent, xlab="", ylab="", main="", cex.axis,
                      cex.lab, lty, legend.inset, legend.cex, lwd){
  label.inset = legend.inset; label.cex=legend.cex;

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot1 <- get4plot(skm); colnames(forplot1) = c("u", "R(u)", "tienext")
  forplot <- data.frame(time = c(0, skm[-nrow(skm),1]), forplot1)
  colnames(forplot) <- c("t", "x", "y", "tienext")

  dat <- data.frame(time, event, group)
  datnew <- survSplit(dat, cut=changepoints, end="time", start="start",
                      event="event", episode="Period")
  fit1 <- coxph(Surv(time, event) ~ group, data = dat)
  fit2 <- coxph(Surv(start, time, event) ~ group:strata(Period), data=datnew)
  HRScsum <- exp(cumsum(coef(fit2)))
  v <- c(0, changepoints, max(time))
  forplot$timeind1 <-  findInterval(forplot$t, v)
  forplot$coxy <- forplot$x^HRScsum[forplot$timeind1]
  forplot <- cbind(forplot, cy=forplot$x^exp(coef(fit1)))

    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    lines(forplot$x, forplot$coxy, col = "black", lty=lty[1], lwd = lwd)
    #lines(forplot$x, forplot$cy, col = "black", lty=lty[1], lwd = lwd)
    lines(forplot$x, forplot$y, col="black", lty=lty[2], lwd = lwd)
    abline(c(0,1), col = "darkgrey", lty=1, lwd = lwd-0.25)

    legend("topleft", c("KM-Based Curve", "Change-point Cox Model"), lty = c(lty[2], lty[1]),
           inset=label.inset, cex=label.cex, bg = "white", bty='n', seg.len = 0.8,
           x.intersp=0.9, y.intersp = 0.85, lwd = lwd)

  res <- list(R_u = forplot1, CPCoxFit = fit2, Coxfit = fit1)
  return(res)

}

