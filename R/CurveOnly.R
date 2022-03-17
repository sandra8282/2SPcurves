#' Plot an ROC curve
#'
#' @param skm passed from ROCsurv
#' @param xlab passed from ROCsurv.
#' @param ylab passed from ROCsurv.
#' @param main passed from ROCsurv.
#'
#' @return A plot of the ROC curve
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @keywords internal
#' @noRd

onlyROC <- function(skm, xlab, ylab, main, cex.axis,
                    cex.lab, lty, label.inset, label.cex, lwd) {

  forplot <- get4plot(skm)
  minx <- min(forplot[,1])
  miny <- min(forplot[,2])

  plot(forplot, type = "l", lty = lty[1], lwd = lwd, #las=1
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
  abline(c(0,1), col = "darkgrey", lty=1, lwd = lwd-0.25)

  colnames(forplot) <- c("u", "R_u", "tienext")
  return(R_u = forplot)

}
