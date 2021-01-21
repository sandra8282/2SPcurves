#' Use incomplete ROC to compute AUC by restricting space
#'
#' @param skm passed from ROCsurv
#' @param silent passed from ROCsurv or btsp
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
#' @keywords internal
#' @noRd

restrictROC <- function(skm, silent, xlab, ylab, main, cex.axis,
                        cex.lab, lty, label.inset, label.cex, lwd) {

  forplot <- get4plot(skm)
  minx <- forplot[nrow(forplot),1]
  area <- 0

  area <- completeROC(skm, silent, xlab, ylab, main, cex.axis,
                      cex.lab, lty, label.inset, label.cex, lwd)

  mina = 0.5*(1 - minx^2)
  maxa = 1-minx
  area_adj = unname(0.5*(1 + (area - mina)/(maxa - mina)))

  return(c(rAUC_unadj = area, rAUC_adj = unname(area_adj)))
}
