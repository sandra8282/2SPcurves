
#' ROC when survival goes to 0 for either group
#'
#' @param skm passed from ROCsurv
#' @param silen passed from ROCsurv or btsp
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
#'
#' @keywords internal
#' @noRd

completeROC <- function(skm, silent, xlab, ylab, main, cex.axis,
                        cex.lab, lty, label.inset, label.cex, lwd) {

  forplot = get4plot(skm)

  area <- 0

  if (silent == FALSE) {
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)

    for (k in 2:nrow(forplot)) {
      coord_new = unname(forplot[k-1,])
      coord_new2 = unname(forplot[k,])
      #figure out areas and shading
      if (forplot[k,2]==forplot[k-1,2]) {#move horizontally
        rect(xright = coord_new[1], ytop = coord_new[2],
             xleft = coord_new2[1], ybottom = 0,
             col = "pink", border = "pink")
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {
        if (forplot[k,1]!=forplot[k-1,1] & forplot[k,2]!=forplot[k-1,2]){
          #area and shading for diagonal
          rect(xright = coord_new[1], ytop = coord_new2[2],
               xleft = coord_new2[1], ybottom = 0,
               col = "pink", border = "pink")
          area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
          polygon(x=c(coord_new[1], coord_new[1], coord_new2[1]),
                  y=c(coord_new[2], coord_new2[2], coord_new2[2]),
                  col = "pink", border = "pink")
          area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
          area = area + area_rectang + area_triang
        }
      }
      segments(x0=coord_new[1], y0=coord_new[2],
               x1=coord_new2[1], y1=coord_new2[2], col="black", lwd = lwd)
    }
    abline(c(0,1), col = "black", lty=3, lwd = lwd-0.25)
    area = unname(area)
    text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
         pos=2, cex = 1)

  } else {

    for (k in 2:nrow(forplot)) {
      coord_new = unname(forplot[k-1,])
      coord_new2 = unname(forplot[k,])
      #figure out areas and shading
      if (forplot[k,2]==forplot[k-1,2]) {#move horizontally
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {
        if (forplot[k,1]!=forplot[k-1,1] & forplot[k,2]!=forplot[k-1,2]){
          #area and shading for diagonal
          area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
          area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
          area = area + area_rectang + area_triang
        }
      }
    }
    area = unname(area)
  }

  return(area)
}
