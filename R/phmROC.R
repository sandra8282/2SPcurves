#' Complete the ROC using proportional hazards model
#'
#' @param skm matrix passed from ROCsurv
#' @param phm a cox proportional hazard model passed from ROCsurv
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object for the treatment group.
#'  \item A survfit object for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @import survival
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @keywords internal
#' @noRd

ph_loglogROC <- function(skm, phm) {



  plot(c(0,1), c(0, 1), type="n", xlab="", ylab="")
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
        cex.main = 1)

  x=y= c(1, rep(NA, nrow(skm)))

  i=1
  while (i <= nrow(skm)) {
    if (skm[i,3]==0 & skm[i,4]==0) {#placebo no ties
      # horizontal move
      x[i+1] = skm[i,2]
      # vertical stay
      if (is.na(y[i])) {y[i+1] = y[i-1]
      } else {y[i+1] = y[i]}
      i=i+1
    } else {
      if (skm[i,3]==1 & skm[i,4]==0) {#drug no ties
        #vertical move
        y[i+1] = skm[i,2]
        #horizontal stay
        if (is.na(x[i])) {x[i+1] = x[i-1]
        } else {x[i+1] = x[i]}
        i=i+1
      } else {
        if (skm[i,4]==1) {#tie
          if (skm[i,2]==0){#horizontal move
            x[i+1] = skm[i,2]
            y[i+1] = skm[i+1,2]
          } else{#vertical move
            x[i+1] = skm[i,2]
            y[i+1] = skm[i+1,2]
          }
          i=i+2
        }
      }
    }
  }

  forplot <- na.omit(cbind(x,y))
  area <- 0

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
               x1=coord_new2[1], y1=coord_new2[2], col="black")
  }
  abline(c(0,1), col = "black", lty=2)
  area = unname(area)
  text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
       pos=2, cex = 1)

  return(area)
}
