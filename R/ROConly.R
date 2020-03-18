#' Use incomplete ROC to compute AUC by restricting space
#'
#' @param skm passed from ROCsurv
#'
#' @return A plot of the ROC curve
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @keywords internal
#' @noRd

onlyROC <- function(skm) {

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
  minx <- min(forplot[,1])
  miny <- min(forplot[,2])
  area <- 0

  plot(NULL, type="n", xlab="", ylab="", las=1,
       xlim=c(0,1), ylim = c(0, 1)) #to make tight axis: xaxs="i", yaxs="i"
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
        cex.main = 1)

  for (k in 2:nrow(forplot)) {
      coord_new = unname(forplot[k-1,])
      coord_new2 = unname(forplot[k,])
      segments(x0=coord_new[1], y0=coord_new[2],
               x1=coord_new2[1], y1=coord_new2[2], col="black")
    }

  abline(c(0,1), col = "red", lty=3)

}
