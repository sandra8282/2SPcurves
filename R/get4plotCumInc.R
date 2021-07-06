#' 2SSP ordered sample ready for plotting
#'
#' @param skm matrix passed from other functions.
#'
#' @return A matrix with x y coordinates ready for plotting ROC points
#'
#' @importFrom stats na.omit
#'
#' @keywords internal
#' @noRd
#'
#'
get4plotCumInc <- function(skm) {
  x=y=tie=c(0, rep(NA, nrow(skm)));
  i=1
  while (i < nrow(skm)) {
    if (skm[i,3]==0 & skm[i,4]==0) {#placebo no ties
      # horizontal move
      x[i+1] = skm[i,2]
      tie[i+1] = 0
      # vertical stay
      if (is.na(y[i])) {y[i+1] = y[i-1]
      } else {y[i+1] = y[i]}
      i=i+1
    } else {
      if (skm[i,3]==1 & skm[i,4]==0) {#drug no ties
        #vertical move
        y[i+1] = skm[i,2]
        tie[i+1] = 0
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
  tie <- ifelse(is.na(tie), 1, tie)
  forplot <- na.omit(cbind(x,y, c(tie[-1],0)))
  colnames(forplot) <- c("x", "y", "tienext")
  return(forplot)
}
