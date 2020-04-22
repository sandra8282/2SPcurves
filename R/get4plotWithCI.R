#' #' ROC ordered sample ready for plotting
#' #'
#' #' @param skm matrix passed from other functions.
#' #'
#' #' @return A matrix with x y coordinates ready for plotting ROC points
#' #'
#' #' @importFrom stats na.omit
#' #'
#' #' @keywords internal
#' #' @noRd
#' #'
#' #'
#' #'
#' get4plotWithCI <- function(skm) {
#'   x=y=c(1, rep(NA, nrow(skm)))
#'   err = c(0, skm[,3]);
#'   i=1
#'   while (i <= nrow(skm)) {
#'     if (skm[i,4]==0 & skm[i,5]==0) {#placebo no ties
#'       # horizontal move
#'       x[i+1] = skm[i,2];
#'       xlow[i+1] = xup[i+1] = x[i];
#'       # vertical stay
#'       if (is.na(y[i])) {
#'         y[i+1] = y[i-1];
#'         ylow[i+1] =  y[i-1] - err[i]; yup[i+1] =  y[i-1] + err[i];
#'       } else {
#'         y[i+1] = y[i];
#'         ylow[i+1] =  y[i] - err[i]; yup[i+1] =  y[i] + err[i];
#'       }
#'       i=i+1
#'     } else {
#'       if (skm[i,4]==1 & skm[i,5]==0) {#drug no ties
#'         #vertical move
#'         y[i+1] = skm[i,2];
#'         ylow[i+1] = yup[i] = y[i];
#'         #horizontal stay
#'         if (is.na(x[i])) {
#'           x[i+1] = x[i-1];
#'           xlow[i+1] =  x[i-1] - err[i]; xup[i+1] =  x[i-1] + err[i];
#'         } else {
#'           x[i+1] = x[i];
#'           xlow[i+1] =  x[i] - err[i]; xup[i+1] =  x[i] + err[i];
#'         }
#'         i=i+1
#'       } else {
#'         if (skm[i,5]==1) {#tie
#'           if (skm[i,4]==0){#horizontal move first
#'             x[i+1] = skm[i,2];
#'             y[i+1] = skm[i+1,2];
#'           } else{#vertical move first
#'             x[i+1] = skm[i,2];
#'             y[i+1] = skm[i+1,2];
#'           }
#'           i=i+2
#'         }
#'       }
#'     }
#'   }
#'
#'   forplot <- na.omit(cbind(x, y, xlow, xup = ifelse(xup > 1, 1, xup),
#'                            ylow, yup = ifelse(yup > 1, 1, yup)))
#'   colnames(forplot) <- c("x", "y", "xlow", "xup", "ylow", "yup")
#'   return(forplot)
#' }
