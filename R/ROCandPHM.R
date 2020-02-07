#' ROC when survival goes to 0 for either group
#'
#' @param skm passed from ROCsurv
#' @param silen passed from ROCsurv or btsp
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
#'
#' @keywords internal
#' @noRd

ROCandPHM <- function(time, event, group) {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")

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
  colnames(forplot) <- c("x", "y")
  area = areaC = sum_sqrres = 0

  plot(NULL, type="n", xlab="", ylab="", las=1,
       xlim=c(0,1), ylim = c(0, 1)) #to make tight axis: xaxs="i", yaxs="i"
  title(main="ROC", xlab="Control Group Survival",
          ylab="Treatment Group Survival",
          cex.main = 1)

  points(forplot[,1], forplot[,2])
  lines(forplot[,1], forplot[,1]^exp(coxfit$coefficients), col="blue")
  abline(c(0,1), col = "red", lty=2)

  #correlations and SSR
  HRcheck <- cor(forplot[,2], forplot[,1]^exp(coxfit$coefficients))
  SSR <- sum((forplot[,2] - forplot[,1]^exp(coxfit$coefficients))^2)

  text(x=0.99, y=0.1,
       labels = paste("rho = ", round(HRcheck, 4), sep=""),
       pos=2)

  return(list(KMres = KMres, SSR = SSR, rho = HRcheck))

}
