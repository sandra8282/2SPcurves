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
  all = matrix(nrow=nrow(forplotCOX), ncol=4)
  all[1,] = rep(1,4)
  colnames(all) <- c("COXs0", "COXs1", "KMs0", "KMs1")

  plot(NULL, type="n", xlab="", ylab="", las=1,
       xlim=c(0,1), ylim = c(0, 1)) #to make tight axis: xaxs="i", yaxs="i"
  title(main="ROC", xlab="Control Group Survival",
          ylab="Treatment Group Survival",
          cex.main = 1)

  for (k in 2:nrow(forplot)) {

    #### FROM KM CURVES ####
    coord_new = unname(forplot[k-1,])
    coord_new2 = unname(forplot[k,])
    # segments(x0=coord_new[1], y0=coord_new[2],
    #          x1=coord_new2[1], y1=coord_new2[2], col="black")
      #figure out areas
          if (forplot[k,2]==forplot[k-1,2]) {
            #move was horizontal
            area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
          } else if (forplot[k,1]!=forplot[k-1,1] & forplot[k,2]!=forplot[k-1,2]){
            #move was diagonal
              area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
              area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
              area = area + area_rectang + area_triang
          }

    #### FROM COX MODEL FITTED CURVES ####
    coord_new_C = unname(forplotCOX[k-1,])
    coord_new2_C = unname(forplotCOX[k,])
    # segments(x0=coord_new_C[1], y0=coord_new_C[2],
    #          x1=coord_new2_C[1], y1=coord_new2_C[2], col="blue")
      #figure out areas
          if (forplotCOX[k,2]==forplotCOX[k-1,2]) {
            #move was horizontal
            areaC = areaC + (coord_new_C[1] - coord_new2_C[1])*(coord_new_C[2])
          } else if (forplotCOX[k,1]!=forplotCOX[k-1,1] & forplotCOX[k,2]!=forplotCOX[k-1,2]){
            #move was diagonal
            area_rectang = (coord_new_C[1] - coord_new2_C[1])*(coord_new2_C[2])
            area_triang = 0.5 * (coord_new_C[1] - coord_new2_C[1]) * (coord_new_C[2] - coord_new2_C[2])
            areaC = areaC + area_rectang + area_triang
          }

    #even out sets
    if (forplot[k, 1]!=forplot[k-1,1]) {
      #there was horizontal or diagonal change
      ind <- which(forplotCOX[,1] < coord_new[1] & forplotCOX[,1] > coord_new2[1])
      all[ind,] <- cbind(matrix(forplotCOX[ind,], ncol=2),
                        matrix(rep(forplot[k,], length(ind)), ncol=2, byrow=TRUE)
                        )
    }
  }

  points(forplot[,1], forplot[,2])
  lines(forplotCOX[,1], forplotCOX[,2], col="blue")
  abline(c(0,1), col = "red", lty=2)

  #correlations and SSR
  HRcheck <- cor(forplot[,2], forplot[,1]^exp(coxfit$coefficients))
  SSR <- sum((all[,4] - all[,2])^2)

  # ind <- max(which(all[,1] < 0.85 & all[,1] > 0.75)[1],
  #           which(all[,3] < 0.85 & all[,3] > 0.75)[1])
  #
  # text(round(all[ind, 1], 1), round(all[ind, 2], 1),
  #      "Cox Estimate", col="blue", pos=1, cex=0.8)
#
#   text(x=0.99, y=0.15,
#        labels = paste("SSR = ", round(SSR, 4), sep=""),
#        pos=2, cex = 0.8)
  text(x=0.99, y=0.1,
       labels = paste("rho = ", round(HRcheck, 4), sep=""),
       pos=2)

  return(list(KMres = KMres, SSR = SSR, rho = HRcheck))

}
