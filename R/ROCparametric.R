getSKM4fit <- function(time, fitsurv, group) {
  fitskm <- cbind(time, fitsurv, group)
  fitskm <- fitskm[!duplicated(fitskm),]
  ties_check <- unique(table(fitskm[,1])) #how many times are repeated
  if (length(ties_check) > 1) {
    ties_times = fitskm[duplicated(fitskm[,1]),1]
    ties_ind <- rep(0, nrow(fitskm))
    ties_ind[which(fitskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(fitskm))}
  fitskm = cbind(fitskm, ties_ind)
  fitskm <- fitskm[order(fitskm[,1], fitskm[,3]),]
  return(fitskm)
}

#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param dist String indicating the distribution.
#' @param xlab passed from ROCsurv
#' @param ylab passed from ROCsurv
#' @param silent passed from ROCsurv
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
#' @import survival
#' @importFrom stats vcov
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom utils capture.output
#'
#' @export

ROCparametric <- function(time, event, group, dist="lognormal",
                          silent, xlab, ylab, main) {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  #obtain coefficients and parameters to calculate survival fcts
  fit = survreg(Surv(time, event) ~ group, dist=dist)
  mu = fit$coefficients[1]; gamma = fit$coefficients[2]; sigma = fit$scale;
  lambda = exp(-unname(mu/sigma)); beta = -unname(gamma)/sigma; alpha = 1/sigma;
  covmat <- vcov(fit)

  #calculate survival fcts
  if (dist=="loglogistic") {
    fitsurv <- (1 + lambda*exp(beta*group)*(time^alpha))^(-1)
    fitskm <- getSKM4fit(time, fitsurv, group)
    forplotfit = get4plot(fitskm)
    newsurv0 <- seq(min(forplotfit[,1]), 0.001, length.out = round(nrow(forplotfit)*0.1, 0))
    forplotfit <- rbind(forplotfit, cbind(newsurv0, rep(0, length(newsurv0))))
    forplotfit[,2] = 1 / (1 + (exp(- gamma/sigma) * (1/forplotfit[,1] - 1)))
  }

  if (dist == "lognormal") {
    fitsurv <- 1 - pnorm((log(time)- mu - gamma*group)/sigma)
    fitskm <- getSKM4fit(time, fitsurv, group)
    forplotfit = get4plot(fitskm)
    newsurv0 <- seq(min(forplotfit[,1]), 0.001, length.out = round(nrow(forplotfit)*0.1, 0))
    forplotfit <- rbind(forplotfit, cbind(newsurv0, rep(0, length(newsurv0))))
    forplotfit[,2] = 1 - pnorm(qnorm(1 - forplotfit[,1]) - gamma/sigma)
  }

   area = 0

   if (silent==TRUE) {
     for (k in 2:nrow(forplotfit)) {
       coord_new = unname(forplotfit[k-1,])
       coord_new2 = unname(forplotfit[k,])
       #figure out areas
       if (forplotfit[k,2]==forplotfit[k-1,2]) {#move horizontally
         area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
       } else {
         if (forplotfit[k,1]!=forplotfit[k-1,1] & forplotfit[k,2]!=forplotfit[k-1,2]){
           #area for diagonal
           area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
           area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
           area = area + area_rectang + area_triang
         }
       }
     }
   } else {
     plot(NULL, type="n", las=1,
          xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
          xlab=xlab, ylab=ylab, main=main, cex.axis = 1.5, cex.lab = 1.5)

     for (k in 2:nrow(forplotfit)) {
       coord_new = unname(forplotfit[k-1,])
       coord_new2 = unname(forplotfit[k,])
       #figure out areas and shading
       if (forplotfit[k,2]==forplotfit[k-1,2]) {#move horizontally
         rect(xright = coord_new[1], ytop = coord_new[2],
              xleft = coord_new2[1], ybottom = 0,
              col = "pink", border = "pink")
         area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
       } else {
         if (forplotfit[k,1]!=forplotfit[k-1,1] & forplotfit[k,2]!=forplotfit[k-1,2]){
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
     }

     points(forplot[,1], forplot[,2], col = "grey50", cex = 0.75)
     lines(forplotfit[,1], forplotfit[,2], lty=1)
     abline(c(0,1), col = "red", lty=2)
     text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
          pos=2, cex = 1)
   }

   colnames(forplotfit) = c("u", "R(u)")
   return(list(fit, area, parametricROC = forplotfit))

}
