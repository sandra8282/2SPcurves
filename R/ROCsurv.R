#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#'
#' @return A plot with the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object containing the survival curve for the treatment group.
#'  \item A survfit object containing the survival curve for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @import survival
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @export

ROCsurv <- function(time, event, group){

  km_placebo <- survfit(Surv(time, event) ~ 1,
                        subset=(group==0), type='kaplan-meier')
  km_drug <- survfit(Surv(time, event) ~ 1,
                     subset=(group==1), type='kaplan-meier')

  skm_p <- cbind(time = summary(km_placebo)$time,
                 surv = summary(km_placebo)$surv)
  skm_p <- cbind(skm_p, type = rep(0, nrow(skm_p)))
  skm_d <- cbind(time = summary(km_drug)$time,
                 surv = summary(km_drug)$surv)
  skm_d <- cbind(skm_d, type = rep(1, nrow(skm_d)))

  skm <- rbind(skm_p, skm_d)
  skm <- skm[order(skm[,1]),]

  plot(c(0,1), c(0, 1), type="n", xlab="", ylab="")
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
        cex.main = 1)

  area = 0
  for (i in 1:nrow(skm)) {
    if(i<2){
      coord_new = c(1, 1)
      #check if drug or placebo
      if (skm[i,3]==0) {#move horizontally
        coord_new2 = c(skm[i,2], 1)
        rect(xright = coord_new[1], ytop = coord_new[2],
             xleft = coord_new2[1], ybottom = 0,
             col = "pink", border = "pink")
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {#move vertically
        coord_new2 = c(skm[i,2], 1)
      }
    } else {
      #check if drug or placebo
      if (skm[i,3]==0) {#move horizontally
        coord_new2 = c(skm[i,2], coord_new[2])
        rect(xright = coord_new[1], ytop = coord_new[2],
             xleft = coord_new2[1], ybottom = 0,
             col = "pink", border = "pink")
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {#move vertically
        coord_new2 = c(coord_new[1], skm[i,2])
      }
      segments(coord_new[1], coord_new[2],
               coord_new2[1], coord_new2[2], col="black")
      coord_new = coord_new2
    }
  }
  abline(c(0,1), col = "black", lty=2)
  area = unname(area)
  text(0.8, 0.15, paste("AUC=", round(area,2), sep=""))

  return(list(control_KaplanMeier = km_placebo,
              treatment_KaplanMeier = km_drug,
              AUC = area))
}
