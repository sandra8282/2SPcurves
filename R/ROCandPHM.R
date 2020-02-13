#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
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
  forplot = get4plot(skm)

  coxfit <- coxph(Surv(time, event) ~ group, ties = "breslow")

  sum_sqrres = 0

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
