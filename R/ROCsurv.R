#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param method Method for calculating AUC if neither of the survival curves reaches 0. See Details.
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
#'
#' @examples
#'
#' # Simulate a clinical trial with 500 subjects that were followed for 2 years where the true hazard
#' # ratio of 0.2 indicates a large treatment effect, assumming a data comes from a Wiebull(0.5, 1.5).
#' n = 500
#' maxt = 2
#' beta = log(0.2)
#'
#' set.seed(28)
#' simdata <- simSurvTrial(size = n, followup = maxt, beta = beta, dist = "weibull")
#'
#' #ROC curve
#' sim_roc <- with(simdata, ROCsurv(time, event, treatment))
#'
#' @export
#'
ROCsurv <- function(time, event, group, method){

  if (missing(method)) {method <- "restrict"}
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

  ties_check <- unique(table(skm[,1]))

  if (length(ties_check) > 1) {
    ties_times = as.integer(names(which(table(skm[,1])>1)))
    ties_ind <- rep(NA, nrow(skm))
    ties_ind[which(skm[,1] %in% ties_times)]=1
    ties_ind[-which(skm[,1] %in% ties_times)]=0
  } else {ties_ind <- rep(0, nrow(skm))}

  skm = cbind(skm, ties_ind)

  mskm <- min(skm[,2])

  if (mskm==0) {result <- completeROC(skm)}

  if (mskm!=0 & method=="restrict") {
      result <- restrictROC(skm)
  }

  # if (mskm!=0 & method=="ph_loglog") {
  #   cox <- coxph(Surv(time, event) ~ group, ties="breslow")
  #   result <- ph_loglogROC(skm, cox)
  # }

  return(list(control_km = km_placebo,
              drug_km = km_drug,
              AUC = result,
              var = NA))

}
