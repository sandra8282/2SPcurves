#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param level The confidence level for the confidence interval of the area under the curve.
#' Must be between 0.50 and 0.99. Default is 0.95. See details.
#' @param method choose
#' @param checkPH to check
#'
#' @return A plot with the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object containing the survival curve for the treatment group.
#'  \item A survfit object containing the survival curve for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @details
#' The methods avaiable are "restrict" or "complete"
#'
#' @examples
#' # Simulate a clinical trial with 300 subjects with 1 year followup and 10% drop-out rate where:
#' # the true hazard ratio of 2 indicates a large treatment effect.
#' # We assumming the data comes from a Wiebull(0.5, 1.5).
#' n = 500
#' maxt = 2
#' true_AUC = 0.8
#' beta = log((1 - true_AUC)/true_AUC)
#'
#' set.seed(28)
#' simdata <- simSurvTrial(size = n, followup = maxt, beta = beta, dist = "weibull",
#'                          params = c(1, 1.5))
#'
#' #ROC curve
#' sim_roc <- with(simdata, ROCsurv(time = time, event = event, group = treatment))
#'
#' @export
#'
ROCsurv <- function(time, event, group, level=NULL, method = NULL, checkPH = NULL){
  #time = dat$ti ; event = dat$di; group = dat$trt; level = 0.95;
  #time = leukemia$Time ; event = leukemia$Event; group = leukemia$Group; level = 0.95;

  all_lengths = c(length(time), length(event), length(group))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group) differs in length from the rest.")
  if ((is.null(method) + is.null(checkPH))==2) {checkPH <- TRUE}
  if (is.null(level)) {level = 0.95}

  if(checkPH == TRUE) {
    result <- ROCandPHM(time, event, group)
    return(result)
  } else {
    KMests <- getKMtab(time, event, group)
  #Point Estimate based on
    if (KMests[[2]]==0) {
      result <- completeROC(KMests[[1]], silent=FALSE)
      return(list(control_km = KMests$km_placebo,
                drug_km = KMests$km_drug,
                AUC = result))
             }
    if (KMests[[2]]!=0 & is.null(method)) {result <- onlyROC(KMests[[1]])
      return(list(control_km = KMests$km_placebo,
                  drug_km = KMests$km_drug))

    } else if (KMests[[2]]!=0 & method=="restrict") {
      result <- restrictROC(KMests[[1]], silent = FALSE)
      return(list(control_km = KMests$km_placebo,
                  drug_km = KMests$km_drug,
                  AUC = result))
    }
  }

  #Calculate bootstrap variance
  #SEandCI <- btsp(time, event, group, method, B = 1000, level)

  # if (mskm!=0 & method=="ph_loglog") {
  #   cox <- coxph(Surv(time, event) ~ group, ties="breslow")
  #   result <- ph_loglogROC(skm, cox)
  # }


              #SEandCI[1],
              #SEandCI[2:3]))

}
