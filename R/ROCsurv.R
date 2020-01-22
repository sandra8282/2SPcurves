#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param method Method for calculating AUC if neither of the survival curves reaches 0. See Details.
#' @param level The confidence level for the confidence interval of the area under the curve.
#' Must be between 0.50 and 0.99. Default is 0.95. See details.
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
#' @details
#' The methods avaiable are "restrict" or "complete"
#'
#' @examples
#' # Simulate a clinical trial with 500 subjects that were followed for 2 years where the true hazard
#' # ratio of 0.2 indicates a large treatment effect, assumming a data comes from a Wiebull(0.5, 1.5).
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
ROCsurv <- function(time, event, group, method, level){

  all_lengths = c(length(time), length(event), length(group))

  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group) differs in length from the rest.")

  if (missing(method)) {method <- "restrict"}
  if (missing(level)) {level = 0.95}

  KMests <- getKMtab(time, event, group)

  #Point Estimate based on
  if (KMests[[2]]==0) {result <- completeROC(KMests[[1]], silent=FALSE)}
  if (KMests[[2]]!=0 & method=="restrict") {result <- restrictROC(KMests[[1]], silent = FALSE)}

  tictoc::tic()
  #Calculate bootstrap variance
  SEandCI <- btsp(time, event, group, method, B = 1000, level)
  tictoc::toc()

  # if (mskm!=0 & method=="ph_loglog") {
  #   cox <- coxph(Surv(time, event) ~ group, ties="breslow")
  #   result <- ph_loglogROC(skm, cox)
  # }

  return(list(control_km = km_placebo,
              drug_km = km_drug,
              AUC = result,
              SEandCI[1],
              SEandCI[2])
         )

}
