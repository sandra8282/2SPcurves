#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param time Time to event or censoring.
#' @param event Vector indicating if the event occurred (event=1) or if the subject was censored (event=0).
#' @param group Vector indicating if the individual was in the treatment arm (group=1) or the control arm (group=0).
#' @param level The confidence level for the confidence interval of the area under the curve.
#' Must be between 0.50 and 0.99. Default is 0.95. See details.
#' @param method Method to be used to obtain a complete curve.
#' @param checkPH logical argument to indicate if user wants to compare an ROC from the Cox Proportional Hazards Model to the KMROC (default is FALSE).
#' @param compare logical argument to indicate if user wants want to compare the KMROC to ROM curves from a Cox model, a loglogistic AFT model and a lognormal AFT model (default is FALSE).
#' @param area logical argument to indicate if user wants an estimate for area under the curve (default is TRUE).
#' @param abtwc logical argument to indicate if area between curves is needed as part of model comparisons (checkPHM or compare must be TRUE).
#' @param silent logical argument, FALSE indicates the user wants ROC plots and TRUE indicates no plots only AUC calculations (default is FALSE).
#' @param xlabel string argument to use as label for the horizontal axis of the ROC curve.
#' @param ylabel string argument to use as label for the vertical axis of the ROC curve.
#' @param main string argument to use as the title of the ROC curve.
#'
#' @return An optional plot with the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object containing the survival curve for the treatment group.
#'  \item A survfit object containing the survival curve for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @details
#' The methods avaiable are "restrict" or "complete"
#'
#' @export
#'
ROCsurv <- function(time, event, group, level=NULL, method = NULL,
                    checkPH = FALSE, compare=FALSE, area=NULL, silent, abtwc,
                    xlabel, ylabel, main){
  #time = dat$ti ; event = dat$di; group = dat$trt; level = 0.95;
  #time = leukemia$Time ; event = leukemia$Event; group = leukemia$Group; level = 0.95;

  all_lengths = c(length(time), length(event), length(group))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group) differs in length from the rest.")
  if ((is.null(method) + is.null(checkPH))==2) {checkPH <- TRUE}
  if (is.null(level)) {level = 0.95}
  if (missing(silent)) {silent=FALSE}
  if (missing(xlabel)) {xlab <- "Control Group Survival"} else {
    xlab <- paste(xlabel, "Survival", sep = " ")}
  if (missing(ylabel)) {ylab <- "Treatment Group Survival"} else {
    ylab <- paste(ylabel, "Survival", sep = " ")}
  if (missing(main)) {main <- ""}

  if(checkPH == TRUE) { #CHECK IF PROPORTIONAL HAZARDS
    if (missing(abtwc)) {abtwc=TRUE}
    result <- ROCandPHM(time, event, group, silent, abtwc, xlab, ylab, main)
    return(result)

  } else if (compare == TRUE) { #COMPARE TO LOGLOGISTIC AND LOGNORMAL
    if (missing(abtwc)) {abtwc=TRUE}
    result <- ROCcompare(time, event, group, silent, abtwc, xlab, ylab, main)
    return(result)

  } else {
    KMests <- getKMtab(time, event, group)

    if (area==FALSE) {
        result <- onlyROC(KMests[[1]], xlab, ylab, main)
        return(list(control_km = KMests$km_placebo,
                drug_km = KMests$km_drug))
    } else {

      #Point Estimate based on
      if (KMests[[2]]==0) {
          result <- completeROC(KMests[[1]], silent, xlab, ylab, main) #time, event, group
          return(list(control_km = KMests$km_placebo,
                      drug_km = KMests$km_drug,
                      AUC = result))}

      if (KMests[[2]]!=0 & is.null(method)) {result <- onlyROC(KMests[[1]])
          return(list(control_km = KMests$km_placebo,
                      drug_km = KMests$km_drug))
      } else if (KMests[[2]]!=0 & method=="restrict") {
          result <- restrictROC(KMests[[1]], silent)

    }

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
