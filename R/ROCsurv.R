#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param time Time to event or censoring.
#' @param event An indicator vector with values of 1 for individuals who had the event occurred or 0 if the participant was censored.
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param method Method to be used to obtain a complete curve.
#' @param checkPH Logical argument to indicate if user wants to compare an ROC from the Cox Proportional Hazards Model to the KMROC (default is FALSE).
#' @param compare Logical argument to indicate if user wants want to compare the KMROC to ROM curves from a Cox model, a loglogistic AFT model and a lognormal AFT model (default is FALSE).
#' @param area Logical argument to indicate if user wants an estimate for area under the curve (default is TRUE).
#' @param abtwc Logical argument to indicate if area between curves is needed as part of model comparisons (checkPHM or compare must be TRUE).
#' @param silent Logical argument, FALSE indicates the user wants ROC plots and TRUE indicates no plots only AUC calculations (default is FALSE).
#' @param xlab String argument for the horizontal axis label of the ROC curve.
#' @param ylab String argument for the vertical axis label of the ROC curve.
#' @param main String argument for the title of the ROC curve.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling the inset of the legend.
#' @param legend.cex Optional graphical parameter for magnification of the legend's text.
#' @param lty Optional graphical parameter to set the type of line to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.
#' @param km Optional logical parameter. If TRUE the result will include survfit object containing survival information for each group.
#'
#' @return A plot of the ROC curve (if \code{silent=FALSE}) and an ROCsurv object containing:
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

###### # @param level The confidence level for the confidence interval of the area under the curve.
###### # Must be between 0.50 and 0.99. Default is 0.95. See details.

ROCsurv <- function(time, event, group, method,
                    checkPH = FALSE, compare=FALSE, area=TRUE, silent=FALSE, abtwc=FALSE,
                    xlab=NULL, ylab=NULL, main=NULL, KM, cex.axis = 1.5, cex.lab = 1.5,
                    legend.inset=0.02, legend.cex=1.5, lty = c(2,1,3), lwd = 1.5){

  #level=NULL,
  #### basic checks for missing parameters
  all_lengths = c(length(time), length(event), length(group))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group) differs in length from the rest.")
  if ((missing(method) + is.null(checkPH))==2) {checkPH <- TRUE}
  #if (is.null(level)) {level = 0.95}
  if (is.null(xlab)) {xlab <- "Control Group Survival"}
  if (is.null(ylab)) {ylab <- "Treatment Group Survival"}
  if (is.null(main)) {main <- ""}
  if (missing(method)) {method <- ""}
  if (missing(KM)) {KM <- FALSE}
  label.inset=legend.inset; label.cex=legend.cex

  ###########################################

  if(checkPH == TRUE) { #CHECK IF PROPORTIONAL HAZARDS
    if (missing(abtwc)) {abtwc=TRUE}
    result <- ROCandPHM(time, event, group, silent, abtwc, xlab, ylab, main, cex.axis = cex.axis,
                        cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                        label.cex = label.cex, lwd = lwd)
    #############################################
  } else {

    if (compare == TRUE) { #COMPARE TO LOGLOGISTIC AND LOGNORMAL
      if (missing(abtwc)) {abtwc=FALSE}
      result <- ROCcompare(time, event, group, silent, abtwc, xlab, ylab, main, cex.axis = cex.axis,
                           cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                           label.cex = label.cex, lwd = lwd)

    #############################################
    } else {
        KMests <- getKMtab(time, event, group)

        if (area==FALSE) {#plot AUC
            result <- onlyROC(KMests[[1]], xlab, ylab, main, cex.axis = cex.axis,
                              cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                              label.cex = label.cex, lwd = lwd)
            #############################################
        } else {
          #Return area for uncensored data
          if (KMests[[2]]==0) {
              result <- completeROC(KMests[[1]], silent, xlab, ylab, main,cex.axis = cex.axis,
                                    cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                                    label.cex = label.cex, lwd = lwd)} #time, event, group

          #Return area for censored data
          if (KMests[[2]]!=0 & method=="") {
            result <- completeROC(skm=KMests[[1]], silent, xlab, ylab, main, cex.axis,
                                  cex.lab, lty, label.inset, label.cex, lwd)
            } else {
              if (KMests[[2]]!=0 & method=="restrict") {
              result <- restrictROC(skm=KMests[[1]], silent, xlab, ylab, main, cex.axis = cex.axis,
                                    cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                                    label.cex = label.cex, lwd = lwd)}}

        }

        if (is.null(result)) {temp <- list(control_km=KMests[[3]], drug_km = KMests[[4]])
          } else {temp <- list(result, control_km=KMests[[3]], drug_km = KMests[[4]])}

        if(KM==TRUE){result <- temp}
    }
  }
  return(result)
}

