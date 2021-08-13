getAUCipw <- function(AKMests, silent, xlab, ylab, main, cex.axis = cex.axis,
                   cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                   label.cex = label.cex, lwd = lwd){
  result <- completeROC(AKMests[[1]], silent, xlab, ylab, main,cex.axis = cex.axis,
                        cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                        label.cex = label.cex, lwd = lwd) #time, event, group
  if (is.null(result)) {temp <- AKMests[[3]]
  } else {temp <- list(AUC = result[[1]], R_u = result[[2]], adjKM=AKMests[[3]])}

  return(temp)
}
#' Create an inverse probability of treatment weight adjusted two-sample survival probability curve for data from observational studies.
#'
#' @description
#' This function creates an inverse probability of treatment weight adjusted two-samle survival probability curve for data from observational studies and has options diagnose the proportional hazard assumption.
#'
#' @param time Time to event or censoring.
#' @param event An indicator vector with values of 1 for individuals who had the event occurred or 0 if the participant was censored.
#' @param group An indicator vector with values of 1 if the the participant was in the exposed group and 0 otherwise.
#' @param weights A vector with weights for each study participant.
#' @param area Logical argument to indicate if user wants an estimate for area under the curve (default is TRUE).
#' @param silent Logical argument, FALSE indicates the user wants ROC plots and TRUE indicates no plots only AUC calculations (default is FALSE).
#' @param checkPH Logical argument to indicate if user wants to compare the nonparametric curve with the curve based on the proportional hazards model (default is FALSE).
#' @param abtwc Logical argument to indicate if area between curves is needed as part of model comparisons (checkPHM or compare must be TRUE).
#' @param xlab String argument for the horizontal axis label of the ROC curve.
#' @param ylab String argument for the vertical axis label of the ROC curve.
#' @param main String argument for the title of the ROC curve.
#' @param CI Logical argument indicating if user wants bootstrap errors and confidence intervals to be calculated for the curve and AUC. Default is FALSE.
#' @param level Numerical argument indicating condifence level for the confidence interval. Default is 0.95. It's ignored if CI=FALSE.
#' @param B Number of bootstrap samples to use for confidence interval and error calculations. Default is 1000. It's ignored if CI=FALSE.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling the inset of the legend.
#' @param legend.cex Optional graphical parameter for magnification of the legend's text.
#' @param lty Optional graphical parameter to set the type of line to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.
#'
#' @return A plot of the ROC curve (if \code{silent=FALSE}) and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object containing the survival curve for the treatment group.
#'  \item A survfit object containing the survival curve for the control group.
#'  \item The area under the curve (AUC).
#'  \item A matrix representation of the two-sample survival probability curve \code{R_u}.
#' }
#'
#'
#' @export
#'

###### # @param level The confidence level for the confidence interval of the area under the curve.
###### # Must be between 0.50 and 0.99. Default is 0.95. See details.

TwoSSPiptw <- function(time, event, group, weights, area = FALSE, CI = FALSE, level = 0.95, B=1000,
                       checkPH = FALSE, abtwc=FALSE, silent = FALSE,
                       xlab=NULL, ylab=NULL, main=NULL, cex.axis = 1.5, cex.lab = 1.5,
                       legend.inset=0.02, legend.cex=1.5, lty = c(2,1,3), lwd = 1.5){

  #### basic checks for missing parameters
  all_lengths = c(length(time), length(event), length(group), length(weights))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group, weights) differs in length from the rest.")
  if (is.null(xlab)) {xlab <- "Not exposed Group Survival"}
  if (is.null(ylab)) {ylab <- "Exposed Group Survival"}
  if (is.null(main)) {main <- ""}
  #if (missing(method)) {method <- ""}
  label.inset=legend.inset; label.cex=legend.cex

  mat <- cbind(time, event, group, weights)
  mat <- na.omit(mat)
  time = mat[,1]; event = mat[,2]; group = mat[,3]; weights = mat[,4];

  ###########################################

  if(checkPH == TRUE) { #CHECK IF PROPORTIONAL HAZARDS
    if (missing(abtwc)) {abtwc=TRUE}
    res <- PHMiptw(time, event, group, weights, abtwc, xlab, ylab, main, cex.axis = cex.axis,
                        cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                        label.cex = label.cex, lwd = lwd)
    #############################################
  }

  if (checkPH == FALSE){

    AKMests <- getAKMtab(time, event, group, weights)

        if (area==FALSE & CI == FALSE) {#plot
          temp <- onlyROC(skm = AKMests[[1]], xlab, ylab, main, cex.axis,
                          cex.lab, lty, label.inset, label.cex, lwd)
          res <- list(adjKM = AKMests$adjkm)
          }

            #############################################
        if (area == TRUE & CI == FALSE) {
          res <- getAUCipw(AKMests, silent, xlab, ylab, main,cex.axis = cex.axis,
                             cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                             label.cex = label.cex, lwd = lwd)
          }
            #############################################
        if (CI == TRUE){
            res <- getAUCipw(AKMests, silent=TRUE, xlab, ylab, main,cex.axis = cex.axis,
                        cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                        label.cex = label.cex, lwd = lwd)
            maindat <- mat
            temp <- btspIPW(res, maindat, B, level, xlab, ylab, main, cex.axis = cex.axis,
                             cex.lab = cex.lab, lty = lty, lwd = lwd)
            res <- list(temp, adjKM=AKMests[[3]])
          }

      }


  return(res)
}
