#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param left Time at the left of the interval (L,R]
#' @param right Time at the right of the interval (L,R]
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param end_follow A numeric value indicating the end of study follow-up. Default is the largest value for \code{right}.
#' @param iterations Number of maximum iterations for the EM algorithm. Default is 1000.
#' @param xlab String argument for the horizontal axis label of the curve.
#' @param ylab String argument for the vertical axis label of the curve.
#' @param main String argument for the title of the curve.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling legend inset.
#' @param legend.cex Optional graphical parameter for magnification of a legend's text.
#' @param lty Optional graphical parameter to set the type of lines to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.
#' @param checkPH Logical argument to indicate if user wants to compare the nonparametric curve with the curve based on the proportional hazards model (default is FALSE).
#' @param checkPO Logical argument to indicate if user wants to compare the nonparametric curve with the curve based on the proportional odds model (default is FALSE).
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item Data frame with the estimated survival curves.
#'  \item icsurv object from application of EM to the control group data.
#'  \item icsurv object from application of EM to the treated group data.
#' }
#'
#' @importFrom Icens EMICM
#' @importFrom data.table data.table
#' @importFrom dplyr %>%
#' @importFrom dplyr full_join
#' @importFrom dplyr lag
#' @importFrom stats complete.cases
#' @importFrom data.table data.table
#'
#'
#' @export
#'

TwoSSPicens <- function(left, right, group, iterations, end_follow, checkPH = FALSE, checkPO = FALSE,
                  xlab=NULL, ylab=NULL, main=NULL, cex.axis = 1.5,
                  cex.lab = 1.5, legend.inset=0.02, legend.cex=1.5, lty = c(1,2,3), lwd = 1.5) {

  if (is.null(xlab)) {xlab <- "Control Group Survival"}
  if (is.null(ylab)) {ylab <- "Treatment Group Survival"}
  if (is.null(main)) {main <- ""}
  if(missing(iterations)) {iterations = 5000}

  all_len <- c(length(left), length(right), length(group))
  if (length(unique(all_len))!=1) {stop("Error: One or more variables (left, right, group) defer in length.")}

  dat <- data.table(left, right, group)
  control <- data.frame(dat[group==0]); trt <- data.frame(dat[group==1])

  # EM calculations for control group #################################################
    NPMLE.control<-EMICM(control[,(1:2)], maxiter = iterations)
    control_pf <- data.frame(L = NPMLE.control$intmap[1,],
                             R = NPMLE.control$intmap[2,],
                             drop = NPMLE.control$pf)
    control_pf <- control_pf[control_pf$drop!=0,]
    control_pf$cumdrop = cumsum(control_pf$drop)

  # EM calculations for treatment group #################################################
    NPMLE.trt<-EMICM(trt[,(1:2)], maxiter = iterations)
    trt_pf <- data.frame(L = NPMLE.trt$intmap[1,],
                         R = NPMLE.trt$intmap[2,],
                         drop = NPMLE.trt$pf)
    trt_pf <- trt_pf[trt_pf$drop!=0,]
    trt_pf$cumdrop = cumsum(trt_pf$drop)

  # 2 sample curve   ##############################################################################################

  res <- getIGroc(npmle_0 = control_pf, npmle_1 = trt_pf,
                  xlab, ylab, main, cex.axis, cex.lab, lwd)

    # 2 sample curve PH and PO models  #########################################################################################
  if (checkPH == TRUE | checkPO == TRUE) {
    res2 <- PHPOicens(checkPH=checkPH, checkPO = checkPO, dat=dat, res=res,
                      legend.inset=legend.inset, legend.cex=legend.cex, lwd=lwd)
      if (is.null(res2$fit_ph)){
        returnobj <- list(two_sample_prob = res[[1]], auc = res[[2]],
                          NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                          fit_po =res2$fit_po)
      } else {
        if (is.null(res2$fit_po)){
          returnobj <- list(two_sample_prob = res[[1]], auc = res[[2]],
                            NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                            fit_ph =res2$fit_ph)
        } else {returnobj <- list(two_sample_prob = res[[1]], auc = res[[2]],
                                  NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                                  fit_ph =res2$fit_ph, fit_po =res2$fit_po)}
      }
  } else {
    returnobj <- list(two_sample_prob = res[[1]], auc = res[[2]],
                      NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt)
  }

  return(returnobj)

}
