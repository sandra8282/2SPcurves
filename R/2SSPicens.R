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
#' @importFrom stats complete.cases
#' @importFrom DescTools AUC
#' @importFrom stats approx
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

  ### survival estimations for missing portions through interpolation
    control_surv <- control_pf %>% mutate(surv = 1- cumdrop, surv0 = lag(surv, default = 1))
    control_surv$R <- ifelse(control_surv$R==Inf, control_surv$R[nrow(control_surv)-1]+control_surv$R[nrow(control_surv)-1]*2, control_surv$R)
    control_surv_new <- data.frame(t=0, s=1)
    for (i in 1:nrow(control_surv)) {
      ti = approx(control_surv[i,1:2], c(control_surv$surv0[i],control_surv$surv[i]), n = 3)
      sfct = data.frame(t = ti$x, s = ti$y)
      control_surv_new = rbind(control_surv_new, sfct)
    }

    trt_surv <- trt_pf %>% mutate(surv = 1- cumdrop, surv0 = lag(surv, default = 1))
    trt_surv$R <- ifelse(trt_surv$R==Inf, trt_surv$R[nrow(trt_surv)-1]+trt_surv$R[nrow(trt_surv)-1]*2, trt_surv$R)
    trt_surv_new <- data.frame(t=0, s=1)
    for (i in 1:nrow(trt_surv)) {
      ti = approx(trt_surv[i,1:2], c(trt_surv$surv0[i],trt_surv$surv[i]), n = 3)
      sfct = data.frame(t = ti$x, s = ti$y)
      trt_surv_new = rbind(trt_surv_new, sfct)
    }


  # add missing portions to two sample curve ######################################################################################
    control_surv_new$group = 0;
    trt_surv_new$group = 1
    surv_new = distinct(rbind(control_surv_new, trt_surv_new))
    surv_new = surv_new[order(surv_new$t, -surv_new$s),]
    KMests <- data.frame(getIntSKM(surv_new))
    res_temp <- data.frame(get4plot(skm = KMests))
    if (res_temp$x[nrow(res_temp)]==0|res_temp$y[nrow(res_temp)]==0){
      res_temp[nrow(res_temp)+1,]=c(0, 0, 1)
    }
    res_temp = distinct(res_temp)
    final_res <-res_temp[,1:2]
    colnames(final_res) = colnames(res)[1:2]
    warn = getOption("warn")
    options(warn=-1)
    auc = AUC(final_res[,1], final_res[,2], method = c("trapezoid"))
    options(warn=warn)

    # 2 sample curve PH and PO models  #########################################################################################
  if (checkPH == TRUE | checkPO == TRUE) {
    res2 <- PHPOicens(checkPH=checkPH, checkPO = checkPO, dat=dat, res=res, lty=lty,
                      legend.inset=legend.inset, legend.cex=legend.cex, lwd=lwd)
      if (is.null(res2$fit_ph)){
        returnobj <- list(two_sample_prob = res, auc = auc,
                          NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                          fit_po =res2$fit_po)
      } else {
        if (is.null(res2$fit_po)){
          returnobj <- list(two_sample_prob = res, auc = auc,
                            NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                            fit_ph =res2$fit_ph)
        } else {returnobj <- list(two_sample_prob = res, auc = auc,
                                  NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt,
                                  fit_ph =res2$fit_ph, fit_po =res2$fit_po)}
      }
  } else {
    lines(res_temp, lty = 2)
    returnobj <- list(two_sample_prob = res, auc = auc,
                      NPMLE_control = NPMLE.control, NPMLE_trt = NPMLE.trt)
  }

  return(returnobj)

}
