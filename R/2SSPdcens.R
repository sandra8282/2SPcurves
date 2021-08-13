#' Create an ROC curve for survival data from clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param o_l Time at the left of the interval (L,R] for the onset time (see details).
#' @param o_r Time at the right of the interval (L,R] for the onset time (see details).
#' @param e_l Time at the left of the interval (L,R] for the event time (see details).
#' @param e_r Time at the right of the interval (L,R] for the event time (see details).
#' @param maxtime Maximum event time in the sample.
#' @param ctime Last time an event could have been observed in the study (length of the study, must be higher than maxtime even if by 1 unit).
#' @param dtype String argument that can take values "dcens" or "DIcens" to indicate either doubly censored or doubly interval censored data (see details).
#' @param group Indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param xlab String argument for the horizontal axis label of the curve.
#' @param ylab String argument for the vertical axis label of the curve.
#' @param main String argument for the title of the curve.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling legend inset.
#' @param legend.cex Optional graphical parameter for magnification of a legend's text.
#' @param lty Optional graphical parameter to set the type of lines to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item R_u: Data frame with the estimated two-sample survival curves.
#'  \item npmle0: NPMLE's (outcome and gap) for the untreated/unexposed group.
#'  \item npmle1: NPMLE's (outcome and gap) for the treated/exposed group.
#' }
#'
## @importFrom dblcens d011
#' @importFrom DescTools AUC
#'
#' @details
#' For doubly censored data as discussed by Turnbull (dtype = "dcens"):
#'   Set o_l = NA and o_r = NA
#'   Set e_l = e_r if the event is exact
#'   Set e_l = -Inf if the event time is left censored and e_r = Inf if it is right censored
#'
#' For doubly interval censored data (dtype = "DIcens) as discussed by De Gruttola and Legakos:
#'   Set o_l = o_r if the onset time is exact, o_l = -Inf if it is left censored, and o_r = Inf if it is right censored.
#'   Set e_l = e_r if the event time is exact, e_l = -Inf if it is left censored, and e_r = Inf if it is right censored.
#'
#' @export
#'

TwoSSPdcens <- function(o_l, o_r,  e_l, e_r, group, maxtime = NULL, ctime = 99, dtype=NULL,
                  xlab=NULL, ylab=NULL, main=NULL, cex.axis = 1.5,
                  cex.lab = 1.5, legend.inset=0.02, legend.cex=1.5, lty = c(1,2,3), lwd = 1.5) {

  if (is.null(xlab)) {xlab <- "Control Group Survival"}
  if (is.null(ylab)) {ylab <- "Treatment Group Survival"}
  if (is.null(main)) {main <- ""}
  if (is.null(dtype)) {stop("Error: must indicate type of double censoring using dtype")}
  if (dtype=="dcens"){
    check1 <- is.na(o_l)&is.na(o_r)
    if (check1!=TRUE) {stop("Error: either o_l and o_r must be NA when dtype='dcens'.")}
    all_len <- c(length(e_l), length(e_r), length(group))
    if (length(unique(all_len))!=1) {stop("Error: One or more variables (e_l, e_r, group) defer in length.")}
    dat <- data.table(e_l, e_r, group)
    control <- data.frame(dat[group==0]); trt <- data.frame(dat[group==1])

    ######################################################################
    ##################### ADD CODE FOR ESTIMATION HERE ###################
    ######################################################################

  }

  if (dtype=="DIcens"){
    all_len <- c(length(o_l), length(o_r), length(e_l), length(e_r), length(group))
    if (length(unique(all_len))!=1) {stop("Error: One or more variables (o_l, o_r, e_l, e_r, group) defer in length.")}
    dat <- data.table(o_l, o_r, e_l, e_r, group)
    dat0 <- data.frame(dat[group==0,])
    dat1 <- data.frame(dat[group==1,])
    DISKM <- getDKMtab(dat0, dat1, dtype=dtype, maxtime, ctime)
    forplot <- get4plot(DISKM$dskm)
    forplot <- forplot[-duplicated(forplot),]
    plot(forplot, type = "l", lty = lty[1], lwd = lwd, #las=1
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    abline(c(0,1), col = "darkgrey", lty=1, lwd = lwd-0.25)
  }

  auc = AUC(forplot[,1], forplot[,2])
  returnobj <- list(R_u = forplot, auc = auc, npmle0 = DISKM$npmle0, npmle1 = DISKM$npmle1)

  return(returnobj)

}
