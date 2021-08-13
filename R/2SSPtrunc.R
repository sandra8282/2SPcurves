#' Create an ROC curve for left truncated right censored data.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param entry The entry time for the participant that determined left truncation.
#' @param exit Time when the event or censoring occurred (must be higher than entry time).
#' @param group An indicator vector with values of 1 if the the participant was in the treatment or exposed group and 0 otherwise.
#' @param event An indicator vector with values of 1 for individuals who had the event occur or 0 if the participant was censored.
#' @param area Logical argument to indicate if user wants an estimate for area under the curve (default is TRUE).
#' @param silent Logical argument, FALSE indicates the user wants ROC plots and TRUE indicates no plots only AUC calculations (default is FALSE).
#' @param CI Logical argument indicating if user wants bootstrap errors and confidence intervals to be calculated for the curve and AUC. Default is FALSE.
#' @param level Numerical argument indicating condifence level for the confidence interval. Default is 0.95. It's ignored if CI=FALSE.
#' @param B Number of bootstrap samples to use for confidence interval and error calculations. Default is 1000. It's ignored if CI=FALSE.
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
#'  \item Data frame with the estimated survival curves.
#'  \item icsurv object from application of EM to the control group data.
#'  \item icsurv object from application of EM to the treated group data.
#' }
#'
#' @export
#'

TwoSSPtrunc <- function(entry, exit, event, group, area=FALSE, silent=FALSE,
                   CI = FALSE, level = 0.95, B=1000,
                   xlab=NULL, ylab=NULL, main=NULL, cex.axis = 1.5, cex.lab = 1.5,
                   legend.inset=0.02, legend.cex=1.5, lty = c(2,1,3), lwd = 1.5){

  #### basic checks for missing parameters
  all_lengths = c(length(entry), length(exit), length(event), length(group))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (entry, exit, event, group) differs in length from the rest.")
  if (is.null(xlab)) {xlab <- "Control Group Survival Probability"}
  if (is.null(ylab)) {ylab <- "Treatment Group Survival Probability"}
  if (is.null(main)) {main <- ""}
  label.inset=legend.inset; label.cex=legend.cex

  mat <- cbind(entry, exit, event, group)
  mat <- na.omit(mat); out <- which(mat[,1]>=mat[,2]);

  if (length(out)>=1) {
    mat <- mat[-out, ]
    print(
    paste("Observations", paste(out, collapse = ', '), "were deleted because entry time >= exit time", sep = " "))}
  entry = mat[,1]; exit = mat[,2]; event = mat[,3]; group = mat[,4];

  ###########################################

    MKMests <- getMKMtab(entry, exit, event, group, c=1, a=0.25)

    unexposed_mkm <- data.frame(rbind(c(0, NA, NA, 1), MKMests[[5]]))
    exposed_mkm <- data.frame(rbind(c(0, NA, NA, 1), MKMests[[6]]))
    if (area==FALSE & CI == FALSE) {#plot
      temp <- onlyROC(skm = MKMests[[1]], xlab, ylab, main, cex.axis,
                      cex.lab, lty, label.inset, label.cex, lwd)
      res <- list(unexposed_km=MKMests[[3]], exposed_km = MKMests[[4]],
                  unexposed_mkm=unexposed_mkm, exposed_mkm = exposed_mkm)
    }

    #############################################
    if (area == TRUE & CI == FALSE) {
      res <- getAUC(MKMests, silent, xlab, ylab, main,cex.axis = cex.axis,
                    cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                    label.cex = label.cex, lwd = lwd)
      res <- list(AUC = res$AUC, R_u = res$R_u,
                  unexposed_km=MKMests[[3]], exposed_km = MKMests[[4]],
                  unexposed_mkm=unexposed_mkm, exposed_mkm = exposed_mkm)
    }
    #############################################
    if (CI == TRUE){
      res <- getAUC(MKMests, silent=TRUE, xlab, ylab, main,cex.axis = cex.axis,
                    cex.lab = cex.lab, lty = lty, label.inset = label.inset,
                    label.cex = label.cex, lwd = lwd)
      maindat <- mat
      temp <- btsp(res, maindat, B, level, xlab, ylab, main, cex.axis = cex.axis,
                   cex.lab = cex.lab, lty = lty, lwd = lwd)
      res <- list(temp, unexposed_km=MKMests[[3]], exposed_km = MKMests[[4]],
                  unexposed_mkm=unexposed_mkm, exposed_mkm = exposed_mkm)
    }

  return(res)
}
