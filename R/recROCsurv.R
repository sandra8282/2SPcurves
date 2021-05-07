#' Create an ROC curve for recurrent event data that meets exchangeability assumption for a clinical trials.
#'
#' @description
#' This function creates ROC curve for survival data from clinical trials and calculates the area under the curve.
#'
#' @param id Unique identifier for each participant.
#' @param episode_j Vector identifying the event number (j=1,..., m_i) for each partipant.
#' @param time Gap time (Y_ij)
#' @param mi Vector with number of episodes/recurrences each participant experienced.
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param xlab String argument for the horizontal axis label of the ROC curve.
#' @param ylab String argument for the vertical axis label of the ROC curve.
#' @param main String argument for the title of the ROC curve.
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
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @details
#' The methods avaiable are "restrict" or "complete"
#'
#' @export
#'
ROCsurvRec <- function(id, episode_j, time, group, mi,
                       xlab=NULL, ylab=NULL, main=NULL, cex.axis = 1.5, cex.lab = 1.5,
                       legend.inset=0.02, legend.cex=1.5, lty = c(2,1,3), lwd = 1.5){

    dat <- data.frame(id, episode_j, time, group, mi)
    dat$delta <- ifelse(mi>1, 1, 0)
    dat <- dat[order(id,episode_j),]

    MWCR <- getKMtabRecurrent(dat)

    res<-onlyROC(MWCR[[1]], xlab, ylab, main, cex.axis = cex.axis,
            cex.lab = cex.lab, lty = lty, lwd = lwd)

}
