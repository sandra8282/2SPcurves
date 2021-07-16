#' Create a two-sample cumulative incidence curve for data from randomized controlled trials.
#'
#' @description
#' This function creates a two-sample cumulative incidence curve for data from randomized controlled trials and has an option to diagnose the proportional hazard assumption.
#'
#' @param time Time to an event of any kind or to censoring.
#' @param event An indicator vector with values of 1, 2, ..., M related to one of M possible causes for an event for an participant that had an event occur or 0 if the participant was censored.
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param xlab String argument for the horizontal axis label of the ROC curve.
#' @param ylab String argument for the vertical axis label of the ROC curve.
#' @param main String argument for the title of the ROC curve.
#' @param rlabels A vector with M strings to label the curves related to each of the M causes for events.
#' @param cex.axis Optional graphical parameter for magnification of axis annotation. See \link[graphics]{par} for more details.
#' @param cex.lab Optional graphical parameter for magnification of x and y labels. See \link[graphics]{par} for more details.
#' @param legend.inset Optional graphical parameter controling the inset of the legend.
#' @param legend.cex Optional graphical parameter for magnification of the legend's text.
#' @param lwd Optional graphical parameter for line width relative to the default. See \link[graphics]{par} for more details.
#'
#' @return A plot of the ROC curve (if \code{silent=FALSE}) and an ROCsurv object containing:
#' \itemize{
#'  \item A Cuminc class object containing survival and cummulative incidence estimates and corresponding standard errors for each group.
#'  \item A matrix representation of the two-sample cumulative incidence curve \code{C_u}.
#' }
#'
#' @import cmprsk
#' @importFrom mstate Cuminc
#'
#'
#' @details
#' The methods avaiable are "restrict" or "complete"
#'
#' @export

TwoSCI <- function(time, event, group, xlab=NULL, ylab=NULL, main=NULL, rlabels,
                      cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5,
                      legend.inset=0.02, legend.cex=1.5){

  ci <- Cuminc(time= time, status = event, group=group)
  nrisktypes = length(unique(event)) - 1
  skm <- ci[, c(2, 1, 4:(4+nrisktypes-1))]
  skm <- skm[order(skm[,1], skm[,2]),]

  ties_check <- unique(table(skm[,1]))
  if (length(ties_check) > 1) {
    ties_times = skm[duplicated(skm[,1]),1]
    ties_ind <- rep(0, nrow(skm))
    ties_ind[which(skm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(skm))}

  skm = cbind(skm, ties_ind)

  list_4plot = list()
  for (skmi in (1:nrisktypes)){
    skmires = matrix(as.numeric(as.matrix(distinct(skm[, c(1, 2+skmi, 2, 5)]))),
                ncol=4)
    list_4plot[[skmi]] = get4plotCumInc(skmires)
  }

  plot(NULL, type="n", las=1,
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
  for (ploti in (1:nrisktypes)){
    lines(list_4plot[[ploti]][,1:2], type="s", lty = ploti, lwd = 2)
  }
  abline(c(0,1), col = "grey", lwd = lwd - 0.25)
  if (missing(rlabels)) {rlabels = as.character(1:nrisktypes)}
  legend("topleft", rlabels, lty = 1:ploti, cex = legend.cex, inset= legend.inset,
         bg = "white", bty='n', seg.len = 1,
         x.intersp=0.9, y.intersp = 0.85, lwd = lwd - 0.25)

  names(list_4plot) <- rlabels
  return(list(cuminc = ci, c_u = list_4plot))

}
