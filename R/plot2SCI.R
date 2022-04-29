#' Plot two-sample cumulative incidence curves
#'
#' @description
#' This function plots a two-sample cumulative incidence curve for data from randomized controlled trials and has an option to diagnose the proportional hazard assumption.
#'
#' @param list_4plot from 2SCI
#' @param list4fitplot from 2SCI
#' @param xlab from 2SCI
#' @param ylab from 2SCI
#' @param main from 2SCI
#' @param rlabels from 2SCI
#' @param cex.axis from 2SCI
#' @param cex.lab from 2SCI
#' @param legend.inset from 2SCI
#' @param legend.cex from 2SCI
#' @param lwd from 2SCI
#' @param lty from 2SCI
#'
#' @importFrom dplyr distinct
#' @importFrom zoo rollmean
#' @importFrom stats ks.test
#' @importFrom Bolstad2 sintegral
#'
#' @noRd
#' @keywords internal
#'

plot2SCI <- function(list_4plot,list4fitplot=NULL, xlab=NULL, ylab=NULL, main=NULL, rlabels,
                     cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = c(2, 3, 5, 6),
                     legend.inset=0.02, legend.cex=1.5, areas, btspind=FALSE){

  if (btspind==FALSE){

  }
  myvals = c(max(list_4plot[[1]][,2]), max(list_4plot[[2]][,2]))
  mxvals = c(max(list_4plot[[1]][,1]), max(list_4plot[[2]][,1]))
  nrisktypes=length(list_4plot)
  par(mar = c(5, 4, 4, 4) + 0.35)
  plot(list_4plot[[1]][,1:2], type="l", lty = lty[1], lwd = 2, #type="n", las=1,
       xlim=c(0,sum(mxvals)), ylim = c(0, sum(myvals)), #to make tight axis: xaxs="i", yaxs="i"
       xlab=paste("Cause 1", xlab), ylab=paste("Cause 1", ylab),
       main=main, cex.axis = cex.axis, cex.lab = cex.lab)
  abline(c(0,1), col = "grey", lwd = lwd - 0.25)
  if (!is.null(list4fitplot)){lines(list4fitplot[[1]], lty = 1, lwd=2, type = "l")}

  # for (i in 2:nrow(list_4plot[[1]])){
  #   coord = list_4plot[[1]][i-1,1:2]
  #   coordnew = list_4plot[[1]][i,1:2]
  #   if (coordnew[2]>coordnew[1]){}
  #   polygon(x=c(coord[1], coord[1], coordnew[1], coordnew[1]),
  #           y=c(coord[2], coord[1], coordnew[1], coordnew[2]),
  #           col = "grey", border = NA)
  # }
  abline(v=mxvals[1])


  par(new = TRUE)
  plot(list_4plot[[2]][,1:2], xlim=rev(c(0,1)), ylim=rev(c(0,1)),
       type="l", lty = lty[2], lwd = 2, axes = FALSE, xlab = "", ylab = "")
  if (!is.null(list4fitplot)){lines(rev(list4fitplot[[2]][,1]), rev(list4fitplot[[2]][,2]),
                           lty = 1, lwd=2, type = "l")}
  axis(side = 4, at = rev(seq(0,1,0.2)), cex.axis = cex.axis)
  axis(side = 3, at = rev(seq(0,1,0.2)), cex.axis = cex.axis)
  mtext(paste("Comp. Risk(s)", ylab), side = 4, line = 3, cex = cex.lab)
  mtext(paste("Comp. Risk(s)", xlab), side = 3, line = 3, cex = cex.lab)

  abline(v=mxvals[2])
  if (!is.null(list4fitplot)){lines(rev(list4fitplot[[2]][,1]), rev(list4fitplot[[2]][,2]),
                           lty = 1, lwd=2, type = "l")}

  axis(side = 4, at = rev(seq(0,1,0.2)), cex.axis = cex.axis)
  axis(side = 3, at = rev(seq(0,1,0.2)), cex.axis = cex.axis)
  mtext(paste("Comp. Risk(s)", ylab), side = 4, line = 3, cex = cex.lab)
  mtext(paste("Comp. Risk(s)", xlab), side = 3, line = 3, cex = cex.lab)

  # for (i in 2:nrow(list_4plot[[2]])){
  #   coord = list_4plot[[2]][i-1,1:2]
  #   coordnew = list_4plot[[2]][i,1:2]
  #   if (coordnew[2]>coordnew[1]){}
  #     polygon(x=c(coord[1], coord[1], coordnew[1], coordnew[1]),
  #             y=c(coord[2], coord[1], coordnew[1], coordnew[2]),
  #             col = "grey", border = NA)
  # }

  # fct <- list_4plot[[1]][,1:2]
  # fct <- data.frame(fct, f1 = fct[,2]-fct[,1])
  # defaultW <- getOption("warn")
  # options(warn = -1)
  # area1 <- Bolstad2::sintegral(fct$x, fct$f1)$int
  # options(warn = defaultW)
  #
  # fct2 <- list_4plot[[2]][,1:2]
  # fct2 <- data.frame(fct2, f2 = fct2[,1]-fct2[,2])
  # defaultW <- getOption("warn")
  # options(warn = -1)
  # area2 <- Bolstad2::sintegral(fct2$x, fct2$f2)$int
  # options(warn = defaultW)

  # names(area1) = c(rlabels[1], paste("ref", rlabels[1], sep = "_"))
  # names(area2) = c(rlabels[2], paste("ref", rlabels[2], sep = "_"))
  #pAAC = list(area1, area2)

  if (missing(rlabels)) {rlabels = as.character(1:nrisktypes)}
    legend("topleft", rlabels, lty = lty[1:length(rlabels)],
           cex = legend.cex, inset= legend.inset,
           bg = "white", bty='n', seg.len = 1,
           x.intersp=0.9, y.intersp = 0.85, lwd = lwd - 0.25)

  return(NULL)
}
