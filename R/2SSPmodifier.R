#' ROC when survival goes to 0 for either group
#'
#' @param time Time to event or censoring.
#' @param event An indicator vector with values of 1 for individuals who had the event occurred or 0 if the participant was censored.
#' @param group An indicator vector with values of 1 if the the participant was in the treatment arm and 0 otherwise.
#' @param modifier An indicator vector with values of 1 for a certain level of the modifier and 0 otherwise. Example: 0 for male and 1 for female.
#' @param mlabels A vector with two strings to label the levels of the modifier.
#' @param area Logical argument to indicate if user wants an estimate for area under the curve. Default is TRUE.
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
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A list with two survfit objects for the treatment and control groups when the modifier=0.
#'  \item A list with two survfit objects for the treatment and control groups when the modifier=1.
#'  \item If area = TRUE, the areas under the curve for the two ROCs in the given plot.
#' }
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @importFrom stats na.omit
#' @import survival
#' @importFrom stats cor
#' @importFrom data.table data.table
#'
#' @export

TwoSSPmodifier <- function(time, event, group, modifier, area=FALSE,
                        mlabels, xlab, ylab, main, cex.axis = 1.5, cex.lab = 1.5,
                        legend.inset=0.02, legend.cex=1.5, lty = c(1,6), lwd = 1.5){

  all_lengths = c(length(time), length(event), length(group), length(modifier))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group, modifier) differs in length from the rest.")

  if (missing(xlab)) {xlab <- "Control Group Survival"}
  if (missing(ylab)) {ylab <- "Treatment Group Survival"}
  if (missing(main)) {main <- ""}


  dat <- data.table(time, event, group, modifier)
  dat0 <- dat[dat$modifier == 0, ]
  dat1 <- dat[dat$modifier == 1, ]

  KMres0 <- getKMtab(dat0$time, dat0$event, dat0$group)
  skm0 <- KMres0[[1]]
  forplot0 = get4plot(skm0)

  KMres1 <- getKMtab(dat1$time, dat1$event, dat1$group)
  skm1 <- KMres1[[1]]
  forplot1 = get4plot(skm1)

  coxfit <- coxph(Surv(time, event) ~ group + modifier + group:modifier)

  plot(NULL, type="n", las=1,
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)

  abline(c(0,1), col = "grey", lty=1, lwd = lwd - 0.25)

  for (k0 in 2:nrow(forplot0)) {
    coord_new = unname(forplot0[k0-1,])
    coord_new2 = unname(forplot0[k0,])
    segments(x0=coord_new[1], y0=coord_new[2],
             x1=coord_new2[1], y1=coord_new2[2], col="black", lty = lty[1], lwd = lwd)
  }

  for (k1 in 2:nrow(forplot1)) {
    coord_new = unname(forplot1[k1-1,])
    coord_new2 = unname(forplot1[k1,])
    segments(x0=coord_new[1], y0=coord_new[2],
             x1=coord_new2[1], y1=coord_new2[2], col="black", lty = lty[2], lwd = lwd)
  }


  legend("topleft", mlabels, lty = lty,
         inset=legend.inset, bg = "white", bty='n', seg.len = 1,
         x.intersp=0.9, y.intersp = 0.85, cex=legend.cex, lwd=lwd)

  return(list(coxfit = coxfit,
              modifier0 = list(control_km = KMres0$km_placebo,
                               drug_km = KMres0$km_drug),
              modifier1 = list(control_km = KMres1$km_placebo,
                               drug_km = KMres1$km_drug)
              )
        )

}
