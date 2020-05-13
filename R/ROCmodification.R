#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric vector with times to event.
#' @param event Vector of indicators with value 1 subject experienced the event ocurred and 0 when censored.
#' @param group Vector of indicators for the treatment group with value 1 for treatment and 0 for control.
#' @param modifier Vector of indicators for the modifier with values of 1 when modifier is present and 0 otherwise.
#' @param mlabels A vector with two strings to label the levels (0, 1) of the modifier.
#' @param area TRUE/FALSE argument to indicate if user wants an estimate for area under the curve.
#' @param xlab passed from ROCsurv
#' @param ylab passed from ROCsurv
#' @param main passed from ROCsurv
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
#'
#' @export

ROCmodifier <- function(time, event, group, modifier, area=FALSE,
                        mlabels, xlab, ylab, main){

  all_lengths = c(length(time), length(event), length(group), length(modifier))
  if (length(unique(all_lengths)) != 1) stop("One or more input vectors (time, event, group, modifier) differs in length from the rest.")

  dat <- data.frame(time, event, group, modifier)
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
       xlab=xlab, ylab=ylab, main=main, cex.axis = 1.5, cex.lab = 1.5)

  for (k0 in 2:nrow(forplot0)) {
    coord_new = unname(forplot0[k0-1,])
    coord_new2 = unname(forplot0[k0,])
    segments(x0=coord_new[1], y0=coord_new[2],
             x1=coord_new2[1], y1=coord_new2[2], col="blue", lty = 1)
  }

  for (k1 in 2:nrow(forplot1)) {
    coord_new = unname(forplot1[k1-1,])
    coord_new2 = unname(forplot1[k1,])
    segments(x0=coord_new[1], y0=coord_new[2],
             x1=coord_new2[1], y1=coord_new2[2], col="black", lty = 2)
  }

  abline(c(0,1), col = "red", lty=3)

  legend("topleft", mlabels, col = c("blue", "black"), lty = c(1,2),
         inset=0.02, bg = "white", bty='n', seg.len = 0.7,
         x.intersp=0.9, y.intersp = 0.85, cex=1.5)
  return(list(coxfit = coxfit,
              modifier0 = list(control_km = KMres0$km_placebo,
                               drug_km = KMres0$km_drug),
              modifier1 = list(control_km = KMres1$km_placebo,
                               drug_km = KMres1$km_drug)
              )
        )

}
