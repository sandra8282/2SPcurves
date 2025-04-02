#' Check Proportional Hazards and/or Proportional Odds Based Curve(s)
#'
#' @param checkPH passed from ic2ss.
#' @param checkPO passed from ic2ss.
#' @param dat passed from ic2ss.
#' @param res passed from ic2ss.
#' @param legend.inset passed from ic2ss.
#' @param legend.cex passed from ic2ss.
#' @param lwd passed from ic2ss.
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item An icenReg object for the treatment group.
#'  \item An icenReg objectfor the control group.
#' }
#'
#' @import graphics
#' @importFrom icenReg ic_sp
#'
#' @keywords internal
#' @noRd


PHPOicens <- function(checkPH, checkPO, dat, res, legend.inset, legend.cex, lwd, lty){
  mu <- seq(min(res$u),1,0.001)
  fit_ph = fit_po = NULL

  if (length(lty)<3) {lty = c(1,3,5)}

  if (checkPH==TRUE){
    fit_ph <- ic_sp(cbind(left, right) ~ group, model = 'ph', bs_samples = 5, data = dat)
    lines(mu, mu^exp(coef(fit_ph)), lty = lty[2], lwd = lwd)
    legendtext <- c("Nonparametric", "Proportional Hazards")
    ltyl <- lty[c(1:2)]
  }

  if (checkPO==TRUE){
    fit_po <- ic_sp(cbind(left, right) ~ group, model = 'po', bs_samples = 5, data = dat)
    oddsu <- mu/(1-mu)
    ebeta <- exp(coef(fit_po))
    lines(mu, ebeta*oddsu/(1+ebeta*oddsu), lty = lty[3], lwd = lwd)
    legendtext <- c("Nonparametric", "Proportional Odds")
    ltyl = lty[c(1,3)]
  }

  if (checkPH==TRUE & checkPO==TRUE){
    ltyl = lty
    legendtext <- c("Nonparametric", "Proportional Hazards", "Proportional Odds")
  }

  legend("topleft", legendtext, lty = ltyl,
         inset=legend.inset, cex=legend.cex, bg = "white", bty='n', seg.len = 0.8,
         x.intersp=0.9, y.intersp = 0.85, lwd = lwd)

  return(list(fit_ph = fit_ph, fit_po = fit_po))
}
