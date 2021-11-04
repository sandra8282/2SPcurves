#' Create a two-sample cumulative incidence curve for data from randomized controlled trials.
#'
#' @description
#' This function creates a two-sample cumulative incidence curve for data from randomized controlled trials and has an option to diagnose the proportional hazard assumption.
#'
#' @param time Time to an event of any kind or to censoring.
#' @param event An indicator vector with values of 0, 1, 2, ..., M indicating one of M possible causes for an event for an participant that had an event occur or 0 if the participant was censored.
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
#' @param checkFG Logical argument to indicate if user wants to compare the nonparametric curve with the curve based on the Fine-Gray model (default is FALSE).
#' @param c_index Logical argument to indicate if user wants an estimate for the concordance for each risk (default is TRUE).
#' @param maxt Duration of the trial.
#' @param silent Logical argument, FALSE indicates the user wants plots and TRUE indicates no plots only calculations (default is FALSE).
#'
#' @return A plot of the ROC curve (if \code{silent=FALSE}) and an ROCsurv object containing:
#' \itemize{
#'  \item A Cuminc class object containing survival and cummulative incidence estimates and corresponding standard errors for each group.
#'  \item List with the two-sample cumulative incidence curves \code{C_u} for each risk type.
#' }
#'
#' @import survival
#' @importFrom dplyr distinct
#'
#' @export

TwoSCI <- function(time, event, group, maxt, xlab=NULL, ylab=NULL, main=NULL, rlabels,
                      cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5,
                      legend.inset=0.02, legend.cex=1.5, checkFG=FALSE, c_index = TRUE, silent){

  id <- 1:length(time)
  mymat <- data.frame(id, time, event, group)
  nrisktypes = length(unique(event)) - 1
  enames = c("censored", rlabels)
  eventf <- factor(event, 0:nrisktypes, labels = enames)
  mymatfg <- data.frame(id, time, eventf, group)
  #mymatfg <- mymatfg[mymatfg$time<maxt,]

  #get cum incidence

  fit <-  survfit(Surv(mymatfg$time, mymatfg$eventf) ~ mymatfg$group)
  sfit <- summary(fit)

  #get C(u)
  skm <- data.frame(time = sfit$time, group = as.numeric(sfit$strata)-1, ci = sfit$pstate[,2:3])
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

  #Check Fine-Gray model
  if (checkFG==TRUE){
    reg_res = predictions = list(); rname = coeffs= c();
        for (i in 1:max(event)){
          temp <- finegray(Surv(time, eventf) ~ ., data=mymatfg, etype=enames[i+1])
          fgfit <- coxph(Surv(temp$fgstart, temp$fgstop, temp$fgstatus) ~ temp$group,
                         weights=temp$fgwt)
          fgsurv <- survfit(fgfit, data.frame(group=c(0,1)))
          #fgmat <- cbind(fgsurv$time, fgsurv$cumhaz)
          reg_res[[i]] <- fgfit
          rname[i] <- enames[i+1]
          coeffs[i] <- fgfit$coefficients
          predictions[[i]] <- data.frame(x = fgsurv$cumhaz[,1], y = fgsurv$cumhaz[,2])
        }
    names(reg_res) <- rname
    list4fitplot <- list()
        for (skmi in (1:nrisktypes)){
          temp <- list_4plot[[skmi]]
          u <- temp[,1] #seq(0, max(temp[,1]), length.out = 50)
          list4fitplot[[skmi]] <- cbind(u, 1-((1-u)^exp(coeffs[[skmi]])))
        }
  }

  if (silent == FALSE) {
    #Plot
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    for (ploti in (1:nrisktypes)){
      lines(list_4plot[[ploti]][,1:2], type="s", lty = ploti+1, lwd = 2)
      if (checkFG==TRUE){lines(list4fitplot[[ploti]], lty = 1, lwd=2, type = "l")}
    }

    abline(c(0,1), col = "grey", lwd = lwd - 0.25)
    if (missing(rlabels)) {rlabels = as.character(1:nrisktypes)}
    legend("topleft", rlabels, lty = 2:(ploti+1), cex = legend.cex, inset= legend.inset,
           bg = "white", bty='n', seg.len = 1,
           x.intersp=0.9, y.intersp = 0.85, lwd = lwd - 0.25)

  }

  if (c_index==TRUE){
    if (missing(maxt)) {stop("Must provide maxt for calculation of c_index")}
    c <- comprsk_c(mymat, rlabels, maxt)
  }

  names(list_4plot) <- rlabels
  if (checkFG==TRUE){
    if (c_index==TRUE){
      return(list(cuminc = sfit, c_u = list_4plot, fits = reg_res, c_u_fit = list4fitplot, c_index = c))
      } else {return(list(cuminc = sfit, c_u = list_4plot, fits = reg_res, c_u_fit = list4fitplot))}
  } else {
    if (c_index==TRUE){
      return(list(cuminc = sfit, c_u = list_4plot, c_index = c))
    } else {return(list(cuminc = sfit, c_u = list_4plot))}}

}
