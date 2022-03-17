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
#' @param checkFG Logical argument to indicate if the user wants to compare the nonparametric curve with the curve based on the Fine-Gray model (default is FALSE).
#' @param c_index Logical argument to indicate if the user wants an estimate for the concordance for each risk (default is TRUE).
#' @param maxt Duration of the trial.
#' @param CI Optional logical argument to indicate if the user wants a bootstrap confidence interval for the curves and the concordance.
#' @param silent Logical argument, FALSE indicates the user wants plots and TRUE indicates no plots only calculations (default is FALSE).
#' @param level A numerical argument to indicate the confidence level for the confidence interval. Default = 0.95.
#' @param B The number of bootstrap samples to use for the confidence interval. Default = 1000.
#' @param lty Optional graphical parameter to set the type of line to use. Can be a number or a vector. See \link[graphics]{par} for more details.

#' @return A plot of the curve (if \code{silent=FALSE}) and an object containing:
#' \itemize{
#'  \item A Cuminc object containing survival and cummulative incidence estimates and corresponding standard errors for each group.
#'  \item List with the two-sample cumulative incidence curves \code{C_u} for each risk type.
#' }
#'
#' @import survival
#' @importFrom dplyr distinct
#' @importFrom Bolstad2 sintegral
#'
#' @export

TwoSCI <- function(time, event, group, maxt=NULL, xlab=NULL, ylab=NULL, main=NULL, rlabels,
                      cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = c(2, 3, 1),
                      legend.inset=0.02, legend.cex=1.5, checkFG=FALSE, c_index = TRUE, CI=FALSE,
                      level=0.95, B){

  if (missing(B)){B=length(time)}

  if (checkFG==TRUE & CI==TRUE) {stop("Cannot support both checkFG=TRUE and CI=TRUE")}

  id <- 1:length(time)
  mymat <- data.frame(id, time, event, group)

  # if(0 %in% mymat$event){
  #
  # }
  nrisktypes = length(unique(event)) - 1
  mymatfg <- data.frame(id, time, group)

  if (0 %in% event){
    enames = c("censored", rlabels)
    mymatfg$eventf <- factor(event, 0:nrisktypes, labels = enames)
  } else {
    eventf <- factor(event, 1:length(unique(event)), labels = rlabels)
    mymatfg$eventf <- factor(eventf, levels = c(levels(eventf), "censored"))
  }

  #mymatfg <- mymatfg[mymatfg$time<maxt,]

  #get cum incidence

  fit <-  survfit(Surv(mymatfg$time, mymatfg$eventf) ~ mymatfg$group)
  sfit <- summary(fit)

  #get C(u)
  skm <- data.frame(time = sfit$time,
                    group = as.numeric(sfit$strata)-1,
                    ci = sfit$pstate[,2:3])
  skm <- skm[order(skm[,1], skm[,2]),]
  list_4plot = list(); abcds = NULL
  for (skmi in (1:nrisktypes)){
    skmires = matrix(as.numeric(as.matrix(distinct(skm[, c(1, 2+skmi, 2)]))),
                ncol=3)
    ties_check <- duplicated(skmires[,1])
    if (length(ties_check) > 1) {
      ties_times = skmires[duplicated(skmires[,1]),1]
      ties_ind <- rep(0, nrow(skmires))
      ties_ind[which(skmires[,1] %in% ties_times)]=1
    } else {ties_ind <- rep(0, nrow(skmires))}
    skmires = cbind(skmires, ties_ind)
    temp = get4plotCumInc(skmires)
    temp = temp[!duplicated(temp),]
    temp = temp[-(duplicated(temp[,1:2])&temp[,2]==0),]
    id <- 1:nrow(temp)
    z <- abs(temp[,2] - temp[,1])
    defaultW <- getOption("warn")
    options(warn = -1)
    abcdi <- sintegral(temp[,1],z)$int
    options(warn = defaultW)
    list_4plot[[skmi]] = temp
    abcds[i] = abcdi
  }

  #Check Fine-Gray model
  if (checkFG==TRUE){
    reg_res = predictions = list(); rname = coeffs= c();
        for (i in 1:max(event)){
          temp <- finegray(Surv(time, eventf) ~ ., data=mymatfg, etype=enames[i+1])
          fgfit <- coxph(Surv(temp$fgstart, temp$fgstop, temp$fgstatus) ~ temp$group,
                         weights=temp$fgwt, control = coxph.control(timefix = FALSE))
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


  if (silent == FALSE & CI == FALSE) {
    ltyploti <- c(5, 3, 6, 1)
    #Plot
    plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)
    for (ploti in (1:nrisktypes)){
      lines(list_4plot[[ploti]][,1:2], lty = ltyploti[ploti], lwd = 2)
      if (checkFG==TRUE){lines(list4fitplot[[ploti]], lty = 1, lwd=2, type = "l")}
    }

    abline(c(0,1), col = "grey", lwd = lwd - 0.25)
    if (missing(rlabels)) {rlabels = as.character(1:nrisktypes)}
    legend("topleft", rlabels, lty = ltyploti[1:length(rlabels)],
           cex = legend.cex, inset= legend.inset,
           bg = "white", bty='n', seg.len = 1,
           x.intersp=0.9, y.intersp = 0.85, lwd = lwd - 0.25)

  }

  if (silent==TRUE & CI==TRUE) {warning("CI=TRUE and silent=TRUE. Will ignore silent=TRUE command and plot curve with confidence interval.")}

  if (c_index==TRUE){
    if (missing(maxt)) {stop("Must provide maxt for calculation of c_index")}
    c <- comprsk_c(mymat, rlabels, maxt)
    if (CI==TRUE){
      BSTPres <- btsp2SCI(res = list_4plot, maindat = mymat, maxt = maxt,
                      nrisktypes, B=B, level, xlab, ylab, rlabels, main,
                      cex.axis = cex.axis, cex.lab = cex.lab, lty = lty,
                      lwd = lwd, bst_c = TRUE, cindex = c)
    }
  } else {
    if (CI==TRUE){
      BSTPres <- btsp2SCI(res = list_4plot, maindat = mymat, maxt = max(mymat[,2]),
                      nrisktypes, B=B, level, xlab, ylab, rlabels, main,
                      cex.axis = cex.axis, cex.lab = cex.lab, lty = lty,
                      lwd = lwd, bst_c = FALSE, cindex = c)
    }
  }

  names(list_4plot) = names(abcds) =rlabels

  if (checkFG==TRUE){
    if (c_index==TRUE){
      return(list(cuminc = sfit, c_u = list_4plot, area.btw.curve.and.diag = abcds,
                  fits = reg_res, c_u_fit = list4fitplot,
                  c_index = c))
      } else {return(list(cuminc = sfit, c_u = list_4plot, fits = reg_res, c_u_fit = list4fitplot))}
  } else {
    if (c_index==TRUE){
      if (CI==TRUE){
        return(list(cuminc = sfit, c_u = BSTPres$C_u, area.btw.curve.and.diag = BSTPres$abcd,
                    c_index = BSTPres$Cindex))
      } else {return(list(cuminc = sfit, c_u = list_4plot,
                          area.btw.curve.and.diag = abcds, c_index = c))}
    } else {return(list(cuminc = sfit, c_u = list_4plot, area.btw.curve.and.diag = abcds))}}

}
