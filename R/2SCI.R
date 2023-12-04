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
#' @param maxt Duration of the trial.
#' @param silent Logical argument, FALSE indicates the user wants plots and TRUE indicates no plots only calculations (default is FALSE).
#' @param B The number of bootstrap samples to use for the confidence interval. Default = 1000.
#' @param lty Optional graphical parameter to set the type of line to use. Can be a number or a vector. See \link[graphics]{par} for more details.
#' @param CIabcd logical argument to indicate if user wants confidence interval for area between the curve and the diagonal (default is TRUE)
#' @param CI Optional logical argument to indicate if the user wants a bootstrap confidence interval for the curves and the concordance.
#' @param level A numerical argument to indicate the confidence level for the confidence interval. Default = 0.95.


#' @return A plot of the curve (if \code{silent=FALSE}) and an object containing:
#' \itemize{
#'  \item A Cuminc object containing survival and cummulative incidence estimates and corresponding standard errors for each group.
#'  \item List with the two-sample cumulative incidence curves \code{C_u} for each risk type.
#' }
#'
#' @import survival
#' @importFrom dplyr distinct
#' @importFrom zoo rollmean
#' @importFrom Bolstad2 sintegral
#@importFrom stats ks.test
#'
#' @export

TwoSCI <- function(time, event, group, maxt=NULL, xlab=NULL, ylab=NULL, main=NULL,
                   rlabels, cex.axis = 1.5, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = c(2, 3, 1),
                      legend.inset=0.02, legend.cex=1.5, checkFG=FALSE, CIabcd=TRUE, CI=FALSE,
                      level=0.95, B){

  if (missing(B)){B=length(time)}
  if (checkFG==TRUE & CI==TRUE) {stop("Cannot support both checkFG=TRUE and CI=TRUE")}
  maxt <- ifelse(is.null(maxt), max(time), maxt)
  mymat <- data.frame(id=1:length(time), time=time, event=event, group=group)

  if (0 %in% event){
    #censoring
    nrisktypes = length(unique(event)) - 1
    mymatfg <- data.frame(id=1:length(time), time=time, group=group)
    enames = c("censored", rlabels)
    mymatfg$eventf <- factor(event, 0:nrisktypes, labels = enames)

  } else {
    #no censoring
    mymatfg <- data.frame(id = 1:(length(time)+2), time = c(time, max(time)*10, max(time)*10), group = c(group, 0, 1))
    nrisktypes = length(unique(event))
    enames = c("censored", rlabels)
    newevent = c(event, 0,0)
    mymatfg$eventf <- factor(newevent, 0:nrisktypes, labels = enames)
  }

  fit <-  survfit(Surv(mymatfg$time, mymatfg$eventf) ~ mymatfg$group)
  sfit <- summary(fit)
  #get C(u)
  skm <- data.frame(time = sfit$time,
                    group = as.numeric(sfit$strata)-1,
                    ci = sfit$pstate[,2:3])

  #get cum incidence
  skm <- skm[order(skm[,1], skm[,2]),]
  ref_skm = list_4plot  = list(); abcds =  NULL #tests =
  for (skmi in (1:nrisktypes)){
    skmires = skm[, c(1, 2+skmi, 2)]
    skmires = skmires[!duplicated(skmires[,2:3]),]
    skmires = skmires[order(skmires$time, skmires$group),]
    ties_check <- duplicated(skmires[,1])
    if (length(ties_check) > 1) {
      ties_times = skmires[duplicated(skmires[,1]),1]
      ties_ind <- rep(0, nrow(skmires))
      ties_ind[which(skmires[,1] %in% ties_times)]=1
    } else {ties_ind <- rep(0, nrow(skmires))}
    skmires = cbind(skmires, ties_ind)
    ref_skm[[skmi]] = skmires
    temp = get4plotCumInc(skmires)
    id <- 1:nrow(temp)
    if (skmi==1){
      defaultW <- getOption("warn")
      options(warn = -1)
      abcdi <- sintegral(temp[,1],temp[,2]-temp[,1])$int
      list_4plot[[skmi]] = temp
      abcds[skmi] = abcdi
      #kstest = ks.test(x=temp[,2], y=temp[,1])
      #tests[skmi] = c(pval = round(kstest$p.value, 6))
      options(warn = defaultW)
    } else {
      defaultW <- getOption("warn")
      options(warn = -1)
      abcdi <- sintegral(temp[,1],temp[,1]-temp[,2])$int
      list_4plot[[skmi]] = temp
      abcds[skmi] = abcdi
      #kstest = ks.test(x=temp[,2], y=temp[,1])
      #tests[skmi] = c(pval = round(kstest$p.value, 6))
      options(warn = defaultW)
    }

  }

  abcds = c(abcds, sum(abcds))
    start = list_4plot[[1]][nrow(list_4plot[[1]]),1:2]
    end = 1-list_4plot[[2]][nrow(list_4plot[[2]]),1:2]
    earea = ((start[2]-start[1])+(end[2]-end[1]))*(end[1]-start[1])/2
  eabcd = abcds[3]+earea

  # ABCDres <- btsp2SCI(res = list_4plot,
  #                     maindat = mymat, nrisktypes=nrisktypes, B=B, level=level,
  #                     xlab, ylab, rlabels, main, CIabcd=TRUE, silent = TRUE)

  #Check Fine-Gray model
  #not for uncensored case
  if (checkFG==TRUE){
    reg_res = predictions = list(); rname = coeffs= c();
        for (i in 1:max(event)){
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

  ################################## plotting
  if (silent == FALSE & CI == FALSE) {
    if (checkFG==TRUE){
      plot2SCI(list_4plot, xlab=xlab, ylab=ylab, main=main, rlabels=rlabels,
             cex.axis = cex.axis, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = lty,
             legend.inset=legend.inset, legend.cex=legend.cex, list4fitplot)

    } else {
      plot2SCI(list_4plot, xlab=xlab, ylab=ylab, main=main, rlabels=rlabels,
                cex.axis = cex.axis, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = lty,
                legend.inset=legend.inset, legend.cex=legend.cex)
    }


  }

    if (checkFG==TRUE){
    results <- list(cuminc = sfit, c_u = list_4plot,
                      area.btw.curve.and.diag = abcds,
                      extrap.area.btw.curve.and.diag = eabcd,
                      fits = reg_res, c_u_fit = list4fitplot)
    } else {
        results <- list(cuminc = sfit, c_u = list_4plot,
                        area.btw.curve.and.diag = abcds,
                        extrap.area.btw.curve.and.diag = eabcd)
        }

  return(results)
}
