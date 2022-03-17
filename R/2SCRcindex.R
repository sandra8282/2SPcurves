

#' Calculate c-index for data from randomized controlled trials with competing risk.
#'
#' @description
#' This function calcultes c-index in precense of competing risks - corresponds to 2SCI curve.
#'
#' @param mymat passed from 2SCI.
#' @param rlabels passed from 2SCI.
#' @param maxt passed from 2SCI.
#'
#' @return A vector or single numeric value with concordance estimates.
#'
#' @import data.table
#' @importFrom dplyr %>%
#'
#' @useDynLib TwoSPC
#' @keywords internal
#' @noRd

comprsk_c <- function(mymat, rlabels, maxt){
  cens <- ifelse(mymat$event==0, 1, 0)
  GKM <- summary(survfit(Surv(mymat$time, cens) ~ 1, type='kaplan-meier'))
  mymat <- mymat[mymat$time<=maxt,]
  a=data.table(mymat); a[,merge:=a$time]
  b=data.table(time = GKM$time, w = GKM$surv); b[,merge:=b$time]
  setkeyv(a,c('merge')); setkeyv(b,c('merge'))
  mymat=b[a,roll=+1]
  indsna <- which(is.na(mymat$w))
  mymat$time[indsna] = 0;
  mymat$w[indsna] = 1;
  mymat <- mymat[,-c(1,3)]
  #n <- nrow(mymat)
  #`%!in%` <- Negate(`%in%`)


  r2f.cindex <- function(mymat){
    out1 <- .Fortran("cindex",
                     n=as.integer(nrow(mymat)),
                     nevent=as.integer(max(mymat$event)),
                     time=as.double(mymat$i.time),
                     event=as.integer(mymat$event),
                     group=as.integer(mymat$group),
                     w=as.double(mymat$w),
                     cind=as.double(0),
                     cwind=as.double(0),
                     NAOK = FALSE, PACKAGE = "TwoSPC")

    c <- out1$cind
    c2 <- out1$cwind
    names(c) = names(c2) = rlabels[1]

    return(c(c=c))
  }

  # c=rep(-10, max(mymat$event));
  # for (e in 1:max(mymat$event)){
  #   dnomcount = numcount = 0
  #   for (i in 1:n){
  #     for (j in 1:n){
  #       if (i==j){dnomcount=dnomcount; numcount=numcount
  #       } else {
  #         temp <- ifelse(mymat$event[i]==e & mymat$event[j]==e & mymat$group[i]==1 & mymat$group[j]==0, 1, 0)
  #         dnomcount = dnomcount + temp
  #         temp2 <- ifelse(mymat$i.time[j] <= mymat$i.time[i] & mymat$event[i]==e & mymat$event[j]==e &
  #                           mymat$group[i]==1 & mymat$group[j]==0, 1,0)
  #         numcount = numcount + temp2
  #       }
  #     }
  #   }
  #   c[e] <- numcount/dnomcount
  # }
  # names(c) <- rlabels
  #
  # c2=c();
  # for (e in 1:max(mymat$event)){
  #   dnomcount = numcount = 0.000000000000000001
  #   for (i in 1:n){
  #     for (j in 1:n){
  #       if (i!=j){
  #         temp <- ifelse(mymat$event[i]==e & mymat$event[j]==e & mymat$group[i]==1 & mymat$group[j]==0, 1, 0)
  #         dnomcount = dnomcount + temp/(mymat$w[i]*mymat$w[j])
  #         temp2 <- ifelse(mymat$i.time[j] <= mymat$i.time[i] & mymat$event[i]==e & mymat$event[j]==e &
  #                           mymat$group[i]==1 & mymat$group[j]==0, 1,0)
  #         numcount = numcount + temp2/(mymat$w[i]*mymat$w[j])
  #       }
  #     }
  #   }
  #   c2[e] <- (numcount-0.000000000000000001)/(dnomcount-0.000000000000000001)
  # }
  # names(c2) <- rlabels

  cres<- r2f.cindex(mymat)

  return(cres)
}

