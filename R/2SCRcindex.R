#`%!in%` <- Negate(`%in%`)

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
  n <- nrow(mymat)
  n1 <- nrow(mymat %>% dplyr::filter(mymat$group == 1))
  n0 <- nrow(mymat %>% dplyr::filter(mymat$group == 0))

  c=NULL; dnom2 = num = 0
  for (e in 1:max(mymat$event)){
    for (i in 1:n1){
      for (j in 1:n0){
        if (i==j){dnomcount=dnomcount; numcount=numcount
        } else {
          if (mymat$i.time[j] <= mymat$i.time[i] & mymat$group[i]==1 & mymat$group[j]==0 &
              mymat$event[i] %!in% c(e, 0) & mymat$event[j] %!in% c(e, 0)){
              denom2 <- denom2+1
          } else {
            temp <- ifelse(mymat$i.time[j] <= mymat$i.time[i] &
                            mymat$group[i]==1 & mymat$group[j]==0 &
                             mymat$event[i]==e & mymat$event[j]==e , 1, 0)
            num = num+temp
          }
        }
      }
    }
    c[e] <- num/(1-dnom2)
  }
  names(c) <- rlabels

  c2=NULL; dnom2 = num = 0
  for (e in 1:max(mymat$event)){
    for (i in 1:n1){
      for (j in 1:n0){
        if (i==j){dnomcount=dnomcount; numcount=numcount
        } else {
          if (mymat$i.time[j] <= mymat$i.time[i] & mymat$group[i]==1 & mymat$group[j]==0 &
              mymat$event[i] %!in% c(e, 0) & mymat$event[j] %!in% c(e, 0)){
            denom2 <- denom2+(1/(mymat$w[i]*mymat$w[j]))
          } else {
            temp <- ifelse(mymat$i.time[j] <= mymat$i.time[i] &
                             mymat$group[i]==1 & mymat$group[j]==0 &
                             mymat$event[i]==e & mymat$event[j]==e , 1, 0)
            num = num+(temp/(mymat$w[i]*mymat$w[j]))
          }
        }
      }
    }
    c2[e] <- num/(1-dnom2)
  }
  names(c2) <- rlabels

  # c=NULL; dnomcount = numcount = 0
  # for (e in 1:max(mymat$event)){
  #   for (i in 1:n1){
  #     for (j in 1:n0){
  #       if (i==j){dnomcount=dnomcount; numcount=numcount
  #         } else {
  #         temp <- ifelse(mymat$event[i]==e & mymat$event[j]==e & mymat$group[i]==1 & mymat$group[j]==0, 1, 0)
  #         dnomcount = dnomcount + temp
  #         temp2 <- ifelse(mymat$i.time[j] <= mymat$i.time[i] & mymat$event[i]==e & mymat$event[j]==e &
  #                               mymat$group[i]==1 & mymat$group[j]==0, 1,0)
  #         numcount = numcount + temp2
  #       }
  #     }
  #   }
  #   c[e] <- numcount/dnomcount
  # }
  # names(c) <- rlabels
  #
  # c2=c(); dnomcount = numcount = 0.000000000000000001
  # for (e in 1:max(mymat$event)){
  #   for (i in 1:n1){
  #     for (j in 1:n0){
  #       if (i!=j){
  #         temp <- ifelse(mymat$event[i]==e & mymat$event[j]==e & mymat$group[i]==1 & mymat$group[j]==0, 1, 0)
  #           dnomcount = dnomcount + temp/(mymat$w[i]*mymat$w[j])
  #           temp2 <- ifelse(mymat$i.time[j] <= mymat$i.time[i] & mymat$event[i]==e & mymat$event[j]==e &
  #                             mymat$group[i]==1 & mymat$group[j]==0, 1,0)
  #           numcount = numcount + temp2/(mymat$w[i]*mymat$w[j])
  #       }
  #     }
  #   }
  #   c2[e] <- (numcount-0.000000000000000001)/(dnomcount-0.000000000000000001)
  # }
  # names(c2) <- rlabels

  return(list(c=c, c_adj = c2))
}


