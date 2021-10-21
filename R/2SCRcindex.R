

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
  #mymat <- mymat[mymat$time<maxt,]

  a=data.table(mymat); a[,merge:=a$time]
  b=data.table(time = GKM$time, w = GKM$surv); b[,merge:=b$time]
  setkeyv(a,c('merge')); setkeyv(b,c('merge'))
  mymat=b[a,roll=+1]
  indsna <- which(is.na(mymat$w))
  mymat$w[indsna] = 1;
  mymat <- mymat[,-c(1,3)]
  n <- nrow(mymat)
  `%!in%` <- Negate(`%in%`)
  colnames(mymat)[3] = "time"
  inds <- which(mymat$time %in% GKM$time)
  mymat$w2 <- mymat$w
  for (ic in inds){
    mymat$w2[ic] <- mymat$w2[ic-1]
  }

  c_h = c_p = c_u = c_w1 = c_w2 = c_w3 = c_u2 = NULL;
  denom = denom2 = denom3 = denom4 = denom5 = denom6 = 0
  num = num2 = num3 = num4 = num5 = num6 = Q = cij = 0
  for (e in 1:max(mymat$event)){
    for (i in 1:n){
      for (j in 1:n){
        if (i!=j){
            if (mymat$event[i]==e & mymat$event[j]==e){
              Q = Q + 1
              cij = cij + ifelse(mymat$group[i]<mymat$group[j] & mymat$time[i]<mymat$time[j], 1,
                                 ifelse(mymat$group[i]>mymat$group[j] & mymat$time[i]>mymat$time[j],1,0))
            }
        if (mymat$time[j]<maxt & mymat$time[i]<maxt){
          temp1 <- ifelse(mymat$group[i]<mymat$group[j]& mymat$event[i]==e & mymat$event[j]==e &
                          mymat$time[i]<maxt, 1, 0)
          temp2 <- ifelse(mymat$group[i]<mymat$group[j] & mymat$time[i]<mymat$time[j]& mymat$time[i]<maxt &
                            mymat$event[i]==e & mymat$event[j]==e, 1, 0)
          denom = denom + temp1; num = num + temp2;
          ####
          denom2 = denom2 + temp1/(mymat$w2[i]^2); num2 = num2 + temp2/(mymat$w2[i]^2)
          #####
          denom3 = denom3 + temp1/(mymat$w[i]*mymat$w[j]); num3 = num3 + temp2/(mymat$w[i]*mymat$w[j])
          #####
          denom4 = denom4 + temp1/(mymat$w[i]*mymat$w2[j]); num4 = num4 + temp2/(mymat$w[i]*mymat$w2[j])
          #####
          denom5 = denom5 + temp1/(mymat$w2[i]*mymat$w[j]); num5 = num5 + temp2/(mymat$w2[i]*mymat$w[j])
          #####
          denom6 = denom6 + temp1/(mymat$w[i]^2); num6 = num6 + temp2/(mymat$w[i]^2)
          #####
          #denom7 = denom7 + temp1/(mymat$w2[i]*mymat$w[j]); num5 = num5 + temp2/(mymat$w2[i]*mymat$w[j])
          }
        }
      }
    }
    c_h[e] <- cij/Q; c_p[e] <- num/denom; c_u[e] <- num2/denom2;
    c_w1[e] <- num3/denom3; c_w2[e] <- num4/denom4; c_w3[e] <- num5/denom5; c_u2[e] = num6/denom6
  }
  names(c_h) = names(c_p) = names(c_u) = names(c_u2) = names(c_w1) = names(c_w2) = names(c_w3) = rlabels

  return(list(c_h=c_h, c_p = c_p, c_u=c_u, c_w1, c_w2, c_w3))
}


