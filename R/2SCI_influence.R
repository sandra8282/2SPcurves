# cshaz <- function(times, df){
#   results = NULL
#   for (t in times[-1]){
#     Yi = ifelse(df$X_i >= t, 1, 0);
#     Ni1 = ifelse(df$X_i<= t & df$epdelta== 1 & Yi==1, 1, 0)
#     Ni2 = ifelse(df$X_i<= t & df$epdelta== 2 & Yi==1, 1, 0)
#     cshaz1 = sum(Ni1)/sum(Yi)
#     cshaz2 = sum(Ni2)/sum(Yi)
#     temp = c(t=t, AtRisk=sum(Yi), cshaz1=cshaz1, cshaz2=cshaz2)
#     results <- rbind(results, temp)
#   }
#   return(results)
# }
#
# part1Inf_ikt <- function(times, df, epsilon){
#   dInf_t = list()
#   for (t in times[-1]){
#     temp = NULL
#     index = which(times==t); tminus = times[index-1]
#     for (i in 1:nrow(df)){
#       subdf = df[i,]
#       subdf$t = t
#       subdf$Yi = ifelse(subdf$X_i >= t, 1, 0);
#       subdf$dNi1 = ifelse(subdf$X_i<= t & subdf$epdelta== epsilon & subdf$Yi, 1, 0)
#       subdf$Mit1 = (subdf$dNi1 - subdf$Yi*subdf$cshaz1)
#       subdf$dIit1 = subdf$survminus*subdf$Mit1*nrow(df)/subdf$Atrisk
#       subdf$i = i
#       temp = rbind(temp, subdf)
#     }
#     dInf_t[[index]] = temp
#   }
#   dInf_t[[1]]=0
#   names(dInf_t) = times
#   return(dInf_t)
# }
#
# part2Inf_ikt <- function(times, df, epsilon){
#   dInf_t = list()
#   for (t in times[-1]){
#     temp = NULL
#     index = which(times==t); tminus = times[index-1]
#     for (i in 1:nrow(df)){
#       subdf = df[i,]
#       subdf$t = t
#       subdf$Yi = ifelse(subdf$X_i >= t, 1, 0);
#       subdf$dNi = ifelse(subdf$X_i<= t & subdf$Yi, 1, 0)
#       subdf$Mit1 = (subdf$dNi - subdf$Yi*subdf$cs)
#       subdf$dIit1 = subdf$survminus*subdf$Mit1*nrow(df)/subdf$Atrisk
#       subdf$i = i
#       temp = rbind(temp, subdf)
#     }
#     dInf_t[[index]] = temp
#   }
#   dInf_t[[1]]=0
#   names(dInf_t) = times
#   return(dInf_t)
# }
#
# library(dplyr)
# ### survival and martingales
# df = data.frame(k = ctn2$trt, X_i = ctn2$time, epdelta = ctn2$status,
#                 delta_i = ifelse(ctn2$status==0, 0, 1))
# df$X_i =  ifelse(df$X_i==0, min(0.0001, min(df$X_i[df$X_i>0])/5), df$X_i)
# times <- unique(c(0, round(df$X_i, 4)))
# times <- sort(times, decreasing = FALSE)
# df0= df %>% filter(k==0) %>% arrange(X_i)
# df1= df %>% filter(k==1) %>% arrange(X_i)
#
# library(survival)
# surv0list <- summary(survfit(Surv(X_i, delta_i)~1, data=df0))
# surv1list <- summary(survfit(Surv(X_i, delta_i)~1, data=df1))
# survtab0 <- data.frame(t = surv0list$time, Atrisk = surv0list$n.risk, survival = surv0list$surv)
# survtab0$survminus <- lag(survtab0$survival, 1, default = 1)
# survtab1 <- data.frame(t = surv1list$time, Atrisk = surv1list$n.risk, survival = surv1list$surv)
# survtab1$survminus <- lag(survtab1$survival, 1, default = 1)
#
# TwoSCIres <- TwoSCI(time = ctn2$time, event = ctn2$status, group = ctn2$trt,
#        xlab="BMT CI", ylab="PBSC CI",
#        rlabels=c("Neutrophil Engraftment", "Pre-engraftment Mortality"), cex.axis = 1.5, cex.lab = 2,
#        lwd = 2, legend.inset=0.01, legend.cex=1.6, silent = FALSE)
# sfit <- TwoSCIres$cuminc
# dfnew <- data.frame(time = sfit$time,
#                   group = as.numeric(sfit$strata)-1,
#                   surv = sfit$pstate[,1],
#                   ci = sfit$pstate[,2:3])
# df0new <- dfnew[dfnew$group==0, ]
# df0new <- dfnew[dfnew$group==0, ]
# df0new <- left_join(df0, df0new, by=c("X_i"="time"))
#
# cshaz0 <- data.frame(cshaz(times, df0new))
# cshaz0 <- cshaz0[cshaz0$AtRisk!=0,]
# df0new <- left_join(df0new, cshaz0, by=c("X_i"="t"))
#
# df1new <- left_join(df1, survtab1, by=c("X_i"="t"))
# df1new <- df1new %>% filter(!is.na(Atrisk))
# cshaz1 <- data.frame(cshaz(times, df1new))
# cshaz1 <- cshaz0[cshaz1$AtRisk!=0,]
# df1new <- left_join(df1new, cshaz1, by=c("X_i"="t"))
#
# library (plyr)
# dInf_i01t <- dInf_ikt(times=times, df=df0new, epsilon=1)[-1]
# Inf_i01t <- rep(NA, nrow(df0new))
# for (i in 1:nrow(df0new)){
#   temp2 <- ldply(lapply(dInf_i01t, function(x) x[x$i==i,]), data.frame)
#   int_spline <- splinefun(
#     x = temp2$t,
#     y = temp2$dIit1,
#     method = "natural"
#   )
#   Inf_i01t[[i]] <- integrate(f = int_spline, lower = 0, upper = max(temp2$t))$value
# }
#
# dInf_i11t <- dInf_ikt(times=times, df=df1new, epsilon=1)[-1]
# Inf_i11t <- rep(NA, nrow(df1new))
# for (i in 1:nrow(df1new)){
#   temp2 <- ldply(lapply(dInf_i11t, function(x) x[x$i==i,]), data.frame)
#   int_spline <- splinefun(
#     x = temp2$t,
#     y = temp2$dIit1,
#     method = "natural"
#   )
#   Inf_i11t[[i]] <- integrate(f = int_spline, lower = 0, upper = max(temp2$t))$value
# }
#
# dInf_i02t <- dInf_ikt(times=times, df=df0new, epsilon=2)[-1]
# Inf_i02t <- rep(NA, nrow(df1new))
# for (i in 1:nrow(df0new)){
#   temp2 <- ldply(lapply(dInf_i02t, function(x) x[x$i==i,]), data.frame)
#   int_spline <- splinefun(
#     x = temp2$t,
#     y = temp2$dIit1,
#     method = "natural"
#   )
#   Inf_i02t[[i]] <- integrate(f = int_spline, lower = 0, upper = max(temp2$t))$value
# }
#
# dInf_i12t <- dInf_ikt(times=times, df=df1new, epsilon=2)[-1]
#
# int_spline<- splinefun(
#   x = tab$t,
#   y = tab$Survival,
#   method = "natural"
# )
#
#
#
