
#library(dplyr)
skm_main <- data.frame(time = sfit$time, n = sfit$n.risk[,1], et = sfit$n.event[,-1],
                  group = as.numeric(sfit$strata)-1, surv = sfit$pstate[,1],
                  ci = sfit$pstate[,2:3])
skm0 <- skm_main[skm_main$group==0,]
  if (skm0$time[1]!=0){rbind(c(0, skm0$n[1], 0, 0, 0, 1, 0, 0), skm0)}
  skm0$hz01 = skm0$et.1/skm0$n
  skm0$hz02 = skm0$et.2/skm0$n
  skm0$sminus0 = dplyr::lag(skm0$surv, 1, 1)
  skm0$ft01 = skm0$sminus0*skm0$hz01
  skm0$ft02 = skm0$sminus0*skm0$hz02
  skm0 <- skm0 %>% select(time, n, et.1, et.2, ci.1, ci.2, ft01, ft02)
  colnames(skm0) <- c("time", "Yi0", "Ni01", "Ni02", "Ft01", "Ft02", "ft01", "ft02")
skm1 <- skm_main[skm_main$group==1,]
  if (skm1$time[1]!=0){rbind(c(0, skm1$n[1], 0, 0, 1, 1, 0, 0), skm1)}
  skm1$hz11 = skm1$et.1/skm1$n
  skm1$hz12 = skm1$et.2/skm1$n
  skm1$sminus1 = dplyr::lag(skm1$surv, 1, 1)
  skm1 <- skm1 %>% select(time, n, et.1, et.2, ci.1, ci.2)
  colnames(skm1) <- c("time", "Yi1", "Ni11", "Ni12","Ft11", "Ft12")
skm_main <- full_join(skm0, skm1, by = "time")
skm_main <- skm_main %>% arrange(time)
for (i in 2:nrow(skm_main)){
  if (is.na(skm_main[i,2])){
    skm_main[i,2] = skm_main[i+1,2]
    skm_main[i,5:8] = skm_main[i-1,5:8]
    skm_main[i,3:4] = 0
    }
  if (is.na(skm_main[i,9])){
    skm_main[i,9] = skm_main[i+1,9]
    skm_main[i,c(12:13)] = skm_main[i-1,c(12:13)]
    skm_main[i,10:11] = 0
  }
}
skm_main$term1 <- (skm_main$Ft11 - skm_main$Ft01)*skm_main$ft01
skm_main$term2 <- (skm_main$Ft02 - skm_main$Ft12)*skm_main$ft02
  int_spline <- splinefun(
    x = skm_main$time,
    y = skm_main$term1,
    method = "natural"
  )
  abcd1 <- integrate(f = int_spline, lower = 0, upper = max(skm_main$time))$value
  int_spline <- splinefun(
    x = skm_main$time,
    y = skm_main$term2,
    method = "natural"
  )
  abcd2 <- integrate(f = int_spline, lower = 0, upper = max(skm_main$time))$value
ABCD = abcd1+abcd2
