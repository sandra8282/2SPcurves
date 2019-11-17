ROCsurv <- function(Time, Event, Group){
  km_placebo <- survfit(Surv(Time, Event) ~ 1, 
                        subset=(Group==0), type='kaplan-meier')
  km_drug <- survfit(Surv(Time, Event) ~ 1, 
                     subset=(Group==1), type='kaplan-meier')
  
  skm_p <- cbind(time = summary(km_placebo)$time, 
                 surv = summary(km_placebo)$surv) 
  skm_p <- cbind(skm_p, type = rep(0, nrow(skm_p)))
  skm_d <- cbind(time = summary(km_drug)$time, 
                 surv = summary(km_drug)$surv) 
  skm_d <- cbind(skm_d, type = rep(1, nrow(skm_d)))
  
  skm <- rbind(skm_p, skm_d)
  skm <- skm[order(skm[,1]),]
  
  plot(c(0,1), c(0, 1), type="n", xlab="", ylab="")
  title(main="ROC", xlab="Survival Placebo Group", 
        ylab="Survival Treatment Group", 
        cex.main = 1)
  area = 0
  for (i in 1:nrow(skm)) {
    if(i<2){
      coord_new = c(1, 1)
      #check if drug or placebo
      if (skm[i,3]==0) {#move horizontally
        coord_new2 = c(skm[i,2], 1)
        rect(xright = coord_new[1], ytop = coord_new[2], 
             xleft = coord_new2[1], ybottom = 0,
             col = "pink", border = "pink")
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {#move vertically
        coord_new2 = c(skm[i,2], 1)
      }
    } else {
      #check if drug or placebo
      if (skm[i,3]==0) {#move horizontally
        coord_new2 = c(skm[i,2], coord_new[2])
        rect(xright = coord_new[1], ytop = coord_new[2], 
             xleft = coord_new2[1], ybottom = 0,
             col = "pink", border = "pink")
        area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
      } else {#move vertically
        coord_new2 = c(coord_new[1], skm[i,2])
      }
      segments(coord_new[1], coord_new[2], 
               coord_new2[1], coord_new2[2], col="black")
      coord_new = coord_new2
    }
  }
  abline(c(0,1), col = "black", lty=2)
  area = unname(area)
  text(0.8, 0.15, paste("Area=", round(area,2), sep=""))  
  
  return(list(control_km = km_placebo, 
              treatment_km = km_drug, 
              area = area))
}

leukemia <- read.csv("ExampleLeukemia.csv")

library(survival)
leukemia_res <- ROCsurv(Time = leukemia$Time, 
        Event = leukemia$Event, 
        Group = leukemia$Group)

plot(km_drug, xlab="Time", ylab="Survival", 
     conf.int = FALSE, col="black")
lines(km_placebo, conf.int = FALSE, col="blue")
legend("topright", c("Treatment", "Placebo"), col = c("black", "blue"),
       lty = 1, cex=0.9, bg = "white", bty='n', seg.len = 0.7,
       x.intersp=0.9, y.intersp = 0.85)