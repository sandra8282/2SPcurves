notiesROC <- function(skm) {

  plot(c(0,1), c(0, 1), type="n", xlab="", ylab="")
  title(main="ROC", xlab="Control Group Survival",
        ylab="Treatment Group Survival",
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
  text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
       pos=2, cex = 1.25)

  return(list(control_KaplanMeier = km_placebo,
              treatment_KaplanMeier = km_drug,
              AUC = area))

}
