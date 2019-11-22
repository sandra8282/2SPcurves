
with_tiesROC <- function(skm) {
plot(c(0,1), c(0, 1), type="n", xlab="", ylab="")
title(main="ROC", xlab="Control Group Survival",
      ylab="Treatment Group Survival",
      cex.main = 1)

x = y = c(1, rep(NA, nrow(skm)))
i=1
while (i < nrow(skm)) {
    if (skm[i,3]==0 & skm[i,4]==0) {#placebo no ties
      x[i+1]=skm[i,2]
      y[i+1]= y[i]
    }
    if (skm[i,3]==1 & skm[i,4]==0) {#drug no ties
        x[i+1]=x[i]
        y[i+1]=skm[i,2]
    }
    if (skm[i,4]==1) {#tie
        if (skm[i,2]==0){
          x[i+1]=skm[i,2]
          y[i+1]=skm[i+1,2]
        } else{
          x[i+1]=skm[i,2]
          y[i+1]=skm[i+1,2]
        }
      i=i+1
    }
  i=i+1
}

forplot <- na.omit(cbind(x,y))
area <- 0

for (k in 2:nrow(forplot)) {
  coord_new = forplot[k-1,]
  coord_new2 = forplot[k,]
  #figure out areas and shading
  if (forplot[k,2]==forplot[k-1,2]) {#move horizontally
      rect(xright = coord_new[1], ytop = coord_new[2],
           xleft = coord_new2[1], ybottom = 0,
           col = "pink", border = "pink")
    area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
    } else {
        if (forplot[k,1]!=forplot[k-1,1] & forplot[k,2]!=forplot[k-1,2]){
         #area and shading for diagonal
            rect(xright = coord_new[1], ytop = coord_new2[2],
                 xleft = coord_new2[1], ybottom = 0,
                 col = "pink", border = "pink")
            area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
            polygon(x=c(coord_new[1], coord_new[1], coord_new2[1]),
                    y=c(coord_new[2], coord_new2[2], coord_new2[2]),
                    col = "pink", border = "pink")
            area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
            area = area + area_rectang + area_triang
        }
      }
  segments(coord_new[1], coord_new[2],
             coord_new2[1], coord_new2[2], col="black")
}
abline(c(0,1), col = "black", lty=2)
area = unname(area)
text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
     pos=2, cex = 1)

return(area)
}
