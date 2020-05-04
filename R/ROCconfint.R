# makegmat1 <- function(vect, b0, b1, sig) {
#   eb <- exp(-(b0 + b1)/sig); eb0 <- exp(-b0/sig); eb1 <- exp(-b1/sig);
#   uinv <- 1/vect[2]; t <- vect[1]; D1<-vect[3]; D2<-vect[4];
#   gmat1 <- ((uinv^2)*((t^(1/sig))*eb0))/(sig*D1)
#   gmat2 <- (uinv - 1)/sig
#   gmat3 <- ((uinv^2 * t^(1/sig) * eb0 * (log(t) - b0)^2) / (sig^2 * D1) ) - (b1/sig^2)*(uinv - 1)
#   gmat = cbind((eb1/D2)*gmat1, (eb1/D2)*gmat2, (eb1/D2)*gmat3)
#   return(gmat)
# }
#
# makegmat2 <- function(vect, b0, b1, sig) {
#   gmat1 <- dnorm(vect[2] - b0/sig) / (sig*vect[3])
#   gmat2 <- 1/sig
#   gmat3 <- -(b1/sig) - ((b0*dnorm(vect[2]-b0/sig))/(sig^2 * vect[3]))
#   gmat = cbind(vect[1]*gmat1, vect[1]*gmat2, vect[1]*gmat3)
#   return(gmat)
# }
#
# #LOGLOGIS
# newt <- (1/newsurv0 - 1)^sigma*exp(mu)
# t2 <- c(fitskm[,1], newt)
# D1 <- (1 + (t2*exp(-mu))^alpha)^2
# D2 <- apply(forplotfit, 2, function(x) (1 + (1/x - 1)*exp(beta)^2))[,1]
# D <- cbind(t = t2, u = forplotfit[,1], D1, D2)
# gmat <- apply(D, 1, function(x) makegmat1(x, b0=mu, b1=gamma, sig=sigma))
# gmat[,c(1, ncol(gmat))]<- c(rep(0,3), gmat[, ncol(gmat)-1])
# var <- apply(gmat, 2, function(x) matrix(x, ncol=3) %*% covmat %*% t(matrix(x, ncol=3)))
# se <- sqrt(var)
# forplotfit <- cbind(forplotfit, se=se)
# low <- forplotfit[,2] - 1.96*se
# up <- forplotfit[,2] + 1.96*se; up <- ifelse(up > 1, 1, up);
# forplotfit <- cbind(forplotfit, lower95 = low, upper95 = up,
#                     ratio = forplotfit[,3]/forplotfit[,2])
# index <- which(forplotfit[,6] > 0.5)
# forplotfit[index,4] <- NA
# forplotfit[index,5] <- NA
# for (i in nrow(na.omit(forplotfit)):1) {
#   diff <- forplotfit[i, 5] - forplotfit[i-1, 5]
#   if (diff < 0) {
#     changeindex = i
#     break
#   }
# }
# forplotfit[changeindex:nrow(forplotfit), 4:5] <- NA
#
# ##LOGNORM
# newt <- exp(qnorm(1-newsurv0) + mu/sigma)
# yi <- log(c(fitskm[,1], newt))
# common <- dnorm(qnorm(1-forplotfit[,1]) - gamma/sigma)
# D <- dnorm(qnorm(1-forplotfit[,1]))
# Dmat <- cbind(common = common, yi = yi, D = D)
# gmat <- apply(Dmat, 1, function(x) makegmat2(x, b0=mu, b1=gamma, sig=sigma))
# gmat[,c(1, ncol(gmat))]<- c(rep(0,3), gmat[, ncol(gmat)-1])
# var <- apply(gmat, 2, function(x) matrix(x, ncol=3) %*% covmat %*% t(matrix(x, ncol=3)))
# se <- sqrt(var)
# forplotfit <- cbind(forplotfit, se=se)
# index <- c(which(forplotfit[,1]==1 & forplotfit[,2]==1), nrow(forplotfit))
# low <- c(rep(1, length(index)-1), forplotfit[-index,2] - 1.96*se[-index], 0)
# up <- c(rep(1, length(index)-1), forplotfit[-index,2] + 1.96*se[-index], 0)
# up <- ifelse(up>1, 1, up)
# forplotfit <- cbind(forplotfit, lower95 = low, upper95 = up,
#                     ratio = forplotfit[,3]/forplotfit[,2])
#
# #PLOTTING
# area = 0
# plot(NULL, type="n", las=1,
#      xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
#      xlab="Control Group Survival", ylab="Treatment Group Survival",
#      cex.axis = 1.25, cex.lab = 1.25)
#
# for (k in 2:nrow(forplotfit)) {
#   coord_new = unname(forplotfit[k-1,])
#   coord_new2 = unname(forplotfit[k,])
#   #figure out areas and shading
#   if (forplotfit[k,2]==forplotfit[k-1,2]) {#move horizontally
#     rect(xright = coord_new[1], ytop = coord_new[2],
#          xleft = coord_new2[1], ybottom = 0,
#          col = "pink", border = "pink")
#     area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
#   } else {
#     if (forplotfit[k,1]!=forplotfit[k-1,1] & forplotfit[k,2]!=forplotfit[k-1,2]){
#       #area and shading for diagonal
#       rect(xright = coord_new[1], ytop = coord_new2[2],
#            xleft = coord_new2[1], ybottom = 0,
#            col = "pink", border = "pink")
#       area_rectang = (coord_new[1] - coord_new2[1])*(coord_new2[2])
#       polygon(x=c(coord_new[1], coord_new[1], coord_new2[1]),
#               y=c(coord_new[2], coord_new2[2], coord_new2[2]),
#               col = "pink", border = "pink")
#       area_triang = 0.5 * (coord_new[1] - coord_new2[1]) * (coord_new[2] - coord_new2[2])
#       area = area + area_rectang + area_triang
#     }
#   }
# }
#
# points(forplot[,1], forplot[,2], col = "grey50", cex = 0.75)
# lines(forplotfit[,1], forplotfit[,2], lty=1)
# lines(forplotfit[,1], forplotfit[,4], lty=2)
# lines(forplotfit[,1], forplotfit[,5], lty=2)
#
#
# colnames(forplotfit) <- c("u", "R(u)", "SE[R(u)]", "lower.95", "upper.95", "ratio")
# abline(c(0,1), col = "red", lty=2)
# text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
#      pos=2, cex = 1)
