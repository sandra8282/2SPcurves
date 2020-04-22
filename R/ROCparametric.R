getSKM4fit <- function(time, fitsurv, group) {
  fitskm <- cbind(time, fitsurv, group)
  fitskm <- fitskm[!duplicated(fitskm),]
  ties_check <- unique(table(fitskm[,1])) #how many times are repeated
  if (length(ties_check) > 1) {
    ties_times = fitskm[duplicated(fitskm[,1]),1]
    ties_ind <- rep(0, nrow(fitskm))
    ties_ind[which(fitskm[,1] %in% ties_times)]=1
  } else {ties_ind <- rep(0, nrow(fitskm))}
  fitskm = cbind(fitskm, ties_ind)
  fitskm <- fitskm[order(fitskm[,1], fitskm[,3]),]
  return(fitskm)
}

makegmat1 <- function(vect, b0, b1, sig) {
  eb <- exp(-(b0 + b1)/sig); eb0 <- exp(-b0/sig); eb1 <- exp(-b1/sig);
  uinv <- 1/vect[2]; t <- vect[1]; D1<-vect[3]; D2<-vect[4];
  gmat1 <- ((uinv^2)*((t^(1/sig))*eb0))/(sig*D1)
  gmat2 <- (uinv - 1)/sig
  gmat3 <- ((uinv^2 * t^(1/sig) * eb0 * (log(t) - b0)^2) / (sig^2 * D1) ) - (b1/sig^2)*(uinv - 1)
  gmat = cbind((eb1/D2)*gmat1, (eb1/D2)*gmat2, (eb1/D2)*gmat3)
  return(gmat)
}

makegmat2 <- function(vect, b0, b1, sig) {
  gmat1 <- dnorm(vect[2] - b0/sig) / (sig*vect[3])
  gmat2 <- 1/sig
  gmat3 <- -(b1/sig) - ((b0*dnorm(vect[2]-b0/sig))/(sig^2 * vect[3]))
  gmat = cbind(vect[1]*gmat1, vect[1]*gmat2, vect[1]*gmat3)
  return(gmat)
}

#' ROC when survival goes to 0 for either group
#'
#' @param time Numeric or character vector of subject's unique identifier (i).
#' @param event Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param group Vector with the lengths of time spent in event of Type I for individual i in episode j.
#' @param dist String indicating the distribution.
#'
#' @return A plot of the ROC curve and an ROCsurv object containing:
#' \itemize{
#'  \item A survfit object for the treatment group.
#'  \item A survfit object for the control group.
#'  \item The area under the curve for the ROC in the given plot.
#' }
#'
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @importFrom stats na.omit
#' @import survival
#' @importFrom stats vcov
#' @importFrom stats dnorm
#' @importFrom stats qnorm
#' @importFrom utils capture.output
#'
#' @export

ROCparametric <- function(time, event, group, dist="lognormal") {

  KMres <- getKMtab(time, event, group)
  skm <- KMres[[1]]
  forplot = get4plot(skm)

  #obtain coefficients and parameters to calculate survival fcts
  fit = survreg(Surv(time, event) ~ group, dist=dist)
  mu = fit$coefficients[1]; gamma = fit$coefficients[2]; sigma = fit$scale;
  lambda = exp(-unname(mu/sigma)); beta = -unname(gamma)/sigma; alpha = 1/sigma;
  covmat <- vcov(fit)

  #calculate survival fcts
  if (dist=="loglogistic") {
    fitsurv <- (1 + lambda*exp(beta*group)*(time^alpha))^(-1)
    fitskm <- getSKM4fit(time, fitsurv, group)
    forplotfit = get4plot(fitskm)
    newsurv0 <- seq(min(forplotfit[,1]), 0.001, length.out = round(nrow(forplotfit)*0.1, 0))
    forplotfit <- rbind(forplotfit, cbind(newsurv0, rep(0, length(newsurv0))))
    forplotfit[,2] = 1 / (1 + (exp(- gamma/sigma) * (1/forplotfit[,1] - 1)))
    newt <- (1/newsurv0 - 1)^sigma*exp(mu)
    t2 <- c(fitskm[,1], newt)
    D1 <- (1 + (t2*exp(-mu))^alpha)^2
    D2 <- apply(forplotfit, 2, function(x) (1 + (1/x - 1)*exp(beta)^2))[,1]
    D <- cbind(t = t2, u = forplotfit[,1], D1, D2)
    gmat <- apply(D, 1, function(x) makegmat1(x, b0=mu, b1=gamma, sig=sigma))
    gmat[,c(1, ncol(gmat))]<- c(rep(0,3), gmat[, ncol(gmat)-1])
    var <- apply(gmat, 2, function(x) matrix(x, ncol=3) %*% covmat %*% t(matrix(x, ncol=3)))
    se <- sqrt(var)
    forplotfit <- cbind(forplotfit, se=se)
    low <- forplotfit[,2] - 1.96*se
    up <- forplotfit[,2] + 1.96*se; up <- ifelse(up > 1, 1, up);
    forplotfit <- cbind(forplotfit, lower95 = low, upper95 = up,
                        ratio = forplotfit[,3]/forplotfit[,2])

  }

  if (dist == "lognormal") {
    fitsurv <- 1 - pnorm((log(time)- mu - gamma*group)/sigma)
    fitskm <- getSKM4fit(time, fitsurv, group)
    forplotfit = get4plot(fitskm)
    newsurv0 <- seq(min(forplotfit[,1]), 1e-22, length.out = round(nrow(forplotfit)*0.15, 0))
    forplotfit <- rbind(forplotfit, cbind(newsurv0, rep(0, length(newsurv0))))
    forplotfit[,2] = 1 - pnorm(qnorm(1 - forplotfit[,1]) - gamma/sigma)
    newt <- exp(qnorm(1-newsurv0) + mu/sigma)
    yi <- log(c(fitskm[,1], newt))
    common <- dnorm(qnorm(1-forplotfit[,1]) - gamma/sigma)
    D <- dnorm(qnorm(1-forplotfit[,1]))
    Dmat <- cbind(common = common, yi = yi, D = D)
    gmat <- apply(Dmat, 1, function(x) makegmat2(x, b0=mu, b1=gamma, sig=sigma))
    gmat[,c(1, ncol(gmat))]<- c(rep(0,3), gmat[, ncol(gmat)-1])
    var <- apply(gmat, 2, function(x) matrix(x, ncol=3) %*% covmat %*% t(matrix(x, ncol=3)))
    se <- sqrt(var)
    forplotfit <- cbind(forplotfit, se=se)
    index <- c(which(forplotfit[,1]==1 & forplotfit[,2]==1), nrow(forplotfit))
    low <- c(rep(1, length(index)-1), forplotfit[-index,2] - 1.96*se[-index], 0)
    up <- c(rep(1, length(index)-1), forplotfit[-index,2] + 1.96*se[-index], 0)
    up <- ifelse(up>1, 1, up)
    forplotfit <- cbind(forplotfit, lower95 = low, upper95 = up,
                        ratio = forplotfit[,3]/forplotfit[,2])
  }

   area = 0
   plot(NULL, type="n", las=1,
         xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
         xlab="Control Group Survival", ylab="Treatment Group Survival",
         cex.axis = 1.25, cex.lab = 1.25)

   for (k in 2:nrow(forplotfit)) {
     coord_new = unname(forplotfit[k-1,])
     coord_new2 = unname(forplotfit[k,])
     #figure out areas and shading
     if (forplotfit[k,2]==forplotfit[k-1,2]) {#move horizontally
       rect(xright = coord_new[1], ytop = coord_new[2],
            xleft = coord_new2[1], ybottom = 0,
            col = "pink", border = "pink")
       area = area + (coord_new[1] - coord_new2[1])*(coord_new[2])
     } else {
       if (forplotfit[k,1]!=forplotfit[k-1,1] & forplotfit[k,2]!=forplotfit[k-1,2]){
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
   }

   points(forplot[,1], forplot[,2], col = "grey50", cex = 0.75)
   lines(forplotfit[,1], forplotfit[,2], lty=1)
   if (dist == "loglogistic") {
        sub4plotfit <- forplotfit[which(forplotfit[,6] < 0.5),]
        lines(sub4plotfit[,1], sub4plotfit[,4], lty=2)
        lines(sub4plotfit[,1], sub4plotfit[,5], lty=2)
   }
   if (dist == "lognormal") {
     lines(forplotfit[,1], forplotfit[,4], lty=2)
     lines(forplotfit[,1], forplotfit[,5], lty=2)
   }

   abline(c(0,1), col = "red", lty=2)
   text(x=0.99, y=0.05, labels = paste("AUC=", round(area,2), sep=""),
        pos=2, cex = 1)

   return(list(fit, area))

}
