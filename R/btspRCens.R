solvena <- function(indi, temp){
    if (temp$tienext[indi[1]-1] == 1){
        xcoords <- temp$ref_u[c(indi[1]-1,
                                indi[length(indi)]+1)]
        ycoords <- temp$R_u[c(indi[1]-1, indi[length(indi)]+1)]
        slope <- diff(ycoords)/diff(xcoords)
        intercept = ycoords[2] - slope*xcoords[2]
        temp$R_u[indi] = slope*temp$ref_u[indi]+intercept
    } else {
        temp$R_u[indi] = temp$R_u[indi[1]-1]}
    return(temp)
}

#' Estimate SE and CI for AUC based on 1000 sample bootstrap
#'
#' @param time passed from ROCsurv
#' @param event passed from ROCsurv
#' @param group passed from ROCsurv
#' @param method passed from ROCsurv
#' @param level passed from ROCsurv
#' @param B number of bootstrap resamples
#'
#' @return SE and CI for indicated level
#'
#' @importFrom stats sd
#' @importFrom graphics segments
#' @importFrom graphics rect
#' @importFrom stats na.omit
#' @import survival
#' @importFrom data.table data.table
#' @importFrom data.table setkeyv
#' @importFrom data.table setorder
#' @importFrom data.table :=
#' @keywords internal
#' @noRd

btsp <- function(res, maindat, B, level, xlab, ylab, main, cex.axis = cex.axis,
                 cex.lab = cex.lab, lty = lty, lwd = lwd) {

    plot(NULL, type="n", las=1,
       xlim=c(0,1), ylim = c(0, 1), #to make tight axis: xaxs="i", yaxs="i"
       xlab=xlab, ylab=ylab, main=main, cex.axis = cex.axis, cex.lab = cex.lab)

    conf.lev = (1-level)/2

     maindat = data.table(maindat)
     g0dat <- data.table(maindat[maindat$group==1,]);
     g1dat <- data.table(maindat[maindat$group==0,]);
     n0 <- nrow(g0dat); g0ids = 1:n0; g0dat$id <- g0ids;
     n1 <- nrow(g1dat); g1ids = 1:n1; g1dat$id <- g1ids;
     rownames(g0dat) = g0ids; rownames(g1dat) = g1ids;
     ref <- res$R_u
     colnames(ref) <- c("ref_u", "ref_Ru", "ref_tienext")

     bstAUC = rep(NA, B)
     bst2SS = vector(mode = "list", length = nrow(ref));
     denom4prop = num4prop = rep(NA, B)

     for (b in 1:B) {
       set.seed(b+7679)
       g0ids_b <- sample(g0ids, n0, replace = TRUE)
       g1ids_b <- sample(g1ids, n1, replace = TRUE)
       btspdat <- data.table(rbind(g0dat[g0ids_b, ], g1dat[g1ids_b, ]))
       KMests <- with(btspdat, getKMtab(time, event, group))
       temp <- completeROC(KMests[[1]], silent = TRUE)
       bstAUC[b] <- temp$AUC
       forplot <- temp$R_u
       lines(forplot[,1:2], col = "bisque")
       temp <- merge(ref[-1,], forplot, by.x = "ref_u", by.y = "u", all = TRUE)
       temp <- temp[order(temp$ref_u, decreasing = TRUE), ]
       inds1 <- which(is.na(temp$R_u)); l <- length(inds1)
       if (l>0){inds1 <- split(inds1, cumsum(c(1, diff(inds1) != 1)))
       temp$tienext[inds1[[length(inds1)]][1]-1]=0}
       while (l > 0) {
           inds <- which(is.na(temp$R_u))
           inds <- split(inds, cumsum(c(1, diff(inds) != 1)))
           temp <- lapply(inds, function(x) solvena(indi = x, temp = temp))
           temp <- temp[[1]]
           l <- length(which(is.na(temp$R_u)))
       }
       temp <- data.frame(rbind(rep(1,5), temp))
       temp <- temp[-which(is.na(temp$ref_Ru)),]
       bst2SS[[b]] <- temp[ , c(1,4)]
       num4prop[b] <- length(which(forplot[,1]>forplot[,2]))
       denom4prop[b] <- nrow(forplot)
     }

     prop <- sum(num4prop)/sum(denom4prop)
     xunik <- unique(ref[,1])
     names(bst2SS) = paste("s", 1:B, sep="")
     se = CIlow = CIup = rep(NA, length(xunik))

     for (xucount in 1:length(xunik)){
       temp2 <- unlist(lapply(bst2SS, function(x) x$R_u[which(x$ref_u==xunik[xucount])]))
       se[xucount] <- sd(temp2)
       CIlow[xucount] <- sort(temp2)[(conf.lev)*length(temp2)]
       CIup[xucount] <- sort(temp2)[(1-conf.lev)*length(temp2)]
     }

     boot <- data.frame(u=xunik, se_Ru = se, CIlow_Ru = CIlow, CIup_Ru = CIup)
     R_u <- merge(ref, boot, by.x = "ref_u", by.y = "u", all.x = TRUE)
     colnames(R_u)[1:2] = c("u", "Ru")
     R_u <- R_u[order(R_u$u, decreasing = TRUE), ]
     lines(R_u$u, R_u$Ru, lty = lty[1], lwd = lwd)
     lines(R_u$u, R_u$CIlow_Ru, lty = 2, col = "darkred", lwd = lwd)
     lines(R_u$u, R_u$CIup_Ru, lty = 2, col = "darkred", lwd = lwd)
     abline(c(0,1), col = "grey", lty=1, lwd = lwd-0.25)

     AUC <- c(AUC = res$AUC,
            CIlow = sort(bstAUC)[(conf.lev)*length(bstAUC)],
            CIup = sort(bstAUC)[(1-conf.lev)*length(bstAUC)])

     return(list(AUC= AUC, R_u = R_u, dprop = prop))
 }
