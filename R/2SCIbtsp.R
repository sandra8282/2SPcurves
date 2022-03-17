solvenaC <- function(indi, temp2){
        if (temp2$tienext[indi[1]-1] == 1){
            xcoords <- temp2$ref_u[c(indi[length(indi)]+1, indi[1]-1)]
            ycoords <- temp2$y[c(indi[length(indi)]+1, indi[1]-1)]
            slope <- diff(ycoords)/diff(xcoords)
            intercept = ycoords[2] - slope*xcoords[2]
            temp2$y[indi] = slope*temp2$ref_u[indi]+intercept
        } else {temp2$y[indi] = temp2$y[indi[1]-1]}
    return(temp2)
}

#' Estimate SE and CI for AUC based on 1000 sample bootstrap
#'
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

btsp2SCI <- function(res, maindat, maxt, nrisktypes, B, level, xlab, ylab, rlabels, main,
                     cex.axis = cex.axis, cex.lab = cex.lab, lty = 1, lwd = lwd,
                     bst_c=NULL, cindex=NULL) {
    ptypes=c(NA, 1)
    darkcolopt =c("black", "grey29")
    lwdopt = c(lwd, 1)

    enames = c("censored", rlabels)

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

     ref <- res
     changecolnam <- function(list_element){
         colnames(list_element) <- c("ref_u", "ref_Cu", "ref_tienext")
         return(list_element)}
     ref <- lapply(ref, function(x) changecolnam(x))

     # for (refi in 1:nrisktypes){
     #     lines(ref[[refi]][,1],ref[[refi]][,2], lty = 2, lwd = lwd) #darkcolopt[[refi]]
     # }

     toreturn = list()
     bstCindex = matrix(rep(NA, B*2), ncol=2)
     for (i in 1:nrisktypes){
         bst2SS = vector(mode = "list", length = B);
         for (b in 1:B) {
             set.seed(b+7679)
             g0ids_b <- sample(g0ids, n0, replace = TRUE)
             g1ids_b <- sample(g1ids, n1, replace = TRUE)
             btspdat <- data.table(rbind(g0dat[g0ids_b, ], g1dat[g1ids_b, ]))
             #####
             if (bst_c==TRUE){bstCindex[b,1] <- comprsk_c(btspdat, rlabels, maxt)[1]}
             #####
             btspdat$eventf <- factor(btspdat$event, 0:nrisktypes, labels = enames)
             fit <-  survfit(Surv(btspdat$time, btspdat$eventf) ~ btspdat$group)
             sfit <- summary(fit)
             #get C(u)
             skm <- data.frame(time = sfit$time, group = as.numeric(sfit$strata)-1, ci = sfit$pstate[,2:3])
             skm <- skm[order(skm[,1], skm[,2]),]
             ties_check <- unique(table(skm[,1]))
             if (length(ties_check) > 1) {
                 ties_times = skm[duplicated(skm[,1]),1]
                 ties_ind <- rep(0, nrow(skm))
                 ties_ind[which(skm[,1] %in% ties_times)]=1
             } else {ties_ind <- rep(0, nrow(skm))}
             skm = cbind(skm, ties_ind)
             list_4plot = bstres = list()
             for (skmi in (1:nrisktypes)){
               skmires = matrix(as.numeric(as.matrix(distinct(skm[, c(1, 2+skmi, 2)]))),
                                ncol=3)
               ties_check <- duplicated(skmires[,1])
               if (length(ties_check) > 1) {
                 ties_times = skmires[duplicated(skmires[,1]),1]
                 ties_ind <- rep(0, nrow(skmires))
                 ties_ind[which(skmires[,1] %in% ties_times)]=1
               } else {ties_ind <- rep(0, nrow(skmires))}
               skmires = cbind(skmires, ties_ind)
               list_4plot[[skmi]] = get4plotCumInc(skmires)
             }
             for (ploti in (1:nrisktypes)){
                 #lines(list_4plot[[ploti]][,1:2], type="s", lty = 1, lwd = 2, col = colopt[ploti])
                 # lines(ref[[ploti]][,1:2])
                 new <- data.frame(ref[[ploti]])
                 outs <- unique(c(which(new$ref_u==0&new$ref_Cu==0), which(duplicated(new))))
                 new <- new[-outs,]
                 new2 <- data.frame(list_4plot[[ploti]])
                 outs <- unique(c(which(new2$x==0&new2$y==0), which(duplicated(new2))))
                 new2 <- new2[-outs,]
                 temp <- merge(new, new2, by.x = "ref_u", by.y = "x", all = TRUE)
                 temp <- temp[order(temp$ref_u, decreasing = FALSE), ]
                 if (is.na(temp$y[1])){temp$y[1]=0; temp$tienext[1]=0}
                 # maxy = which(temp$y==max(temp$y, na.rm = TRUE))
                 # outs <- which(temp$ref_u>temp$ref_u[maxy])
                 # temp = temp[-outs,]
                 inds1 <- which(is.na(temp$y)); l <- length(inds1)
                 if (l>0){
                     inds1 <- split(inds1, cumsum(c(1, diff(inds1) != 1)))
                      temp$tienext[inds1[[length(inds1)]][1]-1]=0}
                 while (l > 0) {
                     inds <- which(is.na(temp$y))
                     inds <- split(inds, cumsum(c(1, diff(inds) != 1)))
                     for (indsi in 1:length(inds)){
                         temp <- solvenaC(indi = inds[[indsi]], temp2 = temp)
                     }
                     # temp <- lapply(inds, function(x) solvenaC(indi = x, temp2 = temp))
                     # temp <- temp[[1]]
                     l <- length(which(is.na(temp$y)))
                     #print(l)
                 }
                 temp <- data.frame(rbind(rep(0,5), temp))
                 temp <- temp[-which(is.na(temp$ref_Cu)),]
                 bstres[[ploti]]=temp
             }
             bst2SS[[b]] <- bstres
         }
     }
        names(bst2SS) = paste("CUMINC", 1:B, sep="")
        xunik = se = CIlow = CIup = boot = C_u = list()
        for (uniki in 1:nrisktypes){
          xunik[[uniki]] <- unique(ref[[uniki]][,1])
          se[[uniki]] = CIlow[[uniki]] = CIup[[uniki]] = rep(NA, length(xunik[[uniki]]))
          for (xucount in 1:length(xunik[[uniki]])){
              temp2 <- unlist(lapply(bst2SS, function(x)
                          x[[uniki]]$y[which(x[[uniki]]$ref_u==xunik[[uniki]][xucount])]))
              se[[uniki]][xucount] <- sd(temp2)
              CIlow[[uniki]][xucount] <- sort(temp2)[(conf.lev)*length(temp2)]
              CIup[[uniki]][xucount] <- sort(temp2)[(1-conf.lev)*length(temp2)]
          }
          boot[[uniki]] <- data.frame(u=xunik[[uniki]], se_Cu = se[[uniki]],
                                      CIlow_Cu = CIlow[[uniki]] , CIup_Cu = CIup[[uniki]] )
          C_u[[uniki]] <- merge(ref[[uniki]], boot[[uniki]], by.x = "ref_u", by.y = "u", all.x = TRUE)
          colnames(C_u[[uniki]])[1:2] = c("u", "Cu")
          C_u[[uniki]] <- C_u[[uniki]][order(C_u[[uniki]]$u, decreasing = FALSE), ]
          # polygon(x=c(C_u[[uniki]]$u, rev(C_u[[uniki]]$u)),
          #         y=c(C_u[[uniki]]$CIlow_Cu, rev(C_u[[uniki]]$CIup_Cu)), col = "grey", border = NA)

          #lines(C_u[[uniki]]$u, C_u[[uniki]]$Cu, lty = 2, lwd = lwd)#, col=darkcolopt[[uniki]]

        }


         abline(c(0,1), lty=1, lwd = lwd-0.5)
         for (refi in 1:nrisktypes){
             uniki=refi
             refindlow = which(diff(C_u[[uniki]]$CIlow_Cu)<0)
             refindup = which(diff(C_u[[uniki]]$CIup_Cu)<0)
             C_u[[uniki]]$CIlow_Cu[refindlow]=NA
             C_u[[uniki]]$CIup_Cu[refindup]=NA
             while (length(refindlow)>=1){
                 C_u[[uniki]]$CIlow_Cu[refindlow] = C_u[[uniki]]$CIlow_Cu[refindlow+1]
                 refindlow = which(diff(C_u[[uniki]]$CIlow_Cu)<0)
             }
             while (length(refindup)>=1){
                 C_u[[uniki]]$CIup_Cu[refindup] = C_u[[uniki]]$CIup_Cu[refindup+1]
                 refindup = which(diff(C_u[[uniki]]$CIup_Cu)<0)
             }

             lines(C_u[[uniki]]$u, C_u[[uniki]]$CIlow_Cu, lty = 1, lwd = lwdopt[uniki], col=darkcolopt[uniki])
             points(C_u[[uniki]]$u, C_u[[uniki]]$CIlow_Cu, pch = ptypes[uniki], cex=0.7, col=darkcolopt[uniki])
             lines(C_u[[uniki]]$u, C_u[[uniki]]$CIup_Cu, lty = 1, lwd = lwdopt[uniki], col=darkcolopt[uniki])
             points(C_u[[uniki]]$u, C_u[[uniki]]$CIup_Cu, pch = ptypes[uniki], cex=0.7, col=darkcolopt[uniki])
             lines(ref[[refi]][,1],ref[[refi]][,2], lty = 1, lwd = lwdopt[uniki], col=darkcolopt[refi])
             points(ref[[refi]][,1],ref[[refi]][,2], pch = ptypes[refi], cex=0.7, col=darkcolopt[refi])
         }
         legend("topleft",
                legend = rlabels, pch = ptypes, lty = 1, cex=1.5, lwd = lwdopt,
                col = darkcolopt, bty="n", x.intersp=0.9, y.intersp = 0.85,)

         if (bst_c==TRUE){
             bstCindex <- bstCindex[,-2]
             finalres <- list(C_u = C_u,
                         Cindex = c(cindex = cindex,
                                    se_cindex = sd(bstCindex),
                                    CIlow = sort(bstCindex)[(conf.lev)*length(bstCindex)],
                                    CIup = sort(bstCindex)[(1-conf.lev)*length(bstCindex)]))
         } else {finalres <- C_u}

return(finalres)
}


