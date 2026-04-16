#' Estimate SE and CI for AUC based on 1000 sample bootstrap
#'
#' @param control_pf passed from TwoSSPicens
#' @param trt_pf passed from TwoSSPicens
#' @param maindat passed from TwoSSPicens
#' @param B number of bootstrap resamples
#' @param level passed from TwoSSPicens
#' @param xlab passed from TwoSSPicens
#' @param ylab passed from TwoSSPicens
#' @param main passed from TwoSSPicens
#' @param cex.axis passed from TwoSSPicens
#' @param cex.lab passed from TwoSSPicens
#' @param legend.inset passed from TwoSSPicens
#' @param legend.cex passed from TwoSSPicens
#' @param lty passed from TwoSSPicens
#' @param lwd passed from TwoSSPicens
#' @param iterations passed from TwoSSPicens
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

btspICEN <- function(control_pf, trt_pf, maindat, B, level, xlab, ylab, main, cex.axis = cex.axis,
                 cex.lab = cex.lab, lty = lty, lwd = lwd, iterations) {

  res <- getIGroc(npmle_0 = control_pf,
                  npmle_1 = trt_pf,
                  xlab, ylab, main, cex.axis=1.25, cex.lab = 1.5, lwd,
                  shade = TRUE)

    conf.lev = (1-level)/2

     maindat = data.table(maindat)
     g0dat <- data.table(maindat[maindat$group==1,]);
     g1dat <- data.table(maindat[maindat$group==0,]);
     n0 <- nrow(g0dat); g0ids = 1:n0; g0dat$id <- g0ids;
     n1 <- nrow(g1dat); g1ids = 1:n1; g1dat$id <- g1ids;
     rownames(g0dat) = g0ids; rownames(g1dat) = g1ids;
     ref <- res
     colnames(ref) <- c("ref_u", "ref_Ru")

     #bstAUC = rep(NA, B)
     bst2SS = vector(mode = "list", length = nrow(ref));

     for (b in 1:B) {
       check = 1
       while (check>0){
         set.seed(b+764579+check)
         g0ids_b <- sample(g0ids, n0, replace = TRUE)
         g1ids_b <- sample(g1ids, n1, replace = TRUE)
         btspdat <- data.table(rbind(g0dat[g0ids_b, ], g1dat[g1ids_b, ]))
         # EM calculations for control group #################################################
         NPMLE.control<-EMICM(btspdat[btspdat$group==0,c("left", "right")], maxiter = iterations)
         bcontrol_pf <- data.frame(L = NPMLE.control$intmap[1,],
                                  R = NPMLE.control$intmap[2,],
                                  drop = NPMLE.control$pf)
         bcontrol_pf <- bcontrol_pf[bcontrol_pf$drop!=0,]
         bcontrol_pf$cumdrop = cumsum(bcontrol_pf$drop)

         # EM calculations for treatment group #################################################
         NPMLE.trt<-EMICM(btspdat[btspdat$group==1,c("left", "right")], maxiter = iterations)
         btrt_pf <- data.frame(L = NPMLE.trt$intmap[1,],
                              R = NPMLE.trt$intmap[2,],
                              drop = NPMLE.trt$pf)
         btrt_pf <- btrt_pf[btrt_pf$drop!=0,]
         btrt_pf$cumdrop = cumsum(btrt_pf$drop)
         res_new <- getIGroc(npmle_0 = bcontrol_pf, npmle_1 = btrt_pf,
                             xlab, ylab, main, cex.axis, cex.lab, lwd, silenceplot=TRUE)
         res_temp <- getINTERPcurve(bcontrol_pf, btrt_pf, res_new, 1000)

         #bstAUC[b] <- temp$AUC
         forplot <- res_temp$res_temp[,1:2]
         colnames(forplot) = colnames(res_new)[1:2]

         #lines(forplot, col = "bisque")
         temp <- merge(ref[,1:2], forplot, by.x = "ref_u", by.y = "u", all = TRUE)
         if (length(table(is.na(temp$ref_Ru)))>1){
           temp <- temp[-which(is.na(temp$ref_Ru)),]
         }
         temp <- temp %>%
           arrange(desc(ref_u), desc(ref_Ru)) %>%   # sort ref_u desc, then ref_Ru desc within each ref_u
           group_by(ref_u) %>% mutate(R_u = sort(R_u, decreasing = TRUE, na.last = TRUE)) %>%   # re-sort R_u descending independently within each ref_u group
           ungroup()
         temp$R_u <- ifelse(temp$R_u<temp$ref_Ru, temp$ref_Ru, temp$R_u)
         temp <- temp %>%
           arrange(desc(ref_u), desc(ref_Ru)) %>%   # sort ref_u desc, then ref_Ru desc within each ref_u
           group_by(ref_u) %>% mutate(R_u = sort(R_u, decreasing = TRUE, na.last = TRUE)) %>%   # re-sort R_u descending independently within each ref_u group
           ungroup()
         l <- length(which(is.na(temp$R_u)))
         while (l > 0) {
           inds <- which(is.na(temp$R_u))
           temp$R_u[inds] <- temp$R_u[inds - 1]
           l <- length(which(is.na(temp$R_u)))
         }
         temp <- temp[!duplicated(temp),]
         check = ifelse(length(table(temp$R_u==1))==1|table(temp$R_u==1)[2]>3, check+1, -50)
       }

       bst2SS[[b]] <- temp
       # num4prop[b] <- length(which(forplot[,1]>forplot[,2]))
       # denom4prop[b] <- nrow(forplot)
        print(b)
     }

     #prop <- sum(num4prop)/sum(denom4prop)
     names(bst2SS) = paste("s", 1:B, sep="")
     se = CIlow = CIup = rep(NA, length(xunik))

     for (xucount in 1:nrow(ref)){
         temp <- unlist(lapply(bst2SS, function(x) x$R_u[x$ref_u==ref[xucount, 1] & x$ref_Ru==ref[xucount, 2]]))
         se[xucount] <- sd(temp)
         CIlow[xucount] <- min(sort(temp)[(conf.lev)*length(temp)], ifelse(xucount>1, CIlow[xucount-1], CIlow[xucount]))
         CIup[xucount] <- min(sort(temp)[(1-conf.lev)*length(temp)], ifelse(xucount>1, CIup[xucount-1], CIup[xucount]))
     }

     boot <- data.frame(u=xunik, se_Ru = se, CIlow_Ru = CIlow, CIup_Ru = CIup)
     for (i in 2:nrow(boot)){
            boot$CIlow_Ru[i] <- ifelse(boot$CIlow_Ru[i]>boot$CIlow_Ru[i-1], boot$CIlow_Ru[i-1], boot$CIlow_Ru[i])
         boot$CIup_Ru[i] <- ifelse(boot$CIup_Ru[i]>boot$CIup_Ru[i-1], boot$CIup_Ru[i-1], boot$CIup_Ru[i])
         }
     R_u <- merge(ref[,1:2], boot, by.x = "ref_u", by.y = "u", all.x = TRUE)
     colnames(R_u)[1:2] = c("u", "Ru")
     R_u <- R_u[order(R_u$u, decreasing = TRUE), ]
     R_u$CIup_Ru <- ifelse(R_u$CIup_Ru<R_u$Ru, R_u$Ru, R_u$CIup_Ru)
     R_u$CIlow_Ru <- ifelse(R_u$CIlow_Ru>R_u$Ru, R_u$Ru, R_u$CIlow_Ru)
     lines(R_u$u, R_u$CIlow_Ru, lty = 2, col = "black", lwd = lwd)
     lines(R_u$u, R_u$CIup_Ru, lty = 2, col = "black", lwd = lwd)

     AUC <- c(AUC = res$AUC,
            CIlow = sort(bstAUC)[(conf.lev)*length(bstAUC)],
            CIup = sort(bstAUC)[(1-conf.lev)*length(bstAUC)])

     return(R_u)
 }
