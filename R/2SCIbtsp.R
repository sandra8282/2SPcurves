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
#' @importFrom zoo rollmean
#' @importFrom Bolstad2 sintegral
#' @keywords internal
#' @noRd

btsp2SCI <- function(res, maindat, ref_skm, nrisktypes, B, level, xlab, ylab, rlabels, main,
                     cex.axis = cex.axis, cex.lab = cex.lab, lty = 1, lwd = lwd,
                     legend.inset=0.01, legend.cex=1.25, silent = TRUE, cens=TRUE) {
    ptypes=c(NA, 1)
    darkcolopt =c("black", "grey29")
    lwdopt = c(lwd, 1)

    enames = c("censored", rlabels)
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

    #ref_times = data.frame(time = reft)
    # for (refi in 1:nrisktypes){
    #     lines(ref[[refi]][,1],ref[[refi]][,2], lty = 2, lwd = lwd) #darkcolopt[[refi]]
    # }

    toreturn = list()
    bstCindex = matrix(rep(NA, B*2), ncol=2)
    bstABCD = matrix(rep(NA, B*2), ncol=2)
    bstedat = matrix(rep(NA, B*4), ncol=4)
    colnames(bstedat) = c("startX", "startY", "endX", "endY")

    for (i in 1:nrisktypes){
        bst2SS = vector(mode = "list", length = B);
        for (b in 1:B) {
            set.seed(b+7679)
            g0ids_b <- sample(g0ids, n0, replace = TRUE)
            g1ids_b <- sample(g1ids, n1, replace = TRUE)
            btspdat <- data.table(rbind(g0dat[g0ids_b, ], g1dat[g1ids_b, ]))

            ##### get cuminc
            if (cens==TRUE){
                btspdat$eventf <- factor(btspdat$event, 0:nrisktypes, labels = enames)
                #get C(u)
            } else {
                #no censoring
                synthdat <- data.frame(id = nrow(maindat)+1:2,
                                       time = rep(max(btspdat$time)*10,2),
                                       event = c(0,0),
                                       group = c(0,1))
                btspdat <- rbind(btspdat, synthdat)
                btspdat$eventf <- factor(btspdat$event, 0:nrisktypes, labels = enames)
            }

            fit <-  survfit(Surv(btspdat$time, btspdat$eventf) ~ btspdat$group)
            sfit <- summary(fit)
            skm <- data.frame(time = sfit$time, group = as.numeric(sfit$strata)-1, ci = sfit$pstate[,2:3])
            bstedat[b,] = c(max(skm$ci.1[skm$group==0]), max(skm$ci.1[skm$group==1]),
                            max(skm$ci.2[skm$group==0]), max(skm$ci.2[skm$group==1]))

            ties_check <- unique(table(skm[,1]))
            if (length(ties_check) > 1) {
                ties_times = skm[duplicated(skm[,1]),1]
                ties_ind <- rep(0, nrow(skm))
                ties_ind[which(skm[,1] %in% ties_times)]=1
            } else {ties_ind <- rep(0, nrow(skm))}
            skm = cbind(skm, ties_ind)
            skm = skm[order(skm$time, skm$group),]
            list_4plot = bstres = list(); abcds = NULL;
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
                temp = get4plotCumInc(skmires)
                temp = temp[!duplicated(temp),]
                outs <- which(duplicated(temp[,1:2])&temp[,2]==0)
                if(length(outs)>=1) {temp = temp[-(duplicated(temp[,1:2])&temp[,2]==0),]}
                id <- 1:nrow(temp)

                if (skmi==1){
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- sintegral(temp[,1],temp[,2])$int - sintegral(temp[,1],temp[,1])$int
                    list_4plot[[skmi]] = temp
                    abcds[skmi] = abcdi
                    #kstest = ks.test(x=temp[,2], y=temp[,1])
                    #tests[skmi] = c(pval = round(kstest$p.value, 6))
                    options(warn = defaultW)
                } else {
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- sintegral(temp[,1],temp[,1])$int - sintegral(temp[,1],temp[,2])$int
                    list_4plot[[skmi]] = temp
                    abcds[skmi] = abcdi
                    #kstest = ks.test(x=temp[,2], y=temp[,1])
                    #tests[skmi] = c(pval = round(kstest$p.value, 6))
                    options(warn = defaultW)
                }

            }
            bstABCD[b,] = abcds
            for (ploti in (1:nrisktypes)){
                new <- data.frame(ref[[ploti]])
                outs <- unique(c(which(new$ref_u==0&new$ref_Cu==0), which(duplicated(new))))
                new <- new[-outs,]
                new2 <- data.frame(list_4plot[[ploti]])
                outs <- unique(c(which(new2$x==0&new2$y==0), which(duplicated(new2))))
                new2 <- new2[-outs,]
                temp <- merge(new, new2, by.x = "ref_u", by.y = "x", all = TRUE)
                temp <- temp[order(temp$ref_u, temp$ref_Cu, decreasing = FALSE), ]
                fix <- which(temp$y<temp$ref_Cu&(temp$y==0))
                temp$y[fix] = temp$y[fix-1]
                if (is.na(temp$y[1])){temp$y[1]=0; temp$tienext[1]=0}

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
                if (length(which(is.na(temp$ref_Cu)))>0){
                    temp <- temp[-which(is.na(temp$ref_Cu)),]
                }
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
            #CIlow[[uniki]][xucount] <- temp2 - qnorm(1-conf.lev)*sd(temp2)
            CIup[[uniki]][xucount] <- sort(temp2)[(1-conf.lev)*length(temp2)]
            #CIlow[[uniki]][xucount] <- temp2 + qnorm(1-conf.lev)*sd(temp2)
        }
        boot[[uniki]] <- data.frame(u=xunik[[uniki]], se_Cu = se[[uniki]],
                                    CIlow_Cu = CIlow[[uniki]] , CIup_Cu = CIup[[uniki]] )
        C_u[[uniki]] <- merge(ref[[uniki]], boot[[uniki]], by.x = "ref_u", by.y = "u", all.x = TRUE)
        colnames(C_u[[uniki]])[1:2] = c("u", "Cu")
        C_u[[uniki]] <- C_u[[uniki]][order(C_u[[uniki]]$u, decreasing = FALSE), ]
    }

    bstABCD = cbind(bstABCD, rowSums(bstABCD))
    abcd = t(apply(bstABCD, 2, function(x) c(estimate = mean(x), se = sd(x),
                                             CIlower = sort(x)[(conf.lev)*length(x)],#CIlower =  mean(x) - qnorm(1-conf.lev)*sd(x),
                                             CIupper = sort(x)[(1-conf.lev)*length(x)]))) #CIupper = mean(x) + qnorm(1-conf.lev)*sd(x)
    bstedat = data.frame(bstedat)
    bstedat$endX = 1-bstedat$endX
    bstedat$endY = 1-bstedat$endY
    bstedat$area = ((bstedat$startY-bstedat$startX)+(bstedat$endY-bstedat$endX))*(bstedat$endX-bstedat$startX)/2
    eareas = bstABCD[,3] + bstedat$area
    eABCD = c(mean(eareas), sd(eareas),
              sort(eareas)[(conf.lev)*length(eareas)], #mean(eareas) - qnorm(1-conf.lev)*sd(eareas), #
              sort(eareas)[(1-conf.lev)*length(eareas)]) #mean(eareas) + qnorm(1-conf.lev)*sd(eareas)) #
    abcd = rbind(abcd, eABCD)
    rownames(abcd) = c(rlabels, "Overall", "Extrapolated")

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
        }}

    finalres <- list(C_u=C_u, abcd=abcd) #, abcd = abcd

    if (silent==FALSE){
        plot2SCI(C_u, xlab=xlab, ylab=ylab, main=main, rlabels=rlabels,
                 cex.axis = cex.axis, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = lty,
                 legend.inset=legend.inset, legend.cex=legend.cex, btspind=TRUE)
    }

    return(finalres)
}

