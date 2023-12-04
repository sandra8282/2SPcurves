solvenaC <- function(indi, temp2){
    if (temp2$tienext[indi[1]-1] == 1){
        xcoords <- temp2$ref_u[c(indi[length(indi)]+1, indi[1]-1)]
        ycoords <- temp2$y[c(indi[length(indi)]+1, indi[1]-1)]
        slope <- diff(ycoords)/diff(xcoords)
        intercept = ycoords[2] - slope*xcoords[2]
        temp2$y[indi] = slope*temp2$ref_u[indi]+intercept
        temp2$tienext[indi] = 1
    } else {
        temp2$y[indi] = temp2$y[indi[1]-1]
        temp2$tienext[indi] = 0
        }
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
#' @importFrom zoo rollmean
#' @importFrom Bolstad2 sintegral
#' @keywords internal
#' @noRd

btsp2SCI <- function(res, maindat, ref_skm, nrisktypes, B, level, xlab, ylab, rlabels, main,
                     cex.axis = 1, cex.lab = 1, lty = 1, lwd = 1, CIabcd = TRUE,
                     legend.inset=0.01, legend.cex=1.25, silent = TRUE) {

    conf.lev = (1-level)/2
    maindat = data.frame(maindat)
    g0dat <- data.frame(maindat[maindat$group==1,]);
    g1dat <- data.frame(maindat[maindat$group==0,]);

    ref <- res
    changecolnam <- function(list_element){
        colnames(list_element) <- c("ref_u", "ref_Cu", "ref_tienext")
        return(list_element)}
    ref <- lapply(ref, function(x) changecolnam(x))

    toreturn = ref
    bstABCD = matrix(rep(NA, B*2), ncol=2)
    bstedat = matrix(rep(NA, B*4), ncol=4)
    colnames(bstedat) = c("startX", "startY", "endX", "endY")
    bst2SS = vector(mode = "list", length = B);

    for (b in 1:B) {
        set.seed(b+76545349)
        btspdat0 <- data.frame(rbind(g0dat[sample(1:nrow(g0dat), replace = TRUE), ],
                                     g1dat[sample(1:nrow(g1dat), replace = TRUE), ]))
        time=btspdat0$time; group=btspdat0$group; event = btspdat0$event;

            ##### get cuminc
            if (0 %in% event){
                #censoring
                btspdat <- data.frame(id=1:length(time), time=time, group=group)
                btspdat$eventf <- factor(ifelse(event==0, "censored",
                                                ifelse(event==1, rlabels[1], rlabels[2])))
            } else {
                #no censoring
                btspdat <- data.frame(id = 1:(length(time)+2),
                                      time = c(time, max(time)*2, max(time)*2),
                                      group = c(group, 0, 1))
                newevent = c(event, 0,0)
                btspdat$eventf <- factor(ifelse(newevent==0, "censored",
                                                ifelse(newevent==1, rlabels[1], rlabels[2])))
            }
            fit <-  survfit(Surv(btspdat$time, btspdat$eventf) ~ btspdat$group)
            sfit <- summary(fit)
            cis <- sfit$pstate
            skm <- data.frame(time = sfit$time,
                              group = as.numeric(sfit$strata)-1,
                              ci = cis[,-1])
            bstedat[b,] = c(max(skm$ci.1[skm$group==0]), max(skm$ci.1[skm$group==1]),
                            max(skm$ci.2[skm$group==0]), max(skm$ci.2[skm$group==1]))

            ### get C_u's, ABCD's and map C_u's to ref
            list_4plot = bstres = list(); abcds_b = rep(NA, nrisktypes);
            for (skmi in (1:nrisktypes)){
                skmires = cbind(time = sfit$time, ci = cis[,skmi+1], group=as.numeric(sfit$strata)-1)
                skmires = skmires[!duplicated(skmires[,2:3]),]
                skmires = skmires[order(skmires[,1], skmires[,3]),]
                ties_check <- which(duplicated(skmires[,1]))
                if (length(ties_check) > 1) {
                    ties_times = skmires[duplicated(skmires[,1]),1]
                    ties_ind <- rep(0, nrow(skmires))
                    ties_ind[which(skmires[,1] %in% ties_times)]=1
                } else {ties_ind <- rep(0, nrow(skmires))}
                skmires = cbind(skmires, ties_ind)
                temp = get4plotCumInc(skmires)
                colnames(temp) = c("x", "y", "tienext")
                id <- 1:nrow(temp)
                if (skmi==1){
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- pracma::trapz(temp[,1],temp[,2]) - pracma::trapz(temp[,1],temp[,1])
                    list_4plot[[skmi]] = temp
                    options(warn = defaultW)
                } else {
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- pracma::trapz(temp[,1],temp[,1])-pracma::trapz(temp[,1],temp[,2])
                    list_4plot[[skmi]] = temp
                    options(warn = defaultW)
                }
                abcds_b[skmi] = abcdi
                if (CIabcd==FALSE){
                    new <- data.frame(ref[[skmi]])
                    new2 <- data.frame(temp)
                    new2$tienext[nrow(new2)]=0
                    unix = unique(new2$x)
                    new2_orig = new2
                    new2 = NULL
                    for (ux in 1:length(unix)){
                      new2sub = subset(new2_orig, new2_orig$x==unix[ux])
                      newy = max(new2sub$y)
                      new2 = rbind(new2,
                                   c(x=unix[ux], y=newy, tienext = new2sub[nrow(new2sub),3]))
                    }
                    new2 = data.frame(new2)
                    tempnew <- merge(new[-1,], new2, by.x="ref_u", by.y = "x", all=TRUE)
                    tempnew <- rbind(rep(0, ncol(tempnew)), tempnew)
                    temp_orig = tempnew
                    temp=NULL
                    cond=as.numeric(length(na.omit(tempnew[,4]))!=nrow(tempnew))
                    while (cond==1){
                        inds = which(is.na(tempnew$y))
                        inds = split(inds, cumsum(c(1, diff(inds) != 1)))
                        all <- lapply(inds, function(x) solvenaC(x, temp2 = tempnew))
                        for(i in 1:length(inds)) {
                            newy = all[[i]]$y[inds[[i]]]
                            tempnew$y[inds[[i]]] = newy
                        }
                        cond <- as.numeric(length(na.omit(tempnew[,4]))!=nrow(tempnew))
                    }
                    temp = tempnew[!is.na(tempnew$ref_Cu),]
                    bstres[[skmi]]=temp
                }
            }
            bstABCD[b,] = abcds_b
            if (CIabcd==FALSE){bst2SS[[b]] = bstres}
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

    if (CIabcd==FALSE){
        names(bst2SS) = paste("CUMINC", 1:B, sep="")
        xunik = se = CIlow = CIup = boot = C_u = list()

        for (uniki in 1:nrisktypes){
            xunik[[uniki]] <- unique(ref[[uniki]][,1])
            se[[uniki]] = CIlow[[uniki]] = CIup[[uniki]] = rep(NA, length(xunik[[uniki]]))
            for (xucount in 1:length(xunik[[uniki]])){
                lc_temp <- lapply(bst2SS, function(x) x[[uniki]])
                temp_all <- lapply(lc_temp, function(x) x[abs(x$ref_u-xunik[[uniki]][xucount])<0.00000001,])
                temp_all <- lapply(temp_all, function(x) x[!(x$ref_Cu==0),])
                temp2 <- unlist(lapply(temp_all, function(x) x$y))
                se[[uniki]][xucount] <- sd(temp2)
                CIlow[[uniki]][xucount] <- sort(temp2)[(conf.lev)*length(temp2)]
                #CIlow[[uniki]][xucount] <- temp2 - qnorm(1-conf.lev)*sd(temp2)
                CIup[[uniki]][xucount] <- sort(temp2)[(1-conf.lev)*length(temp2)]
                #CIlow[[uniki]][xucount] <- temp2 + qnorm(1-conf.lev)*sd(temp2)
            }
            boot[[uniki]] <- data.frame(u=xunik[[uniki]], se_Cu = se[[uniki]],
                                        CIlow_Cu = CIlow[[uniki]] , CIup_Cu = CIup[[uniki]] )
            C_u[[uniki]] <- merge(ref[[uniki]][-1,], boot[[uniki]],
                                  by.x = "ref_u", by.y = "u", all.x = TRUE)
            colnames(C_u[[uniki]])[1:2] = c("u", "Cu")
            C_u[[uniki]] <- C_u[[uniki]][order(C_u[[uniki]]$u, decreasing = FALSE), ]
            C_u[[uniki]] <- rbind(rep(0, ncol(C_u[[uniki]])), C_u[[uniki]])
        }
        finalres <- list(C_u=C_u, abcd=abcd)
        if (silent==FALSE){
            ptypes=c(NA, 1)
            darkcolopt =c("black", "grey29")
            lwdopt = c(lwd, 1)
            plot2SCI(C_u, xlab=xlab, ylab=ylab, main=main, rlabels=rlabels,
                     cex.axis = cex.axis, cex.lab = 1.5, lwd = 1.5, silent=FALSE, lty = lty,
                     legend.inset=legend.inset, legend.cex=legend.cex, btspind=TRUE)
        }
    } else {finalres <- list(abcd=abcd)}
     #, abcd = abcd

    return(finalres)
}

