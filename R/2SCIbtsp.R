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
#' @import dplyr
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
    gdat <- data.table(maindat);
    n=nrow(gdat)
    rownames(gdat) = gdat$id

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
    bst2SS = vector(mode = "list", length = B);

    for (b in 1:B) {
            set.seed(b+7689)
            gids_b <- sample(gdat$id, n, replace = TRUE)
            btspdat <- data.table(gdat[gids_b, ])

            ##### get cuminc
            if (0 %in% btspdat$event){
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

            ###check ties
            ties_check <- unique(table(skm[,1]))
            if (length(ties_check) > 1) {
                ties_times = skm[duplicated(skm[,1]),1]
                ties_ind <- rep(0, nrow(skm))
                ties_ind[which(skm[,1] %in% ties_times)]=1
            } else {ties_ind <- rep(0, nrow(skm))}

            #### get skm info
            skm = cbind(skm, ties_ind)
            skm = skm[order(skm$time, skm$group),]
            list_4plot = bstres = list(); abcds_b = rep(NA, nrisktypes);

            ### get C_u's, ABCD's and map C_u's to ref
            for (skmi in (1:nrisktypes)){
                skmires = skm[, c(1, 2+skmi, 2)]
                ties_check <- duplicated(skmires[,1])
                if (length(ties_check) > 1) {
                    ties_times = skmires[duplicated(skmires[,1]),1]
                    ties_ind <- rep(0, nrow(skmires))
                    ties_ind[which(skmires[,1] %in% ties_times)]=1
                } else {ties_ind <- rep(0, nrow(skmires))}
                skmires = cbind(skmires, ties_ind)
                temp = get4plotCumInc(skmires)
                temp = temp[!duplicated(temp),]
                id <- 1:nrow(temp)
                if (skmi==1){
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- sintegral(temp[,1],temp[,2]-temp[,1])$int
                    list_4plot[[skmi]] = temp
                    options(warn = defaultW)
                } else {
                    defaultW <- getOption("warn")
                    options(warn = -1)
                    abcdi <- sintegral(temp[,1],temp[,1]-temp[,2])$int
                    list_4plot[[skmi]] = temp
                    options(warn = defaultW)
                }
                abcds_b[skmi] = abcdi
                new <- data.frame(ref_skm[,c(1:2, 2+skmi)])
                new <- new[-(new$ref_u==0&new$ref_Cu==0|duplicated(new)),]
                new2 <- data.frame(skm[, c(1, 2+skmi, 2)])
                new2 <- new2[-(new2$x==0&new2$y==0|duplicated(new2)),]
                temp <- merge(new, new2, by.x = "ref_u", by.y = "x", all = TRUE)
                uniquex = unique(temp$ref_u)
                if (TRUE %in% c(FALSE, diff(uniquex) < 1*10^-15)){
                  uniquex = uniquex[-which(c(FALSE, diff(uniquex) < 1*10^-15)==TRUE)]
                }
                uniquex = round(uniquex, 12)
                temp_orig = temp
                temp = NULL
                newy=NULL
                for (tn in 1:length(uniquex)){
                    tempsub <- subset(temp_orig, abs(uniquex[tn]-ref_u)<1*10^-10)
                    tempsub$newy = NA
                    ys <- unique(tempsub$y)
                    if (length(ys)==1|unique(is.na(tempsub$y))|unique(is.na(tempsub$ref_Cu))){
                        tempsub$newy=tempsub$y
                        } else {
                            inds <- apply(tempsub, 1, function(x) which(abs(x[2]-ys)==min(abs(x[2]-ys))))
                            tempsub$newy =ys[inds]
                        }
                    temp=rbind(temp, tempsub)
                }
                temp$y = temp$newy
                temp <- temp[,-6]
                temp$inds <- c(FALSE, diff(temp$ref_u)==0)


                if (length(which(temp$y==0))>0){
                   all0 = 1:max(which(temp$y==0))
                   temp$y[all0] = 0
                   temp$tienext[all0] = 0
                } else {if (is.na(temp$y[1])){temp$y[1]=0; temp$tienext[1]=0}}
                temp = temp[order(temp$ref_u, temp$y),]
                inds1 <- which(is.na(temp$y)); l <- length(inds1)
                if (l>0){
                    inds1 <- split(inds1, cumsum(c(1, diff(inds1) != 1)))
                    temp$tienext[inds1[[length(inds1)]][1]-1]=0
                    }
                while (l > 0) {
                    inds <- which(is.na(temp$y))
                    inds <- split(inds, cumsum(c(1, diff(inds) != 1)))
                    for (indsi in 1:length(inds)){
                        temp <- solvenaC(indi = inds[[indsi]], temp2 = temp)
                    }
                    l <- length(which(is.na(temp$y)))
                }
                if (length(which(is.na(temp$ref_Cu)))>0){
                    temp <- temp[-which(is.na(temp$ref_Cu)),]
                }
                bstres[[skmi]]=temp
            }
            bstABCD[b,] = abcds_b
            bst2SS[[b]] <- bstres
    }

    names(bst2SS) = paste("CUMINC", 1:B, sep="")
    xunik = se = CIlow = CIup = boot = C_u = list()

    for (uniki in 1:nrisktypes){
        xunik[[uniki]] <- unique(ref[[uniki]][,1])
        se[[uniki]] = CIlow[[uniki]] = CIup[[uniki]] = rep(NA, length(xunik[[uniki]]))
        for (xucount in 1:length(xunik[[uniki]])){
            #### get the 1st of all bst2SS
            lc_temp <- lapply(bst2SS, function(x) x[[1]])
            temp_all <- lapply(lc_temp, function(x) x[abs(x$ref_u-xunik[[uniki]][xucount])<0.00000001,])
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

